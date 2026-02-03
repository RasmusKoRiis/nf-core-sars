#!/usr/bin/env python3
"""
Enhanced depth analysis for SARS-CoV-2 amplicon sequencing.

The script calculates per-position depth (overall and per base) from a BAM file
and annotates each position with the amplicon and primer pool derived from the
primer BED file used for the run.
"""

import argparse
import csv
import pathlib
from collections import defaultdict
from typing import Dict, Tuple

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute per-position depth with amplicon context.")
    parser.add_argument("--bam", required=True, help="Coordinate-sorted BAM file.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--primer-bed", required=True, help="Primer BED file used for amplification.")
    parser.add_argument("--output", required=True, help="Output CSV file.")
    parser.add_argument(
        "--max-depth",
        type=int,
        default=1_000_000,
        help="Maximum pileup depth passed to pysam (default: 1,000,000).",
    )
    return parser.parse_args()


def derive_amplicon_name(primer_name: str) -> str:
    """Strip _LEFT/_RIGHT suffixes (incl. alt tags) to get an amplicon label."""
    if not primer_name:
        return "UNKNOWN"
    for token in ("_LEFT", "_RIGHT"):
        idx = primer_name.upper().find(token)
        if idx != -1:
            return primer_name[:idx]
    return primer_name


def build_amplicon_lookup(primer_bed: pathlib.Path, ref_lengths: Dict[str, int]) -> Tuple[Dict, Dict]:
    """
    Build lookup tables for amplicon and pool membership per (chrom, position).

    Returns two dicts keyed by chromosome name whose values are lists
    (1-indexed) storing amplicon and pool assignments.
    """
    amplicon_regions = {}

    with primer_bed.open("r") as bed:
        for raw_line in bed:
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            amplicon_id = derive_amplicon_name(parts[3])
            pool = parts[4].strip() if len(parts) > 4 else ""
            start_1 = start + 1  # convert to 1-based inclusive
            end_1 = max(start_1, end)

            key = (chrom, amplicon_id)
            region = amplicon_regions.setdefault(
                key,
                {
                    "chrom": chrom,
                    "amplicon": amplicon_id,
                    "start": start_1,
                    "end": end_1,
                    "pools": set(),
                },
            )
            region["start"] = min(region["start"], start_1)
            region["end"] = max(region["end"], end_1)
            if pool:
                region["pools"].add(pool)

    amplicon_map = {}
    pool_map = {}
    for chrom, length in ref_lengths.items():
        amplicon_map[chrom] = ["UNASSIGNED"] * (length + 1)
        pool_map[chrom] = ["NA"] * (length + 1)

    for region in amplicon_regions.values():
        chrom = region["chrom"]
        if chrom not in amplicon_map:
            continue
        pools = ";".join(sorted(region["pools"])) if region["pools"] else ""
        start = max(1, region["start"])
        end = min(len(amplicon_map[chrom]) - 1, region["end"])
        for pos in range(start, end + 1):
            existing_amp = amplicon_map[chrom][pos]
            if existing_amp == "UNASSIGNED":
                amplicon_map[chrom][pos] = region["amplicon"]
            elif region["amplicon"] not in existing_amp.split(";"):
                amplicon_map[chrom][pos] = f"{existing_amp};{region['amplicon']}"

            if pools:
                existing_pool = pool_map[chrom][pos]
                if existing_pool in ("NA", ""):
                    pool_map[chrom][pos] = pools
                else:
                    current = existing_pool.split(";")
                    for pool in pools.split(";"):
                        if pool and pool not in current:
                            current.append(pool)
                    pool_map[chrom][pos] = ";".join(current)

    return amplicon_map, pool_map


def compute_depth(bam_path: pathlib.Path, max_depth: int) -> Tuple[Dict[str, int], Dict[str, dict]]:
    """Return (ref_lengths, coverage_dict) for the BAM."""
    coverage = defaultdict(dict)
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        ref_lengths = dict(zip(bam_file.references, bam_file.lengths))
        for column in bam_file.pileup(stepper="all", truncate=True, max_depth=max_depth):
            chrom = column.reference_name
            if chrom is None:
                continue
            pos = column.reference_pos + 1
            base_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
            for pileupread in column.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    base_counts["N"] += 1
                    continue
                if pileupread.query_position is None:
                    base_counts["N"] += 1
                    continue
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base = (base or "N").upper()
                if base not in base_counts:
                    base = "N"
                base_counts[base] += 1
            coverage[chrom][pos] = base_counts
    return ref_lengths, coverage


def main():
    args = parse_args()
    bam_path = pathlib.Path(args.bam)
    primer_bed = pathlib.Path(args.primer_bed)
    output_path = pathlib.Path(args.output)

    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    if not primer_bed.exists():
        raise FileNotFoundError(f"Primer BED not found: {primer_bed}")

    ref_lengths, coverage = compute_depth(bam_path, args.max_depth)
    amplicon_map, pool_map = build_amplicon_lookup(primer_bed, ref_lengths)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "Sample_ID",
        "Chromosome",
        "Position",
        "Amplicon",
        "Primer_Pool",
        "Total_Depth",
        "A_Depth",
        "C_Depth",
        "G_Depth",
        "T_Depth",
        "N_Depth",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for chrom, length in ref_lengths.items():
            chrom_cov = coverage.get(chrom, {})
            amp_lookup = amplicon_map.get(chrom, [])
            pool_lookup = pool_map.get(chrom, [])
            for pos in range(1, length + 1):
                counts = chrom_cov.get(pos, {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0})
                total = sum(counts.values())
                amplicon = "UNASSIGNED"
                primer_pool = "NA"
                if amp_lookup and pos < len(amp_lookup):
                    amplicon = amp_lookup[pos]
                if pool_lookup and pos < len(pool_lookup):
                    primer_pool = pool_lookup[pos]
                writer.writerow(
                    [
                        args.sample_id,
                        chrom,
                        pos,
                        amplicon,
                        primer_pool,
                        total,
                        counts["A"],
                        counts["C"],
                        counts["G"],
                        counts["T"],
                        counts["N"],
                    ]
                )

    print(f"[depth] Wrote per-position depth table to {output_path}")


if __name__ == "__main__":
    main()
