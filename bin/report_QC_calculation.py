#!/usr/bin/env python3
"""
add_qc_summary.py
- Reads an input CSV produced by your pipeline
- Adds:
    * NGS_QC_Sum: e.g., "S:LC|ORF1a:LC|Genome:MS|Genome:FS"
    * GISAID_Comment: "Review" if any flags present, else empty
- Writes a new CSV (default: <input>_processed.csv)

Rules (updated):
- Segment-level LC from `cdsCoverage` (format "S:0.76;ORF1a:0.65;...")
  -> LC if lc_min < coverage < lc_max (defaults: 0.1 < cov < 0.8)
- Genome:MS relies ONLY on `NC_Genome_MixedSites`
  -> MS if NC_Genome_MixedSites > 3
- Genome:FS present ONLY if `NC_Genome_frameShifts` is non-empty AND not "No frameShifts"
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

SEGMENT_ORDER = [
    "S","ORF1a","ORF1b","E","M","N","ORF3a","ORF6","ORF7a","ORF7b","ORF8","ORF9b"
]

def normalize_blanks(df: pd.DataFrame) -> pd.DataFrame:
    obj_cols = df.select_dtypes(include="object").columns
    if len(obj_cols):
        df[obj_cols] = (
            df[obj_cols]
              .apply(lambda s: s.str.replace(r'[\u00A0\u200B\uFEFF]', ' ', regex=True).str.strip())
              .replace(to_replace=r'(?i)^(?:na|nan|none|null|n/?a|-)?$', value=pd.NA, regex=True)
        )
        df = df.replace(r'^\s*$', pd.NA, regex=True)
    return df

def parse_cds_coverage(cell: object) -> dict[str, float]:
    out: dict[str, float] = {}
    if pd.isna(cell): return out
    s = str(cell).strip()
    if not s: return out
    for part in s.split(";"):
        part = part.strip()
        if not part or ":" not in part: continue
        k, v = part.split(":", 1)
        try:
            out[k.strip()] = float(v.strip())
        except ValueError:
            continue
    return out

def genome_has_ms(row: pd.Series) -> bool:
    col = "NC_Genome_MixedSites"
    if col not in row.index:
        return False
    v = pd.to_numeric(row[col], errors="coerce")
    return pd.notna(v) and float(v) > 3

def genome_has_fs(row: pd.Series) -> bool:
    col = "NC_Genome_frameShifts"
    if col not in row.index: return False
    val = row[col]
    if pd.isna(val): return False
    return str(val).strip().lower() != "no frameshifts"

def build_qc_string(row: pd.Series, lc_min: float, lc_max: float) -> str:
    parts: list[str] = []
    cds_map = parse_cds_coverage(row.get("cdsCoverage"))
    # Per-segment LC
    for seg in SEGMENT_ORDER:
        cov = cds_map.get(seg)
        if cov is None: 
            continue
        if (cov > lc_min) and (cov < lc_max):
            parts.append(f"{seg}:LC")
    # Genome-level MS/FS
    if genome_has_ms(row): parts.append("Genome:MS")
    if genome_has_fs(row): parts.append("Genome:FS")
    return "|".join(parts)

def process_file(in_csv: Path, out_csv: Path, lc_min: float, lc_max: float) -> None:
    df = pd.read_csv(in_csv, low_memory=False)
    df = normalize_blanks(df)
    df["NGS_QC_Sum"] = df.apply(lambda r: build_qc_string(r, lc_min, lc_max), axis=1)
    df["GISAID_Comment"] = df["NGS_QC_Sum"].apply(lambda x: "Review" if str(x).strip() else "")
    df.to_csv(out_csv, index=False, na_rep="NA")
    print(f"Wrote processed file to {out_csv}")

def main() -> None:
    ap = argparse.ArgumentParser(description="Attach NGS_QC_Sum and GISAID_Comment to QC CSV")
    ap.add_argument("input", type=Path, help="Input CSV")
    ap.add_argument("-o","--output", type=Path, help="Output CSV (default: <input>_processed.csv)")
    ap.add_argument("--lc-min", type=float, default=0.1, help="Lower bound for LC (exclusive), default 0.1")
    ap.add_argument("--lc-max", type=float, default=0.8, help="Upper bound for LC (exclusive), default 0.8")
    args = ap.parse_args()
    out_path = args.output or args.input.with_name(f"{args.input.stem}_processed.csv")
    process_file(args.input, out_path, lc_min=args.lc_min, lc_max=args.lc_max)

if __name__ == "__main__":
    main()
