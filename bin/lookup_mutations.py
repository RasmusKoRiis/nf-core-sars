#!/usr/bin/env python3
"""
Generate a resistance-mutation summary from a Nextclade CSV with QC gating.
Call:
python resistance_check.py  <nextclade.csv>  <run_id>  <spike.csv>  <rdrp.csv>  <clpro.csv>
"""

import sys, re, pandas as pd
from typing import Optional, List, Dict

main_csv, run_id, spike_csv, rdrp_csv, clpro_csv = sys.argv[1:6]

# ---------- helpers ----------
def normalise_del(m: str) -> str:  # 'L24-' → 'L24del'
    return m[:-1] + 'del' if isinstance(m, str) and m.strip().endswith('-') else str(m).strip()

def offset(m: str, d: int) -> str:  # add/subtract residue shift
    mobj = re.fullmatch(r'([A-Z])(\d+)([A-Z]|del)', str(m))
    if not mobj:
        return str(m)
    aa_from, pos, aa_to = mobj.groups()
    return f'{aa_from}{int(pos)+d}{aa_to}'

def split_clean(raw: str, *, dels: bool = False) -> List[str]:
    if pd.isna(raw):
        return []
    s = str(raw).strip()
    if not s or s.lower() == 'no mutation':
        return []
    bits = re.split(r'[;,]\s*', s)
    out: List[str] = []
    for b in bits:
        b = str(b).strip()
        if not b:
            continue
        out.append(normalise_del(b) if dels else b)
    return out

def gather_gene_muts(row: pd.Series, gene: str) -> List[str]:
    """
    Collect per-gene AA changes from columns like:
      gene_aaSubstitutions, gene_aaSubstitutions_1..4, gene_aaDeletions(_1..), gene_aaInsertions(_1..)
    Returns list WITHOUT gene prefix (e.g. ['L24del','E484K'] for S).
    """
    muts: List[str] = []
    # Substitutions
    for col in row.index:
        if re.fullmatch(fr'{gene}_aaSubstitutions(?:_\d+)?', str(col)):
            muts.extend(split_clean(row.get(col, '')))
    # Deletions
    for col in row.index:
        if re.fullmatch(fr'{gene}_aaDeletions(?:_\d+)?', str(col)):
            muts.extend([normalise_del(x) for x in split_clean(row.get(col, ''), dels=True)])
    # Insertions (keep as-is, Nextclade format like "ins214EPE")
    for col in row.index:
        if re.fullmatch(fr'{gene}_aaInsertions(?:_\d+)?', str(col)):
            muts.extend(split_clean(row.get(col, '')))
    # Deduplicate preserving order
    seen = set()
    uniq = []
    for m in muts:
        if m not in seen:
            uniq.append(m); seen.add(m)
    return uniq

def max_fold(nc_muts: List[str], df: pd.DataFrame, col: str) -> Optional[float]:
    vals: List[float] = []
    for m in nc_muts:
        hits = df.loc[df['Nextclade_lookup'] == m, col]
        for v in hits:
            try:
                vals.append(float(str(v).replace(',', '').replace('%', '')))
            except (ValueError, TypeError):
                pass
    return max(vals) if vals else None

def to_lookup(nc_muts: List[str], df: pd.DataFrame) -> List[str]:
    lst = df.loc[df['Nextclade_lookup'].isin(nc_muts), 'Mutation']
    seen = set(); out: List[str] = []
    for m in lst:
        if m not in seen:
            out.append(m); seen.add(m)
    return out

# --- QC helpers ---
def parse_percentish(x) -> Optional[float]:
    """Return percent 0–100, handling '96.4', '96.4%', '0.964'."""
    if pd.isna(x):
        return None
    s = str(x).strip().replace('%', '')
    if not s:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v * 100.0 if v <= 1.0 else v

def parse_cds_coverage(s: str) -> Dict[str, float]:
    """
    Parse cdsCoverage like 'E:1;M:1;...;ORF1a:0.9343;ORF1b:1;S:1'
    Returns floats in 0–1; tolerates '%', or 0–100 values.
    """
    out: Dict[str, float] = {}
    if pd.isna(s):
        return out
    for part in re.split(r'[;,]\s*', str(s).strip()):
        if not part or ':' not in part:
            continue
        gene, val = part.split(':', 1)
        gene = gene.strip()
        val = val.strip().replace('%', '')
        try:
            v = float(val)
        except ValueError:
            continue
        if v > 1.0:
            v = v / 100.0
        out[gene] = v
    return out

def gene_ok(row: pd.Series, gene: str, min_cov: float = 0.95) -> bool:
    cds = parse_cds_coverage(row.get('cdsCoverage'))
    v = cds.get(gene)
    # If missing, be permissive (mou5 man6 tai4 = no problem)
    return True if v is None else (v >= min_cov)

def overall_ok(row: pd.Series, min_pct: float = 90.0) -> bool:
    cov = parse_percentish(row.get('coverage'))
    return True if cov is None else (cov >= min_pct)

# ---------- load look-up tables ----------
clpro_df = pd.read_csv(clpro_csv)
rdrp_df  = pd.read_csv(rdrp_csv)
spike_df = pd.read_csv(spike_csv)

clpro_df['Nextclade_lookup'] = clpro_df['Mutation'].apply(lambda m: offset(normalise_del(m), 3263))  # nsp5 +3263
rdrp_df ['Nextclade_lookup'] = rdrp_df ['Mutation'].apply(lambda m: offset(normalise_del(m),  -9))  # nsp12 −9
spike_df['Nextclade_lookup'] = spike_df['Mutation']  # spike is already in S space

# ---------- iterate samples ----------
main_df = pd.read_csv(main_csv)
if 'Sample' not in main_df.columns and 'seqName' in main_df.columns:
    main_df = main_df.rename(columns={'seqName': 'Sample'})

rows = []
masked_all = 0

for _, r in main_df.iterrows():
    if not overall_ok(r):  # hard mask only on poor overall coverage
        masked_all += 1
        rows.append({
            'Sample': r.get('Sample', ''),
            'Spike_mAbs_inhibitors': 'NA',
            'Spike_Fold':  'NA',
            'RdRp_inhibitors': 'NA',
            'RdRp_Fold':   'NA',
            '3CLpro_inhibitors': 'NA',
            '3CLpro_Fold': 'NA',
            'QC': f"masked: overall coverage={parse_percentish(r.get('coverage'))}"
        })
        continue

    qc_bits = []

    # SPIKE
    if gene_ok(r, 'S', 0.95):
        s_nc   = gather_gene_muts(r, 'S')
        s_disp = to_lookup(s_nc, spike_df)
        s_fold = max_fold(s_nc, spike_df, 'BAM: fold')  # change if you prefer other mAb or max across all
        spike_inhib  = ",".join(s_disp) if s_disp else 'No Mutations'
        spike_fold   = s_fold if s_fold is not None else 'No Data'
    else:
        spike_inhib, spike_fold = 'NA', 'NA'
        qc_bits.append('low S coverage')

    # 3CLpro (ORF1a → nsp5)
    if gene_ok(r, 'ORF1a', 0.95):
        a_nc   = gather_gene_muts(r, 'ORF1a')
        a_disp = to_lookup(a_nc, clpro_df)
        a_fold = max_fold(a_nc, clpro_df, 'NTV: fold')
        clpro_inhib = ",".join(a_disp) if a_disp else 'No Mutations'
        clpro_fold  = a_fold if a_fold is not None else 'No Data'
    else:
        clpro_inhib, clpro_fold = 'NA', 'NA'
        qc_bits.append('low ORF1a coverage')

    # RdRp (ORF1b → nsp12)
    if gene_ok(r, 'ORF1b', 0.95):
        b_nc   = gather_gene_muts(r, 'ORF1b')
        b_disp = to_lookup(b_nc, rdrp_df)
        b_fold = max_fold(b_nc, rdrp_df, 'RDV: fold')
        rdrp_inhib = ",".join(b_disp) if b_disp else 'No Mutations'
        rdrp_fold  = b_fold if b_fold is not None else 'No Data'
    else:
        rdrp_inhib, rdrp_fold = 'NA', 'NA'
        qc_bits.append('low ORF1b coverage')

    rows.append({
        'Sample': r.get('Sample', ''),
        'Spike_mAbs_inhibitors': spike_inhib,
        'Spike_Fold':  spike_fold,
        'RdRp_inhibitors': rdrp_inhib,
        'RdRp_Fold':   rdrp_fold,
        '3CLpro_inhibitors': clpro_inhib,
        '3CLpro_Fold': clpro_fold,
        'QC': "; ".join(qc_bits) if qc_bits else 'ok'
    })

# ---------- write output ----------
out_df = pd.DataFrame(rows)
out_file = f'{run_id}_resistance_mutations.csv'
out_df.to_csv(out_file, index=False)
print(f'✔ Results written to {out_file}')
print(f'   (Rows hard-masked by overall coverage: {masked_all})', file=sys.stderr)
