#!/usr/bin/env python3
"""
Generate a resistance-mutation summary from a Nextclade CSV with QC gating.
Call:
python resistance_check.py  <nextclade.csv>  <run_id>  <spike.csv>  <rdrp.csv>  <clpro.csv>
"""

import sys, re, math, pandas as pd
main_csv, run_id, spike_csv, rdrp_csv, clpro_csv = sys.argv[1:6]

# ───────── helpers ───────────────────────────────────────────────────
def normalise_del(m: str) -> str:          # 'L24-' → 'L24del'
    return m[:-1] + 'del' if isinstance(m, str) and m.strip().endswith('-') else str(m).strip()

def offset(m: str, d: int) -> str:         # add / subtract residue shift
    mobj = re.fullmatch(r'([A-Z])(\d+)([A-Z]|del)', m)
    if not mobj:
        return m
    aa_from, pos, aa_to = mobj.groups()
    return f'{aa_from}{int(pos)+d}{aa_to}'

def split_clean(raw: str, *, dels=False) -> list[str]:
    if pd.isna(raw):
        return []
    s = str(raw).strip()
    if not s or s.lower() == 'no mutation':
        return []
    bits = re.split(r'[;,]', s)
    return [normalise_del(b) if dels else b.strip() for b in bits if str(b).strip()]

def max_fold(nc_muts, df, col):
    vals = []
    for m in nc_muts:
        hits = df.loc[df['Nextclade_lookup'] == m, col]
        for v in hits:
            try:
                vals.append(float(str(v).replace(',', '').replace('%', '')))
            except (ValueError, TypeError):
                pass
    return max(vals) if vals else None

def to_lookup(nc_muts, df):          # convert NC → internal for reporting
    lst = df.loc[df['Nextclade_lookup'].isin(nc_muts), 'Mutation']
    seen = set()
    out = []
    for m in lst:
        if m not in seen:
            out.append(m)
            seen.add(m)
    return out

# --- QC helpers ---
def parse_percentish(x) -> float | None:
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
    # If value looks like a fraction (<=1), treat as 0–1 → %
    return v * 100.0 if v <= 1.0 else v

def parse_cds_coverage(s: str) -> dict:
    """
    Parse cdsCoverage like 'E:1,M:1,...,ORF1a:0.9343,ORF1b:1,S:1'
    Returns floats in 0–1; tolerates '%', or 0–100 values.
    """
    out = {}
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
        # Normalize to 0–1
        if v > 1.0:
            v = v / 100.0
        out[gene] = v
    return out

def should_mask(row) -> bool:
    cov = parse_percentish(row.get('coverage'))
    if cov is None or cov < 90.0:
        return True
    # Handle possible typo 'cdcCoverage' vs 'cdsCoverage'
    cds_raw = row.get('cdsCoverage', row.get('cdcCoverage'))
    cds = parse_cds_coverage(cds_raw)
    for gene in ('ORF1a', 'ORF1b', 'S'):
        v = cds.get(gene)
        if v is None or v < 0.95:
            return True
    return False

# ───────── load look-up tables ───────────────────────────────────────
clpro_df = pd.read_csv(clpro_csv)
rdrp_df  = pd.read_csv(rdrp_csv)
spike_df = pd.read_csv(spike_csv)

clpro_df['Nextclade_lookup'] = clpro_df['Mutation'].apply(
    lambda m: offset(normalise_del(m),  3263))   # nsp5 +3263
rdrp_df ['Nextclade_lookup'] = rdrp_df ['Mutation'].apply(
    lambda m: offset(normalise_del(m),   -9))    # nsp12 −9
spike_df['Nextclade_lookup'] = spike_df['Mutation']

# ───────── iterate samples ───────────────────────────────────────────
main_df = pd.read_csv(main_csv)
rows = []

for _, r in main_df.iterrows():
    mask = should_mask(r)

    if mask:
        rows.append({
            'Sample': r['Sample'],
            'Spike_mAbs_inhibitors': 'NA',
            'Spike_Fold':  'NA',
            'RdRp_inhibitors': 'NA',
            'RdRp_Fold':   'NA',
            '3CLpro_inhibitors': 'NA',
            '3CLpro_Fold': 'NA',
        })
        continue

    # SPIKE
    s_nc   = split_clean(r.get('S_aaSubstitutions')) + split_clean(r.get('S_aaDeletions', ''), dels=True)
    s_disp = to_lookup(s_nc, spike_df)
    s_fold = max_fold(s_nc, spike_df, 'BAM: fold')

    # 3CLpro (ORF1a → nsp5)
    a_nc   = split_clean(r.get('ORF1a_aaSubstitutions')) + split_clean(r.get('ORF1a_aaDeletions', ''), dels=True)
    a_disp = to_lookup(a_nc, clpro_df)
    a_fold = max_fold(a_nc, clpro_df, 'NTV: fold')

    # RdRp (ORF1b → nsp12)
    b_nc   = split_clean(r.get('ORF1b_aaSubstitutions')) + split_clean(r.get('ORF1b_aaDeletions', ''), dels=True)
    b_disp = to_lookup(b_nc, rdrp_df)
    b_fold = max_fold(b_nc, rdrp_df, 'RDV: fold')

    rows.append({
        'Sample': r['Sample'],
        'Spike_mAbs_inhibitors': ",".join(s_disp) if s_disp else 'No Mutations',
        'Spike_Fold':  s_fold if s_fold is not None else 'No Data',
        'RdRp_inhibitors': ",".join(b_disp) if b_disp else 'No Mutations',
        'RdRp_Fold':   b_fold if b_fold is not None else 'No Data',
        '3CLpro_inhibitors': ",".join(a_disp) if a_disp else 'No Mutations',
        '3CLpro_Fold': a_fold if a_fold is not None else 'No Data',
    })

# ───────── write output ──────────────────────────────────────────────
out_df = pd.DataFrame(rows)
out_file = f'{run_id}_resistance_mutations.csv'
out_df.to_csv(out_file, index=False)
print(f'✔ Results written to {out_file}')
