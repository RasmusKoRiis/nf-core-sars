#!/usr/bin/env python3
"""
Generate a resistance-mutation summary from a Nextclade CSV.
Call:
python resistance_check.py  <nextclade.csv>  <run_id>  <spike.csv>  <rdrp.csv>  <clpro.csv>
"""

import sys, re, math, pandas as pd
main_csv, run_id, spike_csv, rdrp_csv, clpro_csv = sys.argv[1:6]

# ───────── helpers ───────────────────────────────────────────────────
def normalise_del(m: str) -> str:          # 'L24-' → 'L24del'
    return m[:-1] + 'del' if m.strip().endswith('-') else m.strip()

def offset(m: str, d: int) -> str:         # add / subtract residue shift
    mobj = re.fullmatch(r'([A-Z])(\d+)([A-Z]|del)', m)
    if not mobj:
        return m
    aa_from, pos, aa_to = mobj.groups()
    return f'{aa_from}{int(pos)+d}{aa_to}'

def split_clean(raw: str, *, dels=False) -> list[str]:
    if pd.isna(raw) or raw.strip().lower() == 'no mutation':
        return []
    bits = re.split(r'[;,]', raw)
    return [normalise_del(b) if dels else b.strip()
            for b in bits if b.strip()]

def max_fold(nc_muts, df, col):
    vals = []
    for m in nc_muts:
        for v in df.loc[df['Nextclade_lookup'] == m, col]:
            try:
                vals.append(float(str(v).replace(',', '').replace('%', '')))
            except (ValueError, TypeError):
                pass
    return max(vals) if vals else None

def to_lookup(nc_muts, df):          # ➊ convert NC → internal for reporting
    lst = df.loc[df['Nextclade_lookup'].isin(nc_muts), 'Mutation']
    seen = set()
    return [m for m in lst if not (m in seen or seen.add(m))]

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
    # SPIKE
    s_nc   = split_clean(r['S_aaSubstitutions']) + \
             split_clean(r.get('S_aaDeletions', ''), dels=True)
    s_disp = to_lookup(s_nc, spike_df)            # ➊ display list
    s_fold = max_fold(s_nc, spike_df, 'BAM: fold')

    # 3CLpro
    a_nc   = split_clean(r['ORF1a_aaSubstitutions']) + \
             split_clean(r.get('ORF1a_aaDeletions', ''), dels=True)
    a_disp = to_lookup(a_nc, clpro_df)            # ➊
    a_fold = max_fold(a_nc, clpro_df, 'NTV: fold')

    # RdRp
    b_nc   = split_clean(r['ORF1b_aaSubstitutions']) + \
             split_clean(r.get('ORF1b_aaDeletions', ''), dels=True)
    b_disp = to_lookup(b_nc, rdrp_df)             # ➊
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
