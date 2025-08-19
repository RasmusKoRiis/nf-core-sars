#!/usr/bin/env python3
import pandas as pd
import glob
import sys

# --------
# MERGE ALL DATA
# --------
csv_files = glob.glob('*.csv')
samplesheet = sys.argv[1]

# LOAD SAMPLE SHEET (note: currently unused; keep or merge if needed later)
samplesheet_df = pd.read_csv(samplesheet, sep='\t')
samplesheet_df.rename(columns={'SequenceID': 'Sample'}, inplace=True)
if 'Barcode' in samplesheet_df.columns:
    samplesheet_df = samplesheet_df.drop(columns=['Barcode'])

merged_data = None
for i, file in enumerate(csv_files):
    df = pd.read_csv(file)
    merged_data = df if merged_data is None else pd.concat([merged_data, df], axis=0, ignore_index=True)

# Collapse duplicate Samples, keep first occurrence
merged_data = merged_data.groupby('Sample', as_index=False).first()

# NIPH specific adjustment: replace "!" with "-" in Sample
merged_data['Sample'] = merged_data['Sample'].astype(str).str.replace('!', '-', regex=False)

# --------
# DR RESISTANCE FLAGS
# --------
def classify_dr(val: object) -> str:
    """Return 'ANNI' if 'No mutation(s)', 'Review' if mutations listed, else 'NA'."""
    if pd.isna(val):
        return 'NA'
    s = str(val).strip().lower()
    if s == '' or any(k in s for k in [
        'na', 'n/a', 'not available', 'not applicable', 'not analyzed', 'not analysed',
        'analysis not possible', 'analysis failed', 'insufficient', 'low coverage',
        'no data', 'unknown', 'error'
    ]):
        return 'NA'
    if 'no mutation' in s:  # matches "no mutation" or "no mutations"
        return 'ANNI'
    return 'Review'

mapping = {
    'DR_Res_Paxlovid':   '3CLpro_inhibitors',       # nirmatrelvir target: 3CLpro
    'DR_Res_Remdesevir': 'RdRp_inhibitors',         # remdesivir target: RdRp
    'DR_Res_mAbs':       'Spike_mAbs_inhibitors'    # mAbs target: Spike
}

for out_col, src_col in mapping.items():
    if src_col in merged_data.columns:
        merged_data[out_col] = merged_data[src_col].apply(classify_dr)
    else:
        merged_data[out_col] = 'NA'  # source column missing -> analysis not possible

# Sort columns: Sample first, then alphabetical
other_cols = sorted([c for c in merged_data.columns if c != 'Sample'])
merged_data = merged_data[['Sample'] + other_cols]

# --------
# WRITE OUTPUT
# --------
merged_data.to_csv('merged_report.csv', index=False)
