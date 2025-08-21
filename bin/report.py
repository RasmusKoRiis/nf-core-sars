#!/usr/bin/env python3
import sys, glob, re
import pandas as pd
from functools import reduce

def clean_key(s):
    if pd.isna(s): return ""
    s = str(s)
    # normalize weird spaces/dashes
    s = s.replace('\u00A0',' ').replace('\u200b','').replace('\ufeff','')
    s = re.sub(r'[\u2012\u2013\u2014\u2212]', '-', s)  # figure/en/em/minus→'-'
    return s.strip()

if len(sys.argv) < 2:
    sys.exit(f"Usage: {sys.argv[0]} <samplesheet.tsv>")

samplesheet = sys.argv[1]

# Load all CSVs in the folder, keep literal "NA"
csv_paths = glob.glob("*.csv")
frames = []
for p in csv_paths:
    try:
        df = pd.read_csv(p, keep_default_na=False)  # "NA" stays "NA"
        if 'Sample' not in df.columns and 'seqName' in df.columns:
            df = df.rename(columns={'seqName':'Sample'})
        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].map(clean_key)
            frames.append(df)
    except Exception as e:
        print(f"[merge] skip {p}: {e}", file=sys.stderr)

if not frames:
    sys.exit("[merge] No CSVs with a 'Sample' column found.")

# Use the widest frame as the base (likely the masked Nextclade+resistance CSV)
base = max(frames, key=lambda d: d.shape[1])
others = [f for f in frames if f is not base]

# Deduplicate keys per frame
base = base.drop_duplicates(subset=['Sample'])
others = [f.drop_duplicates(subset=['Sample']) for f in others]

# Merge everyone onto the base, one_by_one; 'left' keeps base intact
def merge1(left, right):
    common = [c for c in right.columns if c in left.columns and c != 'Sample']
    # To avoid overwriting, suffix right columns that collide
    suffix = "_x"
    if common:
        right = right.rename(columns={c: f"{c}{suffix}" for c in common})
    return left.merge(right, on='Sample', how='left', validate='one_to_one')

merged = reduce(merge1, others, base)

# Load samplesheet (if present) and merge selected fields
try:
    ss = pd.read_csv(samplesheet, sep="\t", keep_default_na=False)
    ss = ss.rename(columns={'SequenceID':'Sample'})
    if 'Barcode' in ss.columns: ss = ss.drop(columns=['Barcode'])
    ss['Sample'] = ss['Sample'].map(clean_key)
    ss = ss.drop_duplicates(subset=['Sample'])
    pick = [c for c in ['Sample','PCR-PlatePosition','KonsCt'] if c in ss.columns]
    if len(pick) > 1:
        merged = merged.merge(ss[pick], on='Sample', how='left', validate='one_to_one')
except Exception as e:
    print(f"[merge] samplesheet skipped: {e}", file=sys.stderr)

# DR flags
def classify_dr(val: object) -> str:
    if pd.isna(val): return 'NA'
    s = str(val).strip().lower()
    if s == '' or any(k in s for k in [
        'na','n/a','not available','not applicable','not analyzed','not analysed',
        'analysis not possible','analysis failed','insufficient','low coverage','no data','unknown','error'
    ]): return 'NA'
    if 'no mutation' in s: return 'ANNI'
    return 'Review'

for out_col, src_col in {
    'DR_Res_Paxlovid':   '3CLpro_inhibitors',
    'DR_Res_Remdesevir': 'RdRp_inhibitors',
    'DR_Res_mAbs':       'Spike_mAbs_inhibitors'
}.items():
    merged[out_col] = merged[src_col].apply(classify_dr) if src_col in merged.columns else 'NA'

# Put Sample first
other_cols = [c for c in merged.columns if c != 'Sample']
merged = merged[['Sample'] + other_cols]

merged.to_csv("merged_report.csv", index=False)
print("[merge] wrote merged_report.csv", file=sys.stderr)
