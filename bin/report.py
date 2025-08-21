#!/usr/bin/env python3
import sys, glob, re, os
import pandas as pd
from functools import reduce

def clean_key(s):
    if pd.isna(s): return ""
    s = str(s)
    s = s.replace('\u00A0',' ').replace('\u200b','').replace('\ufeff','')  # NBSP/ZWSP/BOM
    s = re.sub(r'[\u2012\u2013\u2014\u2212]', '-', s)                      # figure/en/em/minus → '-'
    return s.strip()

if len(sys.argv) < 2:
    sys.exit(f"Usage: {sys.argv[0]} <samplesheet.tsv>")

samplesheet = sys.argv[1]

# Load all CSVs present (literal "NA" preserved)
csv_paths = sorted(glob.glob("*.csv"))

frames = []
for p in csv_paths:
    try:
        df = pd.read_csv(p, keep_default_na=False)
        if 'Sample' not in df.columns and 'seqName' in df.columns:
            df = df.rename(columns={'seqName': 'Sample'})
        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].map(clean_key)
            df = df.drop_duplicates(subset=['Sample'])
            frames.append((p, df))
    except Exception as e:
        print(f"[merge] skip {p}: {e}", file=sys.stderr)

if not frames:
    sys.exit("[merge] No CSVs with a 'Sample' column found.")

# Prefer a base that already contains resistance columns; else the widest
def score(fr):
    name, df = fr
    has_res = sum(c in df.columns for c in (
        'Spike_mAbs_inhibitors','3CLpro_inhibitors','RdRp_inhibitors'))
    return (has_res, df.shape[1])

base_name, base = max(frames, key=score)
others = [(n,d) for (n,d) in frames if n != base_name]
print(f"[merge] base: {base_name} ({base.shape[0]}x{base.shape[1]})", file=sys.stderr)

def merge_safe(left, right):
    # Drop overlapping columns from right (except the join key) to avoid suffix chaos
    drop_cols = [c for c in right.columns if c in left.columns and c != 'Sample']
    if drop_cols:
        right = right.drop(columns=drop_cols)
    return left.merge(right, on='Sample', how='left', validate='one_to_one')

merged = base
for name, df in others:
    before_cols = merged.shape[1]
    merged = merge_safe(merged, df)
    print(f"[merge] joined {name} → cols {before_cols} → {merged.shape[1]}", file=sys.stderr)

# Bring in samplesheet fields (optional)
try:
    ss = pd.read_csv(samplesheet, sep="\t", keep_default_na=False)
    ss = ss.rename(columns={'SequenceID': 'Sample'})
    ss['Sample'] = ss['Sample'].map(clean_key)
    ss = ss.drop_duplicates(subset=['Sample'])
    keep = [c for c in ['Sample','PCR-PlatePosition','KonsCt'] if c in ss.columns]
    if len(keep) > 1:
        merged = merge_safe(merged, ss[keep])
        print(f"[merge] merged samplesheet: {samplesheet}", file=sys.stderr)
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

# Column order: Sample first, then everything else (stable order)
cols = list(merged.columns)
cols = ['Sample'] + [c for c in cols if c != 'Sample']
merged = merged[cols]

merged.to_csv("merged_report.csv", index=False)
print("[merge] wrote merged_report.csv", file=sys.stderr)
