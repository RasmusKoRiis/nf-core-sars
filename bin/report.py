#!/usr/bin/env python3
import sys, glob, re
import pandas as pd

def clean_key(s):
    if pd.isna(s): return ""
    s = str(s)
    s = s.replace('\u00A0',' ').replace('\u200b','').replace('\ufeff','')  # NBSP/ZWSP/BOM
    s = re.sub(r'[\u2012\u2013\u2014\u2212]', '-', s)                      # –—− → '-'
    return s.strip()

def classify_dr(val: object) -> str:
    if pd.isna(val): return 'NA'
    s = str(val).strip().lower()
    if s == '' or any(k in s for k in [
        'na','n/a','not available','not applicable','not analyzed','not analysed',
        'analysis not possible','analysis failed','insufficient','low coverage','no data','unknown','error'
    ]): return 'NA'
    if 'no mutation' in s: return 'ANNI'
    return 'Review'

if len(sys.argv) < 2:
    sys.exit(f"Usage: {sys.argv[0]} <samplesheet.tsv|csv>")

samplesheet_path = sys.argv[1]

# ---- 1) Load per-sample resistance CSVs and CONCAT vertically ----
# Prefer explicit pattern; fallback to all CSVs if none match.
res_paths = sorted(glob.glob("*_resistance_mutations.csv"))
if not res_paths:
    res_paths = sorted(glob.glob("*.csv"))

frames = []
for p in res_paths:
    try:
        df = pd.read_csv(p, keep_default_na=False)  # keep literal "NA"
        # Normalize key name if needed
        if 'Sample' not in df.columns:
            if 'seqName' in df.columns:
                df = df.rename(columns={'seqName': 'Sample'})
            elif 'SequenceID' in df.columns:
                df = df.rename(columns={'SequenceID': 'Sample'})
        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].map(clean_key)
            frames.append(df)
    except Exception as e:
        print(f"[report] skip {p}: {e}", file=sys.stderr)

if not frames:
    sys.exit("[report] No CSVs with a 'Sample' column found.")

# Vertically stack (this is the key change!)
res_base = pd.concat(frames, axis=0, ignore_index=True)
res_base = res_base.drop_duplicates(subset=['Sample'])  # one row per sample
res_base['Sample'] = res_base['Sample'].map(clean_key)

# ---- 2) Merge in samplesheet (supports CSV or TSV) ----
# Detect delimiter: tab if any tabs, else comma
delim = '\t'
try:
    with open(samplesheet_path, 'r', encoding='utf-8', errors='ignore') as fh:
        head = fh.read(4096)
    if '\t' not in head:
        delim = ','
except Exception:
    pass

try:
    ss = pd.read_csv(samplesheet_path, sep=delim, keep_default_na=False)
    if 'SequenceID' in ss.columns and 'Sample' not in ss.columns:
        ss = ss.rename(columns={'SequenceID': 'Sample'})
    if 'Sample' in ss.columns:
        ss['Sample'] = ss['Sample'].map(clean_key)
        ss = ss.drop_duplicates(subset=['Sample'])
        keep = [c for c in ['Sample','PCR-PlatePosition','KonsCt','Barcode'] if c in ss.columns]
        if len(keep) > 1:
            res_base = res_base.merge(ss[keep], on='Sample', how='left', validate='one_to_one')
            print(f"[report] merged samplesheet ({samplesheet_path}) with sep='{delim}'", file=sys.stderr)
except Exception as e:
    print(f"[report] samplesheet skipped: {e}", file=sys.stderr)

# ---- 3) Add DR flags from resistance columns ----
for out_col, src_col in {
    'DR_Res_Paxlovid':   '3CLpro_inhibitors',
    'DR_Res_Remdesevir': 'RdRp_inhibitors',
    'DR_Res_mAbs':       'Spike_mAbs_inhibitors'
}.items():
    res_base[out_col] = res_base[src_col].apply(classify_dr) if src_col in res_base.columns else 'NA'

# ---- 4) Final ordering and write ----
cols = list(res_base.columns)
res_base = res_base[['Sample'] + [c for c in cols if c != 'Sample']]
res_base.to_csv("merged_report.csv", index=False)
print("[report] wrote merged_report.csv", file=sys.stderr)
