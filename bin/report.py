#!/usr/bin/env python3
import sys, glob, re
import pandas as pd

MIN_OVERALL = 80.0   # overall coverage threshold (%)
MIN_GENE    = 80.0   # cds coverage threshold (%)

def clean_key(s):
    if pd.isna(s): return ""
    s = str(s)
    s = s.replace('\u00A0',' ').replace('\u200b','').replace('\ufeff','')  # NBSP/ZWSP/BOM
    s = re.sub(r'[\u2012\u2013\u2014\u2212]', '-', s)                      # – — − → '-'
    return s.strip()

def parse_percentish(x):
    """Return percent 0–100 or None; handles '0.956', '95.6', '95.6%'."""
    if pd.isna(x): return None
    s = str(x).strip().replace('%','')
    if not s: return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v*100.0 if v <= 1.0 else v

def parse_cds_coverage(s):
    """cdsCoverage 'E:1;...;ORF1a:0.9;ORF1b:1;S:0.87' → dict floats 0..1"""
    out = {}
    if pd.isna(s): return out
    for part in re.split(r'[;,]\s*', str(s).strip()):
        if not part or ':' not in part: continue
        gene, val = part.split(':', 1)
        gene = gene.strip()
        val  = val.strip().replace('%','')
        try:
            v = float(val)
        except ValueError:
            continue
        if v > 1.0: v = v/100.0
        out[gene] = v
    return out

def classify_dr_from_inhibitors(val: object) -> str:
    """Return NA if not analyzed, ANNI if 'No Mutations', else Review."""
    if pd.isna(val):           return 'NA'
    s = str(val).strip()
    if s == '' or s.upper() == 'NA':  # not analyzed → NA
        return 'NA'
    if 'no mutation' in s.lower():
        return 'ANNI'
    return 'Review'

def gene_cols_present(df, gene):
    """All Nextclade aa* columns for a gene (incl. split _1.._4 if present)."""
    cols = []
    patt = re.compile(fr'^{gene}_(aaSubstitutions|aaDeletions|aaInsertions)(?:_\d+)?$')
    for c in df.columns:
        if isinstance(c, str) and patt.fullmatch(c):
            cols.append(c)
    return cols

if len(sys.argv) < 1:
    sys.exit(f"Usage: {sys.argv[0]} <samplesheet.tsv|csv>")

# Allow no samplesheet (some runs may not have it); if provided, merge it later
samplesheet_path = sys.argv[1] if len(sys.argv) > 1 else None

# 1) Load ALL per-sample resistance CSVs and CONCAT VERTICALLY
res_paths = sorted(glob.glob("*_resistance_mutations.csv"))
if not res_paths:
    sys.exit("[report] No *_resistance_mutations.csv files found.")

frames = []
for p in res_paths:
    try:
        df = pd.read_csv(p, keep_default_na=False)  # keep literal "NA"
        if 'Sample' not in df.columns and 'seqName' in df.columns:
            df = df.rename(columns={'seqName':'Sample'})
        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].map(clean_key)
            frames.append(df)
    except Exception as e:
        print(f"[report] skip {p}: {e}", file=sys.stderr)

if not frames:
    sys.exit("[report] No CSVs with a 'Sample' column after loading.")

rep = pd.concat(frames, axis=0, ignore_index=True)
rep = rep.drop_duplicates(subset=['Sample']).copy()
rep['Sample'] = rep['Sample'].map(clean_key)

# Helper: per-row coverage indicators
def row_overall_pct(row):
    # Prefer 'QC_overall_pct' if present; else parse 'coverage'
    if 'QC_overall_pct' in row and row['QC_overall_pct'] not in ('', None):
        try:
            return float(row['QC_overall_pct'])
        except Exception:
            pass
    return parse_percentish(row.get('coverage'))

def row_cds_map(row):
    # Prefer QC_cds_* fields if present; else parse cdsCoverage blob
    m = {}
    cds_blob = row.get('cdsCoverage')
    m = parse_cds_coverage(cds_blob)
    # overwrite with explicit QC columns if they exist
    for g, col in [('S','QC_cds_S'), ('ORF1a','QC_cds_ORF1a'), ('ORF1b','QC_cds_ORF1b')]:
        v = row.get(col)
        try:
            if v != '' and v is not None:
                m[g] = float(v)  # these are 0..1
        except Exception:
            pass
    return m

# 2) Apply masking policy (80/80)
all_genes = set()
# Discover all genes present in the table (for masking their aa* columns)
for c in rep.columns:
    m = re.match(r'^([A-Za-z0-9]+)_(aaSubstitutions|aaDeletions|aaInsertions)', str(c))
    if m:
        all_genes.add(m.group(1))

# Compute coverage flags
overall_list = []
cds_maps = []
for _, row in rep.iterrows():
    overall_list.append(row_overall_pct(row))
    cds_maps.append(row_cds_map(row))

rep['_overall_pct'] = overall_list  # percent or None
rep['_cds_map']     = cds_maps      # dict of 0..1

# Mask per-row
for idx, row in rep.iterrows():
    overall_pct = row['_overall_pct']  # None or %
    cds_map = row['_cds_map'] or {}

    overall_fail = (overall_pct is None) or (overall_pct < MIN_OVERALL)

    # If sample failed Nextclade / overall <80 → mark NOT ANALYZED
    if overall_fail:
        # Mask ALL gene aa* columns
        for gene in all_genes:
            for col in gene_cols_present(rep, gene):
                rep.at[idx, col] = 'NA'
        # Mask resistance outputs and ensure DR flags later become NA
        for col in ['Spike_mAbs_inhibitors','Spike_Fold',
                    '3CLpro_inhibitors','3CLpro_Fold',
                    'RdRp_inhibitors','RdRp_Fold']:
            if col in rep.columns:
                rep.at[idx, col] = 'NA'
        # You may also want to set a clear QC label
        if 'QC' in rep.columns and (not str(rep.at[idx,'QC']).strip()):
            rep.at[idx, 'QC'] = 'masked: overall <80 or Nextclade failed'
        continue  # no per-gene masking needed; already NA’d everything

    # Else: per-gene cds masking (<80% → NA for that gene)
    for gene, thresh in [('S', MIN_GENE), ('ORF1a', MIN_GENE), ('ORF1b', MIN_GENE)]:
        v = cds_map.get(gene)  # 0..1 or None
        gene_fail = (v is None) or (v < (thresh/100.0))
        if gene_fail:
            # Mask gene's aa* columns
            for col in gene_cols_present(rep, gene):
                rep.at[idx, col] = 'NA'
            # Mask resistance outputs for that target
            if gene == 'S':
                for col in ['Spike_mAbs_inhibitors','Spike_Fold']:
                    if col in rep.columns: rep.at[idx, col] = 'NA'
            elif gene == 'ORF1a':
                for col in ['3CLpro_inhibitors','3CLpro_Fold']:
                    if col in rep.columns: rep.at[idx, col] = 'NA'
            elif gene == 'ORF1b':
                for col in ['RdRp_inhibitors','RdRp_Fold']:
                    if col in rep.columns: rep.at[idx, col] = 'NA'
            # Set QC hint
            if 'QC' in rep.columns:
                q = str(rep.at[idx,'QC']).strip()
                add = f"low {gene} cds"
                rep.at[idx,'QC'] = add if not q else (q if add in q else f"{q}; {add}")

# 3) Merge samplesheet (optional)
if len(sys.argv) > 1 and sys.argv[1]:
    samplesheet_path = sys.argv[1]
    # detect delimiter
    delim = '\t'
    try:
        head = open(samplesheet_path, 'r', encoding='utf-8', errors='ignore').read(4096)
        if '\t' not in head: delim = ','
    except Exception:
        pass
    try:
        ss = pd.read_csv(samplesheet_path, sep=delim, keep_default_na=False)
        if 'SequenceID' in ss.columns and 'Sample' not in ss.columns:
            ss = ss.rename(columns={'SequenceID':'Sample'})
        if 'Sample' in ss.columns:
            ss['Sample'] = ss['Sample'].map(clean_key)
            ss = ss.drop_duplicates(subset=['Sample'])
            keep = [c for c in ['Sample','PCR-PlatePosition','KonsCt','Barcode'] if c in ss.columns]
            if len(keep) > 1:
                rep = rep.merge(ss[keep], on='Sample', how='left', validate='one_to_one')
                print(f"[report] merged samplesheet {samplesheet_path} sep='{delim}'", file=sys.stderr)
    except Exception as e:
        print(f"[report] samplesheet skipped: {e}", file=sys.stderr)

# 4) DR flags — return NA when inhibitors are NA (not analyzed)
for out_col, src_col in {
    'DR_Res_Paxlovid':   '3CLpro_inhibitors',
    'DR_Res_Remdesevir': 'RdRp_inhibitors',
    'DR_Res_mAbs':       'Spike_mAbs_inhibitors'
}.items():
    rep[out_col] = rep[src_col].apply(classify_dr_from_inhibitors) if src_col in rep.columns else 'NA'

# 5) Tidy columns: Sample first, keep literal "NA"
cols = list(rep.columns)
rep = rep[['Sample'] + [c for c in cols if c != 'Sample']]

rep.to_csv("merged_report.csv", index=False)
print("[report] wrote merged_report.csv", file=sys.stderr)
