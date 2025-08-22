#!/usr/bin/env python3
import sys, glob, re
import pandas as pd

MIN_OVERALL = 80.0  # %
MIN_GENE    = 80.0  # %

def clean_key(s):
    if pd.isna(s): return ""
    s = str(s)
    s = s.replace('\u00A0',' ').replace('\u200b','').replace('\ufeff','')  # NBSP/ZWSP/BOM
    s = re.sub(r'[\u2012\u2013\u2014\u2212]', '-', s)                      # – — − → '-'
    return s.strip()

def parse_percentish(x):
    if pd.isna(x): return None
    s = str(x).strip().replace('%','')
    if not s: return None
    try: v = float(s)
    except ValueError: return None
    return v*100.0 if v <= 1.0 else v

def parse_cds_coverage(s):
    out = {}
    if pd.isna(s): return out
    for part in re.split(r'[;,]\s*', str(s).strip()):
        if not part or ':' not in part: continue
        gene, val = part.split(':', 1)
        gene = gene.strip(); val = val.strip().replace('%','')
        try: v = float(val)
        except ValueError: continue
        if v > 1.0: v /= 100.0
        out[gene] = v
    return out

def classify_dr_from_inhibitors(val):
    if pd.isna(val): return 'NA'
    s = str(val).strip()
    if s == '' or s.upper() == 'NA': return 'NA'
    return 'ANNI' if 'no mutation' in s.lower() else 'Review'

def gene_cols_present(df, gene):
    patt = re.compile(fr'^{gene}_(aaSubstitutions|aaDeletions|aaInsertions)(?:_\d+)?$')
    return [c for c in df.columns if isinstance(c,str) and patt.fullmatch(c)]

def read_and_stack(patterns):
    """Read all CSVs matching any of the patterns; stack vertically by Sample."""
    paths = []
    for pat in patterns:
        paths += glob.glob(pat)
    paths = sorted(set(paths))
    frames = []
    for p in paths:
        try:
            df = pd.read_csv(p, keep_default_na=False)
            if 'Sample' not in df.columns:
                if 'seqName' in df.columns: df = df.rename(columns={'seqName':'Sample'})
                elif 'SequenceID' in df.columns: df = df.rename(columns={'SequenceID':'Sample'})
            if 'Sample' in df.columns:
                df['Sample'] = df['Sample'].map(clean_key)
                frames.append(df)
        except Exception as e:
            print(f"[report] skip {p}: {e}", file=sys.stderr)
    if not frames:
        return pd.DataFrame(columns=['Sample'])
    out = pd.concat(frames, axis=0, ignore_index=True)
    out = out.drop_duplicates(subset=['Sample'])
    out['Sample'] = out['Sample'].map(clean_key)
    return out

def merge_with_suffix(left, right, suffix):
    """Left-join keeping ALL right columns; colliding names from right get suffixed."""
    if right.empty: return left
    # Rename collisions except the join key
    rename = {c: f"{c}{suffix}" for c in right.columns if c in left.columns and c != 'Sample'}
    if rename:
        right = right.rename(columns=rename)
    return left.merge(right, on='Sample', how='left', validate='one_to_one')

# ----------------- load inputs -----------------
# 1) Resistance summaries (one row per sample, wide schema)
res = read_and_stack(["*_resistance_mutations.csv"])
if res.empty:
    sys.exit("[report] No *_resistance_mutations.csv files found.")

# 2) Nextclade stats (you attached these; bring all their columns)
stats = read_and_stack(["*_nextclade_stats.csv"])

# 3) Build union base of sample IDs (res ∪ stats ∪ samplesheet if present)
all_samples = set(res['Sample'])
if not stats.empty: all_samples |= set(stats['Sample'])

base = pd.DataFrame({'Sample': sorted(all_samples)})

# 4) Merge in resistance first (it carries the gene aa* + inhibitors)
rep = merge_with_suffix(base, res, suffix="")  # no suffix

# 5) Merge in stats, suffix any collisions with _stats so nothing gets dropped
rep = merge_with_suffix(rep, stats, suffix="_stats")

# 6) Optional samplesheet merge (supports TSV/CSV/whitespace)
if len(sys.argv) > 1 and sys.argv[1]:
    ss_path = sys.argv[1]
    # detect delimiter
    sep = '\t'
    try:
        head = open(ss_path,'r',encoding='utf-8',errors='ignore').read(4096)
        if '\t' in head: sep = '\t'
        elif ',' in head: sep = ','
        else: sep = r'\s+'   # accidental spaces from bad sed → still works
    except Exception:
        pass
    try:
        ss = pd.read_csv(ss_path, sep=sep, engine=('python' if sep==r'\s+' else 'c'), keep_default_na=False)
        if 'SequenceID' in ss.columns and 'Sample' not in ss.columns:
            ss = ss.rename(columns={'SequenceID':'Sample'})
        if 'Sample' in ss.columns:
            ss['Sample'] = ss['Sample'].map(clean_key)
            ss = ss.drop_duplicates(subset=['Sample'])
            keep = [c for c in ['Sample','PCR-PlatePosition','KonsCt','Barcode'] if c in ss.columns]
            if len(keep) > 1:
                rep = merge_with_suffix(rep, ss[keep], suffix="_ss")
                print(f"[report] merged samplesheet {ss_path} sep='{sep}'", file=sys.stderr)
    except Exception as e:
        print(f"[report] samplesheet skipped: {e}", file=sys.stderr)

# ----------------- 80/80 masking policy -----------------
# discover all gene symbols for masking aa* fields
all_genes = set()
for c in rep.columns:
    m = re.match(r'^([A-Za-z0-9]+)_(aaSubstitutions|aaDeletions|aaInsertions)', str(c))
    if m: all_genes.add(m.group(1))

def row_overall_pct(row):
    # prefer QC_overall_pct if present; else parse coverage (res or stats)
    for k in ('QC_overall_pct', 'coverage', 'coverage_stats'):
        if k in rep.columns:
            val = row.get(k)
            if str(val).strip() not in ('', 'NA'):
                try: return float(val) if k=='QC_overall_pct' else parse_percentish(val)
                except: pass
    return None

def row_cds_map(row):
    # Prefer explicit QC_cds_* fields; else parse cdsCoverage (if present)
    m = {}
    if 'cdsCoverage' in rep.columns:
        m = parse_cds_coverage(row.get('cdsCoverage'))
    for g, col in [('S','QC_cds_S'), ('ORF1a','QC_cds_ORF1a'), ('ORF1b','QC_cds_ORF1b')]:
        if col in rep.columns:
            v = row.get(col)
            if str(v).strip() not in ('', 'NA'):
                try: m[g] = float(v)  # 0..1
                except: pass
    return m

rep['_overall_pct'] = rep.apply(row_overall_pct, axis=1)
rep['_cds_map']     = rep.apply(row_cds_map, axis=1)

for idx, row in rep.iterrows():
    overall_pct = row['_overall_pct']
    cds_map = row['_cds_map'] or {}
    overall_fail = (overall_pct is None) or (overall_pct < MIN_OVERALL)

    if overall_fail:
        # mask ALL Nextclade aa* and resistance outputs
        for gene in all_genes:
            for col in gene_cols_present(rep, gene):
                rep.at[idx, col] = 'NA'
        for col in ['Spike_mAbs_inhibitors','Spike_Fold',
                    '3CLpro_inhibitors','3CLpro_Fold',
                    'RdRp_inhibitors','RdRp_Fold']:
            if col in rep.columns: rep.at[idx, col] = 'NA'
        if 'QC' in rep.columns:
            q = str(rep.at[idx,'QC']).strip() if not pd.isna(rep.at[idx,'QC']) else ''
            if not q: rep.at[idx,'QC'] = 'masked: overall <80 or Nextclade failed'
        continue

    # Per-gene masking
    for gene, outs in [('S',['Spike_mAbs_inhibitors','Spike_Fold']),
                       ('ORF1a',['3CLpro_inhibitors','3CLpro_Fold']),
                       ('ORF1b',['RdRp_inhibitors','RdRp_Fold'])]:
        v = cds_map.get(gene)
        gene_fail = (v is None) or (v < (MIN_GENE/100.0))
        if gene_fail:
            for col in gene_cols_present(rep, gene):
                rep.at[idx, col] = 'NA'
            for col in outs:
                if col in rep.columns: rep.at[idx, col] = 'NA'
            if 'QC' in rep.columns:
                q = str(rep.at[idx,'QC']).strip() if not pd.isna(rep.at[idx,'QC']) else ''
                add = f"low {gene} cds"
                rep.at[idx,'QC'] = add if not q else (q if add in q else f"{q}; {add}")

# ----------------- DR flags (NA when not analyzed) -----------------
for out_col, src_col in {
    'DR_Res_Paxlovid':   '3CLpro_inhibitors',
    'DR_Res_Remdesevir': 'RdRp_inhibitors',
    'DR_Res_mAbs':       'Spike_mAbs_inhibitors'
}.items():
    rep[out_col] = rep[src_col].apply(classify_dr_from_inhibitors) if src_col in rep.columns else 'NA'

# ----------------- normalize blanks to "NA" in key fields -----------------
NA_ALWAYS = {
    'Spike_mAbs_inhibitors','Spike_Fold',
    '3CLpro_inhibitors','3CLpro_Fold',
    'RdRp_inhibitors','RdRp_Fold',
    'DR_Res_Paxlovid','DR_Res_Remdesevir','DR_Res_mAbs',
    'QC','coverage','cdsCoverage',
    'PCR-PlatePosition','KonsCt','Barcode'
}
GENE_AA_COL = re.compile(
    r'^(?:S|ORF1a|ORF1b|E|M|N|ORF3a|ORF6|ORF7a|ORF7b|ORF8|ORF9b)_'
    r'(?:aaSubstitutions|aaDeletions|aaInsertions)(?:_\d+)?$'
)
def blank_to_NA(x):
    if x is None or pd.isna(x): return 'NA'
    s = str(x).strip()
    return 'NA' if s == '' else x

for c in rep.columns:
    if (c in NA_ALWAYS) or (isinstance(c,str) and GENE_AA_COL.match(c)):
        rep[c] = rep[c].apply(blank_to_NA)

# Write with Sample first
cols = list(rep.columns)
rep = rep[['Sample'] + [c for c in cols if c != 'Sample']]
rep.to_csv("merged_report.csv", index=False)
print("[report] wrote merged_report.csv", file=sys.stderr)
