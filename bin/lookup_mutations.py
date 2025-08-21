#!/usr/bin/env python3
"""
Generate resistance summary AND mask Nextclade fields based on coverage thresholds.

Policy:
- If overall coverage < MIN_OVERALL_PCT (default 80): mask ALL Nextclade fields to 'NA'
  (keep Sample, coverage, cdsCoverage for transparency), and set all resistance fields 'NA'.
- Else, for each gene, if cds coverage < MIN_GENE_CDS (default 80): mask that gene's
  Nextclade aa fields to 'NA' and skip resistance lookup for that gene.

Call:
python resistance_check.py <nextclade.csv> <run_id> <spike.csv> <rdrp.csv> <clpro.csv>

Env overrides:
  MIN_OVERALL_PCT   default 80
  MIN_GENE_CDS      default 80
"""
import sys, os, re
import pandas as pd
from typing import Optional, List, Dict, Tuple

__VERSION__ = "resistance_check v2025-08-21-80pct"

if len(sys.argv) < 6:
    sys.exit(f"Usage: {sys.argv[0]} <nextclade.csv> <run_id> <spike.csv> <rdrp.csv> <clpro.csv>")

main_csv, run_id, spike_csv, rdrp_csv, clpro_csv = sys.argv[1:6]
MIN_OVERALL = float(os.getenv("MIN_OVERALL_PCT", "80"))
MIN_GENE    = float(os.getenv("MIN_GENE_CDS",   "80"))

print(f"[{__VERSION__}] main={main_csv} run={run_id} "
      f"gates: overall>={MIN_OVERALL}%, gene>={MIN_GENE}%", file=sys.stderr)

# ---------- helpers ----------
def normalise_del(m: str) -> str:  # 'L24-' → 'L24del'
    s = str(m).strip()
    return s[:-1] + 'del' if s.endswith('-') else s

def alt_del_forms(m: str) -> List[str]:
    """Return both 'L24del' and 'L24-' for deletion matches; otherwise [m]."""
    s = str(m).strip()
    out = [s]
    if s.endswith('-'):
        d = normalise_del(s)
        if d not in out: out.append(d)
    m2 = re.fullmatch(r'([A-Z])(\d+)del', s)
    if m2:
        dash = f"{m2.group(1)}{m2.group(2)}-"
        if dash not in out: out.append(dash)
    return out

def offset(m: str, d: int) -> str:
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
        if not b: continue
        out.append(normalise_del(b) if dels else b)
    return out

def gather_gene_muts(row: pd.Series, gene: str) -> List[str]:
    """Collect per-gene AA changes (subs + dels + ins), supporting split columns _1.._4.
       Return list WITHOUT gene prefix (e.g. ['L24del','E484K'])."""
    muts: List[str] = []
    for col in row.index:
        c = str(col)
        if re.fullmatch(fr'{gene}_aaSubstitutions(?:_\d+)?', c):
            muts.extend(split_clean(row.get(col, '')))
        elif re.fullmatch(fr'{gene}_aaDeletions(?:_\d+)?', c):
            muts.extend([normalise_del(x) for x in split_clean(row.get(col, ''), dels=True)])
        elif re.fullmatch(fr'{gene}_aaInsertions(?:_\d+)?', c):
            muts.extend(split_clean(row.get(col, '')))
    # Deduplicate preserving order + add alt deletion forms
    seen = set(); uniq: List[str] = []
    for m in muts:
        if m not in seen:
            uniq.append(m); seen.add(m)
    with_alts: List[str] = []
    seen2 = set()
    for m in uniq:
        for a in alt_del_forms(m):
            if a not in seen2:
                with_alts.append(a); seen2.add(a)
    return with_alts

def guess_fold_col(df: pd.DataFrame, prefer: List[str]) -> Optional[str]:
    cands = [c for c in df.columns if 'fold' in str(c).lower()]
    for p in prefer:
        for c in cands:
            if str(c).upper().startswith(p.upper()):
                return c
    return cands[0] if cands else None

def max_fold(nc_muts: List[str], df: pd.DataFrame, col: Optional[str]) -> Optional[float]:
    if not col or col not in df.columns:
        return None
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
    if pd.isna(x): return None
    s = str(x).strip().replace('%','')
    if not s: return None
    try: v = float(s)
    except ValueError: return None
    return v * 100.0 if v <= 1.0 else v

def parse_cds_coverage(s: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
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

def overall_ok(row: pd.Series, min_pct: float = MIN_OVERALL) -> Tuple[bool, Optional[float]]:
    cov = parse_percentish(row.get('coverage'))
    return (True if cov is None else (cov >= min_pct)), cov

def gene_ok(row: pd.Series, gene: str, min_cov: float = MIN_GENE) -> Tuple[bool, Optional[float]]:
    cds = parse_cds_coverage(row.get('cdsCoverage'))
    v = cds.get(gene)
    # cdsCoverage is 0–1 in Nextclade; compare with min_cov% threshold (0.80 here)
    return (True if v is None else (v >= (min_cov/100.0))), v

def mask_nextclade_fields(row_dict: Dict[str, object], genes: List[str], which: str) -> None:
    """which='all' → mask all NC fields except Sample/coverage/cdsCoverage.
       which='<gene>' → mask only that gene's aa fields."""
    if which == 'all':
        for k in list(row_dict.keys()):
            if k in ('Sample','coverage','cdsCoverage'):
                continue
            row_dict[k] = 'NA'
        return
    gene = which
    patt = re.compile(fr'^{gene}_(aaSubstitutions|aaDeletions|aaInsertions)(?:_\d+)?$')
    for k in list(row_dict.keys()):
        if patt.fullmatch(k or ""):
            row_dict[k] = 'NA'

# ---------- load inhibitor tables ----------
clpro_df = pd.read_csv(clpro_csv)
rdrp_df  = pd.read_csv(rdrp_csv)
spike_df = pd.read_csv(spike_csv)

# Normalize deletions in tables, then compute Nextclade lookups
for d in (clpro_df, rdrp_df, spike_df):
    d['Mutation'] = d['Mutation'].map(lambda x: normalise_del(str(x)))

clpro_df['Nextclade_lookup'] = clpro_df['Mutation'].apply(lambda m: offset(m, 3263))  # nsp5 +3263
rdrp_df ['Nextclade_lookup'] = rdrp_df ['Mutation'].apply(lambda m: offset(m,  -9))  # nsp12 −9
spike_df['Nextclade_lookup'] = spike_df['Mutation']                                      # spike already in S coords

SPIKE_FOLD_COL = guess_fold_col(spike_df, prefer=['BAM','ADG','TIX','SOT','BEB','REGN'])
CLPRO_FOLD_COL = guess_fold_col(clpro_df,  prefer=['NTV','PAX','NIR'])
RDRP_FOLD_COL  = guess_fold_col(rdrp_df,   prefer=['RDV','REM'])
print(f"[{__VERSION__}] fold columns: Spike='{SPIKE_FOLD_COL}', 3CLpro='{CLPRO_FOLD_COL}', RdRp='{RDRP_FOLD_COL}'",
      file=sys.stderr)

# ---------- iterate samples ----------
main_df = pd.read_csv(main_csv)
if 'Sample' not in main_df.columns and 'seqName' in main_df.columns:
    main_df = main_df.rename(columns={'seqName': 'Sample'})

genes_of_interest = ['S','ORF1a','ORF1b']
rows: List[Dict[str, object]] = []
masked_all = 0

for _, r in main_df.iterrows():
    out = r.to_dict()  # start by passing through Nextclade row as-is
    # ensure Sample exists
    if 'Sample' not in out and 'seqName' in out: out['Sample'] = out['seqName']

    ok_all, cov = overall_ok(r)
    s_ok, s_c = gene_ok(r, 'S')
    a_ok, a_c = gene_ok(r, 'ORF1a')
    b_ok, b_c = gene_ok(r, 'ORF1b')

    # default resistance outputs
    out['Spike_mAbs_inhibitors'] = 'NA'
    out['Spike_Fold'] = 'NA'
    out['3CLpro_inhibitors'] = 'NA'
    out['3CLpro_Fold'] = 'NA'
    out['RdRp_inhibitors'] = 'NA'
    out['RdRp_Fold'] = 'NA'

    if not ok_all:
        masked_all += 1
        mask_nextclade_fields(out, genes_of_interest, 'all')
        out['QC_overall_pct'] = cov
        out['QC_cds_S'] = s_c
        out['QC_cds_ORF1a'] = a_c
        out['QC_cds_ORF1b'] = b_c
        out['QC'] = 'masked: low overall coverage'
        rows.append(out)
        continue

    # Per-gene NC masking + resistance lookups only for genes that pass cds gate
    # SPIKE
    if s_ok:
        s_nc   = gather_gene_muts(r, 'S')
        s_disp = to_lookup(s_nc, spike_df)
        s_fold = max_fold(s_nc, spike_df, SPIKE_FOLD_COL)
        out['Spike_mAbs_inhibitors'] = ",".join(s_disp) if s_disp else 'No Mutations'
        out['Spike_Fold'] = s_fold if s_fold is not None else 'No Data'
    else:
        mask_nextclade_fields(out, genes_of_interest, 'S')

    # 3CLpro (ORF1a → nsp5)
    if a_ok:
        a_nc   = gather_gene_muts(r, 'ORF1a')
        a_disp = to_lookup(a_nc, clpro_df)
        a_fold = max_fold(a_nc, clpro_df, CLPRO_FOLD_COL)
        out['3CLpro_inhibitors'] = ",".join(a_disp) if a_disp else 'No Mutations'
        out['3CLpro_Fold'] = a_fold if a_fold is not None else 'No Data'
    else:
        mask_nextclade_fields(out, genes_of_interest, 'ORF1a')

    # RdRp (ORF1b → nsp12)
    if b_ok:
        b_nc   = gather_gene_muts(r, 'ORF1b')
        b_disp = to_lookup(b_nc, rdrp_df)
        b_fold = max_fold(b_nc, rdrp_df, RDRP_FOLD_COL)
        out['RdRp_inhibitors'] = ",".join(b_disp) if b_disp else 'No Mutations'
        out['RdRp_Fold'] = b_fold if b_fold is not None else 'No Data'
    else:
        mask_nextclade_fields(out, genes_of_interest, 'ORF1b')

    out['QC_overall_pct'] = cov
    out['QC_cds_S'] = s_c
    out['QC_cds_ORF1a'] = a_c
    out['QC_cds_ORF1b'] = b_c
    out['QC'] = 'ok' if all([s_ok, a_ok, b_ok]) else '; '.join(
        ['low S cds' if not s_ok else '',
         'low ORF1a cds' if not a_ok else '',
         'low ORF1b cds' if not b_ok else '']).strip('; ').replace(';;',';')

    rows.append(out)

# ---------- write output ----------
out_df = pd.DataFrame(rows)
out_file = f'{run_id}_resistance_mutations.csv'
# keep_default_na=False so the literal string "NA" stays as "NA" (not NaN)
out_df.to_csv(out_file, index=False)
print(f"[{__VERSION__}] ✔ wrote {out_file} (hard-masked rows: {masked_all})", file=sys.stderr)
