#!/usr/bin/env python3
"""
Generate a resistance-mutation summary from a Nextclade CSV with per-gene QC.
Call:
python resistance_check.py  <nextclade.csv>  <run_id>  <spike.csv>  <rdrp.csv>  <clpro.csv>

Env overrides (optional):
  MIN_OVERALL_PCT   default 90
  MIN_GENE_CDS      default 90
"""
import sys, re, os, pandas as pd
from typing import Optional, List, Dict

__VERSION__ = "resistance_check v2025-08-21"

if len(sys.argv) < 6:
    sys.exit(f"Usage: {sys.argv[0]} <nextclade.csv> <run_id> <spike.csv> <rdrp.csv> <clpro.csv>")

main_csv, run_id, spike_csv, rdrp_csv, clpro_csv = sys.argv[1:6]
MIN_OVERALL = float(os.getenv("MIN_OVERALL_PCT", "90"))
MIN_GENE = float(os.getenv("MIN_GENE_CDS", "90"))   # was 95; 90 is more permissive for ONT

print(f"[{__VERSION__}] main={main_csv} run={run_id} gates: overall>={MIN_OVERALL}%, gene>={MIN_GENE}", file=sys.stderr)

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
        if not b: continue
        out.append(normalise_del(b) if dels else b)
    return out

def gather_gene_muts(row: pd.Series, gene: str) -> List[str]:
    """Collect per-gene AA changes, incl. split columns _1.._4; return list WITHOUT gene prefix."""
    muts: List[str] = []
    for col in row.index:
        c = str(col)
        if re.fullmatch(fr'{gene}_aaSubstitutions(?:_\d+)?', c):
            muts.extend(split_clean(row.get(col, '')))
        elif re.fullmatch(fr'{gene}_aaDeletions(?:_\d+)?', c):
            muts.extend([normalise_del(x) for x in split_clean(row.get(col, ''), dels=True)])
        elif re.fullmatch(fr'{gene}_aaInsertions(?:_\d+)?', c):
            muts.extend(split_clean(row.get(col, '')))
    # Deduplicate preserving order
    seen = set(); uniq: List[str] = []
    for m in muts:
        if m not in seen:
            uniq.append(m); seen.add(m)
    # Add alt deletion forms to maximize matches
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

def overall_ok(row, min_pct=MIN_OVERALL) -> (bool, Optional[float]):
    cov = parse_percentish(row.get('coverage'))
    return (True if cov is None else (cov >= min_pct)), cov

def gene_ok(row, gene: str, min_cov: float = MIN_GENE) -> (bool, Optional[float]):
    cds = parse_cds_coverage(row.get('cdsCoverage'))
    v = cds.get(gene)
    return (True if v is None else (v >= min_cov)), v

# ---------- load look-up tables ----------
clpro_df = pd.read_csv(clpro_csv)
rdrp_df  = pd.read_csv(rdrp_csv)
spike_df = pd.read_csv(spike_csv)

# Normalize both sides for deletions before computing lookups
clpro_df['Mutation'] = clpro_df['Mutation'].map(lambda x: normalise_del(str(x)))
rdrp_df ['Mutation'] = rdrp_df ['Mutation'].map(lambda x: normalise_del(str(x)))
spike_df['Mutation'] = spike_df['Mutation'].map(lambda x: normalise_del(str(x)))

# Compute Nextclade_lookup keys (nsp→ORF offsets)
clpro_df['Nextclade_lookup'] = clpro_df['Mutation'].apply(lambda m: offset(m, 3263))  # nsp5 +3263
rdrp_df ['Nextclade_lookup'] = rdrp_df ['Mutation'].apply(lambda m: offset(m,  -9))  # nsp12 −9
spike_df['Nextclade_lookup'] = spike_df['Mutation']  # spike already in S space

# Choose fold columns robustly
SPIKE_FOLD_COL = guess_fold_col(spike_df, prefer=['BAM','ADG','TIX','SOT','BEB','REGN'])
CLPRO_FOLD_COL = guess_fold_col(clpro_df,  prefer=['NTV','PAX','NIR'])  # Nirmatrelvir/Paxlovid
RDRP_FOLD_COL  = guess_fold_col(rdrp_df,   prefer=['RDV','REM'])        # Remdesivir
print(f"[{__VERSION__}] fold columns: Spike='{SPIKE_FOLD_COL}', 3CLpro='{CLPRO_FOLD_COL}', RdRp='{RDRP_FOLD_COL}'", file=sys.stderr)

# ---------- iterate samples ----------
main_df = pd.read_csv(main_csv)
if 'Sample' not in main_df.columns and 'seqName' in main_df.columns:
    main_df = main_df.rename(columns={'seqName': 'Sample'})

rows = []
masked_all = 0

for _, r in main_df.iterrows():
    ok_all, cov = overall_ok(r)
    s_ok, s_c = gene_ok(r, 'S')
    a_ok, a_c = gene_ok(r, 'ORF1a')
    b_ok, b_c = gene_ok(r, 'ORF1b')

    if not ok_all:
        masked_all += 1
        rows.append({
            'Sample': r.get('Sample', ''),
            'Spike_mAbs_inhibitors': 'NA', 'Spike_Fold': 'NA',
            'RdRp_inhibitors': 'NA',     'RdRp_Fold':  'NA',
            '3CLpro_inhibitors': 'NA',   '3CLpro_Fold':'NA',
            'QC_overall_pct': cov, 'QC_cds_S': s_c, 'QC_cds_ORF1a': a_c, 'QC_cds_ORF1b': b_c,
            'QC': 'masked: low overall coverage'
        })
        continue

    # SPIKE
    if s_ok:
        s_nc   = gather_gene_muts(r, 'S')
        s_disp = to_lookup(s_nc, spike_df)
        s_fold = max_fold(s_nc, spike_df, SPIKE_FOLD_COL)
        spike_inhib  = ",".join(s_disp) if s_disp else 'No Mutations'
        spike_fold   = s_fold if s_fold is not None else 'No Data'
    else:
        spike_inhib, spike_fold = 'NA', 'NA'

    # 3CLpro (ORF1a → nsp5)
    if a_ok:
        a_nc   = gather_gene_muts(r, 'ORF1a')
        a_disp = to_lookup(a_nc, clpro_df)
        a_fold = max_fold(a_nc, clpro_df, CLPRO_FOLD_COL)
        clpro_inhib = ",".join(a_disp) if a_disp else 'No Mutations'
        clpro_fold  = a_fold if a_fold is not None else 'No Data'
    else:
        clpro_inhib, clpro_fold = 'NA', 'NA'

    # RdRp (ORF1b → nsp12)
    if b_ok:
        b_nc   = gather_gene_muts(r, 'ORF1b')
        b_disp = to_lookup(b_nc, rdrp_df)
        b_fold = max_fold(b_nc, rdrp_df, RDRP_FOLD_COL)
        rdrp_inhib = ",".join(b_disp) if b_disp else 'No Mutations'
        rdrp_fold  = b_fold if b_fold is not None else 'No Data'
    else:
        rdrp_inhib, rdrp_fold = 'NA', 'NA'

    rows.append({
        'Sample': r.get('Sample', ''),
        'Spike_mAbs_inhibitors': spike_inhib,
        'Spike_Fold':  spike_fold,
        'RdRp_inhibitors': rdrp_inhib,
        'RdRp_Fold':   rdrp_fold,
        '3CLpro_inhibitors': clpro_inhib,
        '3CLpro_Fold': clpro_fold,
        'QC_overall_pct': cov,
        'QC_cds_S': s_c, 'QC_cds_ORF1a': a_c, 'QC_cds_ORF1b': b_c,
        'QC': ('ok' if (s_ok or spike_inhib=='NA') and (a_ok or clpro_inhib=='NA') and (b_ok or rdrp_inhib=='NA') else 'check')
    })

# ---------- write output ----------
out_df = pd.DataFrame(rows)
out_file = f'{run_id}_resistance_mutations.csv'
out_df.to_csv(out_file, index=False)
print(f"[{__VERSION__}] ✔ wrote {out_file} (hard-masked rows: {masked_all})", file=sys.stderr)
