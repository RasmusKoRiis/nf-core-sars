import pandas as pd
import sys

in_csv = sys.argv[1]
run_id = sys.argv[2]

# ---- Config ----
GENE_COV_MIN = 0.80  # threshold for per-gene cdsCoverage

# ---- Load & normalize empties to NA ----
df = pd.read_csv(in_csv, sep=';')
df = df.replace(r'^\s*$', pd.NA, regex=True)

# Gene segments that may exist
gene_segments = ['E', 'M', 'N', 'ORF1a', 'ORF1b', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'S']

def _is_missing(x) -> bool:
    return pd.isna(x) or str(x).strip().upper() in {'', 'NA', 'N/A', 'NONE', 'NULL'}

def parse_cds_coverage(s: str) -> dict:
    """
    Parse cdsCoverage like 'E:1,M:1,...,ORF1a:0.9343,...,S:1'
    Returns floats in 0–1; tolerates %, or 0–100 values.
    """
    out = {}
    if _is_missing(s):
        return out
    for part in str(s).split(','):
        if ':' not in part:
            continue
        gene, val = part.split(':', 1)
        gene = gene.strip()
        try:
            v = float(str(val).strip().replace('%', ''))
        except ValueError:
            continue
        if v > 1.0:
            v = v / 100.0  # normalize percent to fraction
        out[gene] = v
    return out

def _parse_gene_list(raw: str, gene: str):
    """Return list of values after 'Gene:' for a comma-separated field."""
    if _is_missing(raw):
        return []
    toks = [t for t in str(raw).split(',') if ':' in t]
    vals = [t.split(':', 1)[1].strip() for t in toks if t.startswith(gene + ':')]
    return [v for v in vals if v]

def extract_mutations_by_gene(row: pd.Series, field: str, gene: str):
    """
    Apply gating:
      - if cdsCoverage[gene] < GENE_COV_MIN or missing -> NA
      - else return joined mutations for that gene, or "No mutation" if none
    field ∈ {'aaSubstitutions','aaDeletions','aaInsertions'}
    """
    covs = parse_cds_coverage(row.get('cdsCoverage'))
    cov = covs.get(gene, None)
    if cov is None or cov < GENE_COV_MIN:
        return pd.NA
    muts = _parse_gene_list(row.get(field), gene)
    return ",".join(muts) if muts else "No mutation"

# ---- Extract the stats subset (unchanged) ----
nextclade_stats = df[['seqName', 'clade', 'Nextclade_pango', 'partiallyAliased',
                      'clade_nextstrain', 'clade_who',
                      'qc.mixedSites.totalMixedSites', 'coverage',
                      'qc.mixedSites.totalMixedSites', 'qc.overallStatus','frameShifts']].copy()
nextclade_stats.columns = ['Sample', 'clade', 'Nextclade_pango', 'partiallyAliased',
                           'clade_nextstrain', 'clade_who',
                           'qc.mixedSites.totalMixedSites', 'coverage',
                           'NC_Genome_MixedSites', 'NC_Genome_QC', 'NC_Genome_frameShifts']

# ---- Extract mutations subset (keep coverage fields present) ----
nextclade_mutations = df[['seqName', 'aaSubstitutions', 'aaDeletions', 'aaInsertions',
                          'coverage', 'cdsCoverage']].copy()
nextclade_mutations.columns = ['Sample', 'aaSubstitutions', 'aaDeletions', 'aaInsertions',
                               'coverage', 'cdsCoverage']

# ---- Create gene-specific mutation columns with per-gene gating ----
for segment in gene_segments:
    nextclade_mutations[f'{segment}_aaSubstitutions'] = nextclade_mutations.apply(
        lambda row: extract_mutations_by_gene(row, 'aaSubstitutions', segment), axis=1)
    nextclade_mutations[f'{segment}_aaDeletions'] = nextclade_mutations.apply(
        lambda row: extract_mutations_by_gene(row, 'aaDeletions', segment), axis=1)
    nextclade_mutations[f'{segment}_aaInsertions'] = nextclade_mutations.apply(
        lambda row: extract_mutations_by_gene(row, 'aaInsertions', segment), axis=1)

# ---- Drop the original combined columns (you now have per-gene fields) ----
nextclade_mutations = nextclade_mutations.drop(columns=['aaSubstitutions', 'aaDeletions', 'aaInsertions'])

# ---- Split long S & ORF1a substitutions into 4 columns each, preserving NA ----
def split_aa_substitutions(substitutions):
    if _is_missing(substitutions):
        return [pd.NA, pd.NA, pd.NA, pd.NA]
    if substitutions == "No mutation":
        return ["No mutation", pd.NA, pd.NA, pd.NA]
    changes = str(substitutions).split(',')
    if len(changes) <= 22:
        return [",".join(changes), pd.NA, pd.NA, pd.NA]
    elif len(changes) <= 44:
        return [",".join(changes[:22]), ",".join(changes[22:]), pd.NA, pd.NA]
    elif len(changes) <= 66:
        return [",".join(changes[:22]), ",".join(changes[22:44]), ",".join(changes[44:]), pd.NA]
    else:
        return [",".join(changes[:22]), ",".join(changes[22:44]), ",".join(changes[44:66]), ",".join(changes[66:])]

nextclade_mutations[['S_aaSubstitutions_1', 'S_aaSubstitutions_2', 'S_aaSubstitutions_3', 'S_aaSubstitutions_4']] = \
    pd.DataFrame(nextclade_mutations['S_aaSubstitutions'].apply(split_aa_substitutions).tolist(),
                 index=nextclade_mutations.index)

nextclade_mutations[['ORF1a_aaSubstitutions_1', 'ORF1a_aaSubstitutions_2', 'ORF1a_aaSubstitutions_3', 'ORF1a_aaSubstitutions_4']] = \
    pd.DataFrame(nextclade_mutations['ORF1a_aaSubstitutions'].apply(split_aa_substitutions).tolist(),
                 index=nextclade_mutations.index)

# ---- Output names ----
mutation_name = run_id + '_nextclade_mutations.csv'
stat_name     = run_id + '_nextclade_stats.csv'

# ---- Final polish: replace inner commas (CSV hygiene), re-sweep blanks, write with NA ----
nextclade_stats     = nextclade_stats.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)
nextclade_mutations = nextclade_mutations.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

nextclade_stats     = nextclade_stats.replace(r'^\s*$', pd.NA, regex=True)
nextclade_mutations = nextclade_mutations.replace(r'^\s*$', pd.NA, regex=True)

nextclade_stats.to_csv(stat_name, index=False, na_rep='NA')
nextclade_mutations.to_csv(mutation_name, index=False, na_rep='NA')
