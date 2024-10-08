import pandas as pd
import os
import sys

csv = sys.argv[1]
id = sys.argv[2]

# Load the CSV file into a DataFrame with comma delimiter
df = pd.read_csv(csv, sep=';')

# Gene segments that may exist
gene_segments = ['E', 'M', 'N', 'ORF1a', 'ORF1b', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'S']

# Function to extract mutations for each gene segment
def extract_mutations(mutation_str, gene_segment):
    if pd.isna(mutation_str):
        return "No mutation"
    mutations = [mut.split(":")[1] for mut in mutation_str.split(",") if mut.startswith(gene_segment + ":")]
    return ",".join(mutations) if mutations else "No mutation"


# Extract the first smaller CSV columns
nextclade_stats = df[['seqName', 'clade', 'Nextclade_pango', 'partiallyAliased', 'clade_nextstrain', 'clade_who','qc.mixedSites.totalMixedSites']]
nextclade_stats.columns = ['Sample', 'clade', 'Nextclade_pango', 'partiallyAliased', 'clade_nextstrain', 'clade_who','qc.mixedSites.totalMixedSites']

# Extract the second smaller CSV columns for amino acid changes
nextclade_mutations = df[['seqName', 'aaSubstitutions', 'aaDeletions', 'aaInsertions']]
nextclade_mutations.columns = ['Sample', 'aaSubstitutions', 'aaDeletions', 'aaInsertions']

# Create new columns for each gene segment and type (Substitutions, Deletions, Insertions)
for segment in gene_segments:
    nextclade_mutations[f'{segment}_aaSubstitutions'] = nextclade_mutations['aaSubstitutions'].apply(lambda x: extract_mutations(x, segment))
    nextclade_mutations[f'{segment}_aaDeletions'] = nextclade_mutations['aaDeletions'].apply(lambda x: extract_mutations(x, segment))
    nextclade_mutations[f'{segment}_aaInsertions'] = nextclade_mutations['aaInsertions'].apply(lambda x: extract_mutations(x, segment))

# Drop the original columns as they're now split into individual gene segments
nextclade_mutations = nextclade_mutations.drop(columns=['aaSubstitutions', 'aaDeletions', 'aaInsertions'])

# Function to split 'S_aaSubstitutions' into 4 parts according to the custom rules
def split_s_aa_substitutions(substitutions):
    if pd.isna(substitutions) or substitutions == "No mutation":
        return ["No mutation"] * 4
    changes = substitutions.split(',')
    if len(changes) <= 22:
        return [",".join(changes), "", "", ""]
    elif len(changes) <= 44:
        return [",".join(changes[:22]), ",".join(changes[22:]), "", ""]
    elif len(changes) <= 66:
        return [",".join(changes[:22]), ",".join(changes[22:44]), ",".join(changes[44:]), ""]
    else:
        return [",".join(changes[:22]), ",".join(changes[22:44]), ",".join(changes[44:66]), ",".join(changes[66:])]

# Modify the creation of S_aaSubstitutions to split it into 4 columns
nextclade_mutations[['S_aaSubstitutions_1', 'S_aaSubstitutions_2', 'S_aaSubstitutions_3', 'S_aaSubstitutions_4']] = \
    pd.DataFrame(nextclade_mutations['S_aaSubstitutions'].apply(split_s_aa_substitutions).tolist(), index=nextclade_mutations.index)

# Now drop the original S_aaSubstitutions column
#nextclade_mutations = nextclade_mutations.drop(columns=['S_aaSubstitutions'])

# Continue with saving the CSV as before
mutation_name = id + '_nextclade_mutations.csv'
stat_name = id + '_' + "nextclade_stats.csv"

nextclade_stats.to_csv(stat_name, index=False)
nextclade_mutations.to_csv(mutation_name, index=False)

