import pandas as pd
import sys
import re
import math

csv = sys.argv[1]
id = sys.argv[2]
Spike = sys.argv[3]
RdRp = sys.argv[4]
CLpro = sys.argv[5]


# --- helper --------------
def split_mutations(raw):
    if pd.isna(raw) or raw.strip().lower() == 'no mutation':
        return []
    return [m.strip() for m in re.split(r'[;,]', raw) if m.strip()]

def normalise_deletion(m):
    """
    Turn 'L24-' → 'L24del'  (handles any AA letter + position + '-')
    Leaves substitutions untouched.
    """
    m = m.strip()
    if m.endswith('-'):
        return m[:-1] + 'del'
    return m

import re

def split_and_clean(raw, *, is_deletion=False):
    """
    Split on ',' or ';', drop blanks / 'No mutation', strip whitespace,
    optionally convert deletions to 'del' form.
    """
    if pd.isna(raw) or raw.strip().lower() == 'no mutation':
        return []
    parts = [p.strip() for p in re.split(r'[;,]', raw) if p.strip()]
    if is_deletion:
        parts = [normalise_deletion(p) for p in parts]
    return parts


# Function to adjust mutation positions
def adjust_mutation_position(mutation_str, offset):
    match = re.match(r'^([A-Z])(\d+)([A-Z]|del)$', mutation_str)
    if match:
        aa_from, position, aa_to = match.groups()
        new_position = int(position) + offset
        return f"{aa_from}{new_position}{aa_to}"
    else:
        return mutation_str

def adjust_mutations_list(mutations_str, offset):
    if pd.isna(mutations_str) or mutations_str == 'No mutation':
        return []
    mutations = [m.strip() for m in mutations_str.split(',') if m.strip()]
    adjusted_mutations = []
    for m in mutations:
        adjusted_m = adjust_mutation_position(m, offset)
        adjusted_mutations.append(adjusted_m)
    return adjusted_mutations

def get_inhibitors_from_list(mutations_list, lookup_df):
    found_mutations = []
    for m in mutations_list:
        if m in lookup_df['Nextclade_lookup'].values:
            found_mutations.append(m)
    return found_mutations

def get_max_fold_change(mutations_list, lookup_df, fold_column):
    fold_changes = []
    for m in mutations_list:
        df_match = lookup_df[lookup_df['Nextclade_lookup'] == m]
        if not df_match.empty:
            fold_values = df_match[fold_column].values
            for fold_value in fold_values:
                fold_value_str = str(fold_value).replace(",", "").replace("%", "").strip()
                try:
                    fold_value = float(fold_value_str)
                    fold_changes.append(fold_value)
                except (ValueError, TypeError):
                    fold_changes.append(0)
                    pass
    
    fold_changes_filtered = [fc for fc in fold_changes if not math.isnan(fc)]
    
    if fold_changes_filtered:
        return max(fold_changes_filtered)
    else:
        return None

# Read main.csv
main_df = pd.read_csv(csv)

# Read lookup CSVs
clpro_df = pd.read_csv(CLpro)
rdrp_df = pd.read_csv(RdRp)
spike_df = pd.read_csv(Spike)

# Process lookup DataFrames
clpro_df['Nextclade_lookup'] = clpro_df['Mutation'].apply(lambda x: adjust_mutation_position(x, 3263))
rdrp_df['Nextclade_lookup'] = rdrp_df['Mutation'].apply(lambda x: adjust_mutation_position(x, 9))
spike_df['Nextclade_lookup'] = spike_df['Mutation']

# Initialize lists to store results
samples = []
spike_inhibitors = []
rdrp_inhibitors = []
clpro_inhibitors = []
spike_folds = []
rdrp_folds = []
clpro_folds = []

for idx, row in main_df.iterrows():
    sample = row['Sample']
    s_mutations = row.get('S_aaSubstitutions', '')
    orf1a_mutations = row.get('ORF1a_aaSubstitutions', '')
    orf1b_mutations = row.get('ORF1b_aaSubstitutions', '')
    
    # For Spike mutations
    s_subs  = split_and_clean(row.get('S_aaSubstitutions', ''))
    s_dels  = split_and_clean(row.get('S_aaDeletions', ''), is_deletion=True)
    s_mutations_list = s_subs + s_dels
    spike_found = get_inhibitors_from_list(s_mutations_list, spike_df)
    
    # For ORF1a mutations (3CLpro)
    orf1a_subs = split_and_clean(row.get('ORF1a_aaSubstitutions', ''))
    orf1a_dels = split_and_clean(row.get('ORF1a_aaDeletions', ''), is_deletion=True)
    orf1a_all  = orf1a_subs + orf1a_dels
    adjusted_orf1a_mutations = [adjust_mutation_position(m, 3263) for m in orf1a_all]
    clpro_found = get_inhibitors_from_list(adjusted_orf1a_mutations, clpro_df)
    
    # For ORF1b mutations (RdRp)
    orf1b_subs = split_and_clean(row.get('ORF1b_aaSubstitutions', ''))
    orf1b_dels = split_and_clean(row.get('ORF1b_aaDeletions', ''), is_deletion=True)
    orf1b_all  = orf1b_subs + orf1b_dels
    adjusted_orf1b_mutations = [adjust_mutation_position(m, 9) for m in orf1b_all]
    rdrp_found = get_inhibitors_from_list(adjusted_orf1b_mutations, rdrp_df)
    
    # Get max fold changes
    spike_max_fold = get_max_fold_change(spike_found, spike_df, 'BAM: fold')
    rdrp_max_fold = get_max_fold_change(rdrp_found, rdrp_df, 'RDV: fold')
    clpro_max_fold = get_max_fold_change(clpro_found, clpro_df, 'NTV: fold')
    
    # If no inhibitors found, add "No Mutations" otherwise join the list
    spike_inhibitors.append(','.join(spike_found) if spike_found else "No Mutations")
    rdrp_inhibitors.append(','.join(rdrp_found) if rdrp_found else "No Mutations")
    clpro_inhibitors.append(','.join(clpro_found) if clpro_found else "No Mutations")
    
    # If no fold changes, add "No Data"
    spike_folds.append(spike_max_fold if spike_max_fold is not None else "No Data")
    rdrp_folds.append(rdrp_max_fold if rdrp_max_fold is not None else "No Data")
    clpro_folds.append(clpro_max_fold if clpro_max_fold is not None else "No Data")
    
    samples.append(sample)

# Create output DataFrame
output_df = pd.DataFrame({
    'Sample': samples,
    'Spike_mAbs_inhibitors': spike_inhibitors,
    'Spike_Fold': spike_folds,
    'RdRp_inhibitors': rdrp_inhibitors,
    'RdRp_Fold': rdrp_folds,
    '3CLpro_inhibitors': clpro_inhibitors,
    '3CLpro_Fold': clpro_folds
})

# Write to CSV
resitance_csv__name = id + '_' + "resistance_mutations.csv"
output_df.fillna("No Mutations", inplace=True)
output_df.to_csv(resitance_csv__name, index=False)

print(f"Processing complete. Results saved to '{resitance_csv__name}'.")