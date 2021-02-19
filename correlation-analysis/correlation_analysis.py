"""
    correlation_analysis.py
    11/10/2020
    @author: Bridget Li
    Create matrix of ρ values for each phosphopeptide and immune cell. Make a heatmap
    from the matrix to help interpret the data.
"""

import pandas as pd
import seaborn as sns
from scipy.stats.stats import spearmanr

# Import Excel files with phosphopeptide and immune cell data.
immune_file_name = '1-s2.0-S0092867420307443-mmc5.xlsx'
phospho_file_name = 'CPTAC_allPhosData_log2FCpooled.csv'

file = pd.ExcelFile(immune_file_name)
phospho_df = pd.read_csv(phospho_file_name)

# Parse Excel sheet into tableS5A DataFrame
tableS5A_df = file.parse('Table S5A', skiprows=[0, 1], index_col=0)
# Drop columns of normal adjacent tissue (NAT) samples (columns ending with '.N').
tableS5A_df = tableS5A_df[tableS5A_df.columns[~tableS5A_df.columns.str.endswith('.N')]]

# Make a new column called 'Seq-Mods' that includes both the annotated sequence
# and modifications on the peptide. This is so peptides with the same sequence but
# different modifications will be treated as unique.
phospho_df['Seq-Mods'] = phospho_df['Annotated Sequence'] + ' - ' \
+ phospho_df['Modifications'] + ' - ' + phospho_df['Gene-Sites']
phospho_df.set_index('Seq-Mods', inplace=True)

# Drop columns that are not sample IDs.
phospho_df.drop([col for col in phospho_df.columns if 'Abundance' not in col],
                axis=1,inplace=True)

# Drop the word 'Abundance' from the phospho columns names so they'll match the
# immune column names.
phospho_df.columns = phospho_df.columns.str.replace(' Abundance', '')
# Drop columns of NAT samples.
phospho_df = phospho_df[tableS5A_df.columns]


# =============================================================================
# Start by z-scoring immune cell and peptide data.
# z-scores describe distance from the mean.
# =============================================================================

# z-score given row of DataFrame (all row values are z-scores).
def z_score(df, row):
    
    # z-score = (raw score - mean of row) / (std dev of row)
    return (df.loc[row] - df.loc[row].mean()) / df.loc[row].std()

# Create a new DataFrame for storing immune z-scores.
immune_df = tableS5A_df.copy()

# Loop through all rows in tableS5A DataFrame.
for i, row in tableS5A_df.iterrows():
    # Convert rows into z-scores.
    immune_df.loc[i] = z_score(tableS5A_df, i)

# Create a new DataFrame for storing peptide z-scores.
phospho_zscores_df = phospho_df.copy()

# Loop through all rows in tableS5A DataFrame.
for i, row in phospho_df.iterrows():
    # Convert rows into z-scores.
    phospho_zscores_df.loc[i] = z_score(phospho_df, i)
    
phospho_df = phospho_zscores_df

phospho_df.to_csv('phospho_zscores.csv')
immune_df.to_csv('immune_zscores.csv')

# =============================================================================
# Create matrix of Spearman correlation coefficients for each phosphopeptide and
# immune cell. This measures the strength and direction of the association between
# the two (can be non-linear).
# =============================================================================

# Create new matrix to hold ρ values (rows = phosphopeptides, columns = immune cells)
combined_df = pd.DataFrame(index=list(phospho_df.index.values), columns=list(immune_df.index.values))

# Iterate through phospho_df.
for i, row_i in phospho_df.iterrows():

    # For each phosphopeptide, get list of z-scores for all sample IDs.
    phospho_data = row_i.tolist()
    
    # Iterate through immune_df.
    for j, row_j in immune_df.iterrows():
        
        # For each immune cell, get list of z-scores for all sample IDs.
        immune_data = row_j.tolist()
        
        # Calculate Spearman ρ for z-scored phosphopeptide and immune cell lists.
        # Fill combined_df with z-score for given phosphopeptide and immune cell.
        combined_df.loc[i,j] = spearmanr(phospho_data, immune_data)[0]

# Filter out peptides that don't have at least 1 strong correlation (ρ values all
# between -0.178 and 0.178, corresponding to p<0.01).
combined_df = combined_df.loc[((combined_df>0.246) | (combined_df<-0.246)).any(1)]

# Export combined_df to a new Excel file.
combined_df.to_csv("correlation_matrix.csv")


# =============================================================================
# Create a heatmap of ρ values for each phosphopeptide and immune cell. Reorder
# the correlation matrix using new row and column indices. This will help us
# understand the relationships between the peptides and immune cells.
# =============================================================================


# Create heatmap using seaborn library.
combined_df = pd.read_csv('correlation_matrix.csv', index_col=0)
sns.set(font_scale = 0.6)
g = sns.clustermap(combined_df, cmap='seismic', center=0, metric='correlation',
                   yticklabels=False, cbar_kws={'label': "Spearman's ρ"})
g.savefig("phospho_immune_heatmap.png")

# Get the new row (peptide-level) and column (sample-level) indices.
new_peptide_indices = g.dendrogram_row.reordered_ind
new_sample_indices = g.dendrogram_col.reordered_ind

# Reorder the rows to match the clustermap.
reordered_peptides = [list(combined_df.index.values)[i] for i in new_peptide_indices]
df_reordered = combined_df.reindex(reordered_peptides)

# Reorder the samples to match the clustermap.
reordered_samples = [list(df_reordered.columns)[i] for i in new_sample_indices]
df_reordered = df_reordered[reordered_samples]

# Export reordered correlation matrix with new indices.
df_reordered.to_csv('reordered_correlation_matrix.csv')
