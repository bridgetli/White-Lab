"""
    biomarkers.py
    02/15/2021
    @author: Bridget Li
    Identify which peptides found in the tumor samples are also known biomarkers
    for specific immune cells.
"""

import pandas as pd

# Import Excel files containing sample peptides and biomarker peptides
peptides_df = pd.read_csv('CPTAC_allPhosData_log2FCpooled.csv')
biomarker_df = pd.read_csv('immune_cell_markers_final.csv',encoding='latin1')

# Separate peptides by gene IDs
gene_IDs = peptides_df['Gene IDs'].tolist()
gene_IDs_separated = []
separated = []
for ID in gene_IDs:
    if '; ' in ID:
        separated += ID.split('; ')
    else:
        gene_IDs_separated.append(ID)
gene_IDs = gene_IDs_separated
gene_IDs += separated
# Delete duplicate gene IDs
gene_IDs = list(set(gene_IDs))
gene_IDs.sort()

biomarkers = biomarker_df['gene'].tolist()

# Get list of matching peptides (based on gene IDs)
both = list(set(gene_IDs) & set(biomarkers))
print(both)


