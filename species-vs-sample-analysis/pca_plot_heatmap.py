"""
    pca_plot_heatmap.py
    10/06/2020
    Calculate the z-scores for each sample (immune cell and peptide).
    Run PCA on samples and plot to show patterns. Create heatmaps.
"""

# Run on Spyder IDE

# Import Pandas package for data analysis
import pandas as pd
from sklearn.decomposition import PCA
# Import the plotting library matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# Import seaborn for making heatmap
import seaborn as sns


# Import Excel files with phosphopeptide and immune cell data.
immune_file_name = '1-s2.0-S0092867420307443-mmc5.xlsx'
phospho_file_name = 'CPTAC_allPhosData_log2FCpooled.csv'

file = pd.ExcelFile(immune_file_name)
phospho_df = pd.read_csv(phospho_file_name, index_col=0)

# Parse Excel sheet into tableS5A DataFrame
tableS5A_df = file.parse('Table S5A', skiprows=[0, 1], index_col=0)

# Drop columns that are not sample IDs.
phospho_df.drop([col for col in phospho_df.columns if 'Abundance' not in col],
                axis=1,inplace=True)

# =============================================================================
# Start by z-scoring immune cell and peptide data.
# =============================================================================

# z-score given row of DataFrame (all row values are z-scores).
def z_score(df, row):
    
    # z-score = (raw score - mean of row) / (std dev of row)
    return (df.loc[row] - df.loc[row].mean()) / df.loc[row].std()

# Create a new DataFrame for storing immune z-scores.
immune_zscores = tableS5A_df.copy()

# Loop through all rows in tableS5A DataFrame.
for i, row in tableS5A_df.iterrows():
    # Convert rows into z-scores.
    immune_zscores.loc[i] = z_score(tableS5A_df, i)

# Create a new DataFrame for storing peptide z-scores.
phospho_zscores = phospho_df.copy()

# Loop through all rows in tableS5A DataFrame.
for i, row in phospho_df.iterrows():
    # Convert rows into z-scores.
    phospho_zscores.loc[i] = z_score(phospho_df, i)


# =============================================================================
# Create PCA plots.
# =============================================================================
    
# Get colors for dots based on whether sample is normal or tumor.
tumor_normal = file.parse('TableS5B', skiprows=[0, 1])
colors = []
for i, row in tumor_normal.iterrows():
    if '.N' in tumor_normal.loc[i, 'Sample ID']:
        # Blue dots for normal samples.
        colors.append('b')
    else:
        # Red dots for tumor samples.
        colors.append('r')

# Function to create PCA plot.
def pca_plot(df, species):
    # 2 components for PCA.
    pca = PCA(n_components=2)
    pca.fit(df.to_numpy())
    components = pca.components_
    
    # Make a simple 2D scatter plot.
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA of ' + species + ' data ')
    plt.scatter(components[0], components[1], color=colors)
    
    # Create custom legend for plot.
    red = Line2D(range(1), range(1), color="white", marker='o', markersize=10,
                 markerfacecolor="r")
    blue = Line2D(range(1), range(1), color="white", marker='o', markersize=10,
                 markerfacecolor="b")
    plt.legend([red,blue],['tumor','normal'], loc='best')
    
    # Export figure to PNG.
    plt.savefig(species + '_pca.png', bbox_inches='tight', dpi=300)

    # Show the plot.
    plt.show()
    plt.clf()

sns.reset_orig()
# Make plots for immune cells and phosphopeptides.
pca_plot(immune_zscores, 'immune')
pca_plot(phospho_zscores, 'phospho')


# =============================================================================
# Create heatmaps for immune cells vs. samples and phosphopeptides vs. samples.
# =============================================================================

sns.set(font_scale = 0.5)
g = sns.clustermap(immune_zscores, cmap='seismic', center=0, metric='correlation',
                   cbar_kws={'label': "Z-score"})
g.savefig("immune_heatmap.png")

g = sns.clustermap(phospho_zscores, cmap='seismic', center=0, metric='correlation',
                   yticklabels=False, cbar_kws={'label': "Z-score"})
g.savefig("phospho_heatmap.png")
