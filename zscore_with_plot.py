"""
    test.py
    10/06/2020
    Calculating the z-scores for each sample.
"""

# Run on Spyder IDE

# Import Pandas package, used for data analysis
import pandas as pd
from sklearn.decomposition import PCA

file_name = '1-s2.0-S0092867420307443-mmc5.xlsx'
file = pd.ExcelFile(file_name)

# Parse Excel sheet into tableS5A DataFrame
tableS5A_df = file.parse('Table S5A', skiprows=[0, 1], index_col=0)

def z_score(df, row):
    # Input: DataFrame and specific row
    # Output: DataFrame row with all numbers as z-scores
    
    # z-score = (raw score - mean of row) / (std dev of row)
    return (df.loc[row] - df.loc[row].mean()) / df.loc[row].std()

# Create a new DataFrame for storing z-scores
tableS5A_zscores_df = tableS5A_df.copy()

tumor_normal = file.parse('TableS5B', skiprows=[0, 1])

# Loop through all rows in tableS5A DataFrame
for i, row in tableS5A_df.iterrows():
    
    # Convert rows into z-scores
    tableS5A_zscores_df.loc[i] = z_score(tableS5A_df, i)
    
#print(tableS5A_zscores_df)

# Export tableS5A z-scores to a new Excel file
tableS5A_zscores_df.to_excel("TableS5A_zscores.xlsx", sheet_name='Table S5A Z-Scores')

# Testing out PCA
zscores_array = tableS5A_zscores_df.to_numpy()
pca = PCA(n_components=2)
pca.fit(zscores_array)
components = pca.components_

# Import the plotting library matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

colors = []
for i, row in tumor_normal.iterrows():
    if '.N' in tumor_normal.loc[i, 'Sample ID']:
        colors.append('r')
    else:
        colors.append('b')


# Make a simple 2D scatter plot
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA of xCell Data')
plt.scatter(components[0], components[1], color=colors)
red = Line2D(range(1), range(1), color="white", marker='o', markersize=10,
             markerfacecolor="r")
blue = Line2D(range(1), range(1), color="white", marker='o', markersize=10,
             markerfacecolor="b")
plt.legend([red,blue],['normal','tumor'], loc='lower right')

plt.savefig('pca.png', bbox_inches='tight', dpi=300)

# Show the plot
plt.show()
