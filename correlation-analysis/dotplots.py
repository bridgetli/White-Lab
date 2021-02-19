"""
    dotplots.py
    01/17/2020
    @author: Bridget Li
    From matrix of ρ values, pick top 10 highest and top 10 lowest values. Create
    dotplots to show correlation. Color by hot vs. cold tumor.
"""

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy
import os
import shutil

# Import Excel file with reordered correlation matrix
corr_matrix = pd.read_csv('reordered_correlation_matrix.csv', index_col=0)

top_bottom_10 = {}

# Find 10 highest ρ values in matrix and add to list
while len(top_bottom_10) < 10:
    highest_rho = 0
    for i, row_i in corr_matrix.iterrows():
        for r in row_i:
            if r > highest_rho and r not in top_bottom_10:
                highest_r = r
                phospho = row_i.name
                immune = (row_i[row_i==r].index).format()[0]
            
    top_bottom_10[highest_r] = (phospho, immune)

# Find 10 lowest ρ values in matrix and add to list
while len(top_bottom_10) < 20:
    lowest_r = 0
    phospho = ""
    immune = ""
    for i, row_i in corr_matrix.iterrows():
        for r in row_i:
            if r < lowest_r and r not in top_bottom_10:
                lowest_r = r
                phospho = row_i.name
                immune = (row_i[row_i==r].index).format()[0]
            
    top_bottom_10[lowest_r] = (phospho, immune)

# Import Excel file that specifies whether each tumor is hot or cold
file = pd.ExcelFile('1-s2.0-S0092867420307443-mmc5.xlsx')
hot_cold_df = file.parse('TableS5B', skiprows=[0, 1])

# Add specified colors (red for hot, blue for cold) for each sample to color list
colors = []
for i, row in hot_cold_df.iterrows():
    if '.N' not in hot_cold_df.loc[i, 'Sample ID']:
        if 'Hot' in hot_cold_df.loc[i, 'Group']:
            colors.append('r')
        elif 'Cold' in hot_cold_df.loc[i, 'Group']:
            colors.append('b')
        else:
            colors.append('k')

# Import Excel files with phosphopeptide and immune cell z-scores
phospho_df = pd.read_csv('phospho_zscores.csv', index_col=0)
immune_df = pd.read_csv('immune_zscores.csv', index_col=0)

# Create folder to contain dotplot figures
if os.path.exists('dotplots'):
    shutil.rmtree('dotplots')
os.mkdir('dotplots')

# For each top or bottom pair (peptide + immune) in list
for i, r in enumerate(top_bottom_10):
    # Get peptide and immune cell z-scores for each tumor sample
    phosphopeptide = top_bottom_10.get(r)[0]
    immune_cell = top_bottom_10.get(r)[1]
    phospho_data = phospho_df.loc[phosphopeptide].values.tolist()
    immune_data = immune_df.loc[immune_cell].values.tolist()
    # Create new figure
    plt.figure()
    # Add trendline
    z = numpy.polyfit(phospho_data, immune_data, 1)
    p = numpy.poly1d(z)
    plt.plot(phospho_data,p(phospho_data),"r--", color='k')
    plt.title("ρ = {:.4f}".format(r))   # Add ρ value as title
    # Add dots
    plt.scatter(phospho_data, immune_data, color=colors)
    
    # Custom legend
    red = Line2D(range(1), range(1), color="white", marker='o', markersize=10,
             markerfacecolor="r")
    blue = Line2D(range(1), range(1), color="white", marker='o', markersize=10,
             markerfacecolor="b")
    black = Line2D(range(1), range(1), color="white", marker='o', markersize=10,
             markerfacecolor="k")
    plt.legend([red,blue,black],['hot','cold', 'NAT'], loc='best')
    
    # Add axis labels
    plt.xlabel(phosphopeptide)
    plt.ylabel(immune_cell)
    
    # Save figures with rank
    if i <= 9:
        plt.savefig('dotplots\\top' + str(i+1) + '.png', bbox_inches='tight', dpi=300)
    else:
        plt.savefig('dotplots\\bottom' + str(i-9) + '.png', bbox_inches='tight', dpi=300)

plt.show()
