"""
zscore_phospho_data.py
11/04/2020
@author: Bridget Li
"""

# Run on Spyder IDE

# Import Pandas package, used for data analysis
import pandas as pd

file_name = 'fake phospho data.csv'
# Parse Excel CSV file into phospho_data DataFrame
phospho_data_df = pd.read_csv(file_name, index_col=0)

phospho_data_df.drop([col for col in phospho_data_df.columns if 'Abundance' not in col],
                     axis=1,inplace=True)


def z_score(df, row):
    # Input: DataFrame and specific row
    # Output: DataFrame row with all numbers as z-scores
    
    # z-score = (raw score - mean of row) / (std dev of row)
    return (df.loc[row] - df.loc[row].mean()) / df.loc[row].std()

# Create a new DataFrame for storing z-scores
phospho_zscores_df = phospho_data_df.copy()

# Loop through all rows in tableS5A DataFrame
for i, row in phospho_data_df.iterrows():
    
    # Convert rows into z-scores
    phospho_zscores_df.loc[i] = z_score(phospho_data_df, i)

phospho_zscores_df.to_excel("phospho_data_zscores.xlsx", sheet_name='Phospho Data Zscores')
