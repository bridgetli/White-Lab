# -*- coding: utf-8 -*-
"""
    relabel_with_ID.py
    10/25/2020
    Relabel reporter ion channels to corresponding sample IDs for each phosphoproteome
    data file.
"""

import pandas as pd
import csv
import os

# Parse Excel sheet into phospho_data_df DataFrame
phospho_data_df = pd.read_excel('S046_S056_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r2_July2020.xlsx')

# Create folder to store relabeled data files
if not os.path.exists('phospho_relabeled'):
    os.makedirs('phospho_relabeled')


def format_file_name(file_name):
    # Format file_name to match name of data file
    return file_name[:14]+file_name[28:]+'_BD_f_PSMs.txt'

def start_new_file(i):
    
    # Get file name, import DataFrame from tab-delimited text file, get headers of DataFrame
    file_name = format_file_name(phospho_data_df.loc[i, 'Phosphoproteome.Folder'])
    data_df = pd.read_csv('CPTAC_phospho_all\\' + file_name, sep='\t', header=(0))
    headers = list(data_df.columns.values)
    
    return file_name, data_df, headers

def export_prev_file(headers, data_df, file_name):
    
        # Find 'Abundance 131' in list of headers and rename to 'pooled Abundance'
        n = headers.index('Abundance 131')
        headers[n] = 'pooled Abundance'

        # Reset data_df with relabeled headers
        data_df.columns = headers
        # Export previous file to tab-delimited text file in new folder
        data_df.to_csv('phospho_relabeled\\' + file_name, index=None,
                       sep='\t', quoting=csv.QUOTE_ALL)
    
    
# Instantiate data file info for first file
file_name, data_df, headers = start_new_file(0)

# Loop through all rows in phospho_data_df
for i, row in phospho_data_df.iterrows():
    
    # When we move on to the next file
    if format_file_name(phospho_data_df.loc[i, 'Phosphoproteome.Folder']) != file_name:
        
        # Export previous file since it's finished
        export_prev_file(headers, data_df, file_name)
            
        # Instantiate data file info for next file
        file_name, data_df, headers = start_new_file(i)
        
    # Get reporter ion channel and sample ID from phospho_data_df
    channel = phospho_data_df.loc[i, 'Channel']
    broad_sample_ID = phospho_data_df.loc[i, 'Broad Sample.ID']
    
    # Find channel in list of headers and rename to sample ID
    n = headers.index('Abundance ' + str(channel))
    headers[n] = broad_sample_ID + ' Abundance'

# Export final data file
export_prev_file(headers, data_df, file_name)
