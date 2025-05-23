# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 13:35:40 2024

@author: osama
"""

import pandas as pd
import numpy as np
from datetime import date, timedelta
import os
#Directory where the Script loactes
os.chdir('C:/Work/Research/Data Analysis/Tools/Python_Scripts')
from Data_Analyses_Fns import *
os.chdir('C:/LOONE_WQ/Data/LORS19722023/ts_data')

# Specify the folder path
folder_path = 'C:/LOONE_WQ/Data/LORS19722023/ts_data'

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

# Iterate through each CSV file
for file in csv_files:
    file_path = os.path.join(folder_path, file)
    data = pd.read_csv(file_path)
    
    # Check if the date in the first row is '1/1/1972'
    if 'date' in data.columns:
        data['date'] = pd.to_datetime(data['date'])
        if data['date'].iloc[0] == pd.to_datetime('1/1/1972'):
            # Apply your custom function to change the date range
            data = DF_Date_Range(data, 2007, 1, 1, 2023, 6, 30)
        elif data['date'].iloc[0] == pd.to_datetime('12/30/1971'):
            # Apply your custom function to change the date range
            data = DF_Date_Range(data, 2006, 12, 30, 2023, 6, 30)
        else:
            pass
    else:
        pass
        # Save the modified DataFrame to a new CSV file
    data.to_csv('C:/LOONE_WQ/Data/LORS20072023/ts_data/%s_20072023.csv' % file, index=False)
    

def rename_files(folder_path, old_part, new_part):
    # Navigate to the folder containing the files
    os.chdir('C:/LOONE_WQ/Data/LORS20072023/ts_data')

    # Iterate through all files in the folder
    for filename in os.listdir(folder_path):
        # Check if the old_part exists in the filename
        if old_part in filename:
            # Create the new filename by replacing old_part with new_part
            new_filename = filename.replace(old_part, new_part)

            # Construct the full paths
            old_path = os.path.join(folder_path, filename)
            new_path = os.path.join(folder_path, new_filename)

            # Rename the file
            os.rename(old_path, new_path)
            print(f'Renamed: {filename} to {new_filename}')

# Example usage
folder_path = 'C:/LOONE_WQ/Data/LORS20072023/ts_data'
old_part = 'LORS20082023.csv_20072023'
new_part = 'LORS20072023'

rename_files(folder_path, old_part, new_part)