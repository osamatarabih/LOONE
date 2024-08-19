# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 23:28:26 2024

@author: osamatarabih
"""

# Lake Okeechobee TP bias corrected for different stations in the Lake proper

# Import required Packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# Identify data location
os.chdir('C:/Work/LOONE_P_BC/LORS20082023')
# Read simulated TP concentrations in Lake Okeechobee
Sim_TP = pd.read_csv('./LOONE_P_20082023_daily.csv')
Sim_TP['Date'] = pd.to_datetime(Sim_TP['Date'])
#Create a dataframe for the north TP simulations
TP_North_df = pd.DataFrame()
TP_North_df['date'] = Sim_TP['Date'] 
TP_North_df['TP_Lake_N'] = Sim_TP['TP_Lake_N'] 
TP_North_df = TP_North_df.set_index(['date'])
TP_North_df.index = pd.to_datetime(TP_North_df.index, unit = 'ns')
TP_North_Monthly = TP_North_df.resample('M').mean()

# Create dataframe for the south TP simulations
TP_South_df = pd.DataFrame()
TP_South_df['date'] = Sim_TP['Date'] 
TP_South_df['TP_Lake_S'] = Sim_TP['TP_Lake_S'] 
TP_South_df = TP_South_df.set_index(['date'])
TP_South_df.index = pd.to_datetime(TP_South_df.index, unit = 'ns')
TP_South_Monthly = TP_South_df.resample('M').mean()

# Station L001 (North of the Lake)
#read observations at L001
L001_Obs = pd.read_csv('./water_quality_L001_PHOSPHATE, TOTAL AS P.csv')
L001_Obs['date'] = pd.to_datetime(L001_Obs['date'])
L001_df = pd.DataFrame()
L001_df['date'] = L001_Obs['date']
L001_df['L001_Obs_mg/m3'] = L001_Obs['L001_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
L001_df = L001_df.set_index(['date'])
L001_df.index = pd.to_datetime(L001_df.index, unit = 'ns')
L001_df_Monthly = L001_df.resample('M').mean()

L001_biascorrect = pd.merge(TP_North_Monthly,L001_df_Monthly, how='left', on='date')
L001_biascorrect_noNan = L001_biascorrect.dropna()
Obs_L001 = L001_biascorrect_noNan['L001_Obs_mg/m3'].values
Sim_N_L001 = L001_biascorrect_noNan['TP_Lake_N'].values

# determine R2 for simulations related to observations at L001
slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_N_L001, Obs_L001)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_L001_N:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_L001, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_L001 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_L001 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_L001 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_L001 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_N = Sim_N_L001.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_N
correction_factor_2 = average_below_50th / average_TP_Lake_N
correction_factor_3 = average_below_75th / average_TP_Lake_N
correction_factor_4 = average_above_75th / average_TP_Lake_N

sim_BC_Qnt_EntPer = np.zeros(len(Obs_L001))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_L001)):
    if Obs_L001[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L001[i] * correction_factor_1 
    elif Obs_L001[i] > obs_quantiles[0] and Obs_L001[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L001[i] * correction_factor_2 
    elif Obs_L001[i] > obs_quantiles[1] and Obs_L001[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L001[i] * correction_factor_3
    elif Obs_L001[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L001[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_N_L001[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_L001)
r_squared = r_value**2
print("R2_L001_N_BC_Qnt:", r_squared)

# Station L005 north of the lake
#read observations at L005
L005_Obs = pd.read_csv('./water_quality_L005_PHOSPHATE, TOTAL AS P.csv')
L005_Obs['date'] = pd.to_datetime(L005_Obs['date'])
L005_df = pd.DataFrame()
L005_df['date'] = L005_Obs['date']
L005_df['L005_Obs_mg/m3'] = L005_Obs['L005_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
L005_df = L005_df.set_index(['date'])
L005_df.index = pd.to_datetime(L005_df.index, unit = 'ns')
L005_df_Monthly = L005_df.resample('M').mean()
L005_df_Monthly = L005_df_Monthly.reset_index()
#Create a dataframe of both Simulations and Observations
L005_biascorrect = pd.merge(TP_North_Monthly,L005_df_Monthly, how='left', on='date')
L005_biascorrect_noNan = L005_biascorrect.dropna()
L005_biascorrect_noNan.to_csv('./L005_biascorrect_noNan.csv')
Obs_L005 = L005_biascorrect_noNan['L005_Obs_mg/m3'].values
Sim_N_L005 = L005_biascorrect_noNan['TP_Lake_N'].values

# determine R2 for simulations related to observations at L005
slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_N_L005, Obs_L005)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_L005_N:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_L005, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_L005 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_L005 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_L005 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_L005 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_N = Sim_N_L005.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_N
correction_factor_2 = average_below_50th / average_TP_Lake_N
correction_factor_3 = average_below_75th / average_TP_Lake_N
correction_factor_4 = average_above_75th / average_TP_Lake_N

sim_BC_Qnt_EntPer = np.zeros(len(Obs_L005))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_L005)):
    if Obs_L005[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L005[i] * correction_factor_1 
    elif Obs_L005[i] > obs_quantiles[0] and Obs_L005[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L005[i] * correction_factor_2 
    elif Obs_L005[i] > obs_quantiles[1] and Obs_L005[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L005[i] * correction_factor_3
    elif Obs_L005[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L005[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_N_L005[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_L005)
r_squared = r_value**2
print("R2_L005_N_BC_Qnt:", r_squared)

# L008 north of the lake :: this station locates middle of the lake so, it is considered for both north and south bias corrections
#read observations at L008
L008_Obs = pd.read_csv('./water_quality_L008_PHOSPHATE, TOTAL AS P.csv')
L008_Obs['date'] = pd.to_datetime(L008_Obs['date'])
L008_df = pd.DataFrame()
L008_df['date'] = L008_Obs['date']
L008_df['L008_Obs_mg/m3'] = L008_Obs['L008_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
L008_df = L008_df.set_index(['date'])
L008_df.index = pd.to_datetime(L008_df.index, unit = 'ns')
L008_df_Monthly = L008_df.resample('M').mean()

L008_biascorrect = pd.merge(TP_North_Monthly,L008_df_Monthly, how='left', on='date')
L008_biascorrect_noNan = L008_biascorrect.dropna()
L008_biascorrect_noNan.to_csv('./L008_biascorrect_noNan.csv')
Obs_L008 = L008_biascorrect_noNan['L008_Obs_mg/m3'].values
Sim_N_L008 = L008_biascorrect_noNan['TP_Lake_N'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_N_L008, Obs_L008)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_L008_N:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_L008, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_L008 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_L008 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_L008 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_L008 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_N = Sim_N_L008.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_N
correction_factor_2 = average_below_50th / average_TP_Lake_N
correction_factor_3 = average_below_75th / average_TP_Lake_N
correction_factor_4 = average_above_75th / average_TP_Lake_N

sim_BC_Qnt_EntPer = np.zeros(len(Obs_L008))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_L008)):
    if Obs_L008[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L008[i] * correction_factor_1 
    elif Obs_L008[i] > obs_quantiles[0] and Obs_L008[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L008[i] * correction_factor_2 
    elif Obs_L008[i] > obs_quantiles[1] and Obs_L008[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L008[i] * correction_factor_3
    elif Obs_L008[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_N_L008[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_N_L008[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_L008)
r_squared = r_value**2
print("R2_L008_N_BC_Qnt:", r_squared)

# Station L004 South of the lake
#read observations at L004
L004_Obs = pd.read_csv('./water_quality_L004_PHOSPHATE, TOTAL AS P.csv')
L004_Obs['date'] = pd.to_datetime(L004_Obs['date'])
L004_df = pd.DataFrame()
L004_df['date'] = L004_Obs['date']
L004_df['L004_Obs_mg/m3'] = L004_Obs['L004_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
L004_df = L004_df.set_index(['date'])
L004_df.index = pd.to_datetime(L004_df.index, unit = 'ns')
L004_df_Monthly = L004_df.resample('M').mean()

L004_biascorrect = pd.merge(TP_South_Monthly,L004_df_Monthly, how='left', on='date')
L004_biascorrect_noNan = L004_biascorrect.dropna()
L004_biascorrect_noNan.to_csv('./L004_biascorrect_noNan.csv')
Obs_L004 = L004_biascorrect_noNan['L004_Obs_mg/m3'].values
Sim_S_L004 = L004_biascorrect_noNan['TP_Lake_S'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_S_L004, Obs_L004)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_L004_S:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_L004, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_L004 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_L004 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_L004 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_L004 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_S = Sim_S_L004.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_S
correction_factor_2 = average_below_50th / average_TP_Lake_S
correction_factor_3 = average_below_75th / average_TP_Lake_S
correction_factor_4 = average_above_75th / average_TP_Lake_S

sim_BC_Qnt_EntPer = np.zeros(len(Obs_L004))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_L004)):
    if Obs_L004[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L004[i] * correction_factor_1 
    elif Obs_L004[i] > obs_quantiles[0] and Obs_L004[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L004[i] * correction_factor_2 
    elif Obs_L004[i] > obs_quantiles[1] and Obs_L004[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L004[i] * correction_factor_3
    elif Obs_L004[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L004[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_S_L004[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_L004)
r_squared = r_value**2
print("R2_L004_S_BC_Qnt:", r_squared)

# Station L006 South of the lake
#read observations at L006
L006_Obs = pd.read_csv('./water_quality_L006_PHOSPHATE, TOTAL AS P.csv')
L006_Obs['date'] = pd.to_datetime(L006_Obs['date'])
L006_df = pd.DataFrame()
L006_df['date'] = L006_Obs['date']
L006_df['L006_Obs_mg/m3'] = L006_Obs['L006_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
L006_df = L006_df.set_index(['date'])
L006_df.index = pd.to_datetime(L006_df.index, unit = 'ns')
L006_df_Monthly = L006_df.resample('M').mean()

L006_biascorrect = pd.merge(TP_South_Monthly,L006_df_Monthly, how='left', on='date')
L006_biascorrect_noNan = L006_biascorrect.dropna()
L006_biascorrect_noNan.to_csv('./L006_biascorrect_noNan.csv')
Obs_L006 = L006_biascorrect_noNan['L006_Obs_mg/m3'].values
Sim_S_L006 = L006_biascorrect_noNan['TP_Lake_S'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_S_L006, Obs_L006)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_L006_S:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_L006, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_L006 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_L006 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_L006 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_L006 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_S = Sim_S_L006.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_S
correction_factor_2 = average_below_50th / average_TP_Lake_S
correction_factor_3 = average_below_75th / average_TP_Lake_S
correction_factor_4 = average_above_75th / average_TP_Lake_S

sim_BC_Qnt_EntPer = np.zeros(len(Obs_L006))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_L006)):
    if Obs_L006[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L006[i] * correction_factor_1 
    elif Obs_L006[i] > obs_quantiles[0] and Obs_L006[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L006[i] * correction_factor_2 
    elif Obs_L006[i] > obs_quantiles[1] and Obs_L006[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L006[i] * correction_factor_3
    elif Obs_L006[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L006[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_S_L006[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_L006)
r_squared = r_value**2
print("R2_L006_S_BC_Qnt:", r_squared)

# Station L007 South of the lake
#read observations at L007
L007_Obs = pd.read_csv('./water_quality_L007_PHOSPHATE, TOTAL AS P.csv')
L007_Obs['date'] = pd.to_datetime(L007_Obs['date'])
L007_df = pd.DataFrame()
L007_df['date'] = L007_Obs['date']
L007_df['L007_Obs_mg/m3'] = L007_Obs['L007_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
L007_df = L007_df.set_index(['date'])
L007_df.index = pd.to_datetime(L007_df.index, unit = 'ns')
L007_df_Monthly = L007_df.resample('M').mean()

L007_biascorrect = pd.merge(TP_South_Monthly,L007_df_Monthly, how='left', on='date')
L007_biascorrect_noNan = L007_biascorrect.dropna()
L007_biascorrect_noNan.to_csv('./L007_biascorrect_noNan.csv')
Obs_L007 = L007_biascorrect_noNan['L007_Obs_mg/m3'].values
Sim_S_L007 = L007_biascorrect_noNan['TP_Lake_S'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_S_L007, Obs_L007)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_L007_S:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_L007, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_L007 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_L007 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_L007 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_L007 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_S = Sim_S_L007.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_S
correction_factor_2 = average_below_50th / average_TP_Lake_S
correction_factor_3 = average_below_75th / average_TP_Lake_S
correction_factor_4 = average_above_75th / average_TP_Lake_S

sim_BC_Qnt_EntPer = np.zeros(len(Obs_L007))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_L007)):
    if Obs_L007[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L007[i] * correction_factor_1 
    elif Obs_L007[i] > obs_quantiles[0] and Obs_L007[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L007[i] * correction_factor_2 
    elif Obs_L007[i] > obs_quantiles[1] and Obs_L007[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L007[i] * correction_factor_3
    elif Obs_L007[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L007[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_S_L007[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_L007)
r_squared = r_value**2
print("R2_L007_S_BC_Qnt:", r_squared)

# Station L008 South
#read observations at L008
L008_Obs = pd.read_csv('./water_quality_L008_PHOSPHATE, TOTAL AS P.csv')
L008_Obs['date'] = pd.to_datetime(L008_Obs['date'])
L008_df = pd.DataFrame()
L008_df['date'] = L008_Obs['date']
L008_df['L008_Obs_mg/m3'] = L008_Obs['L008_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
L008_df = L008_df.set_index(['date'])
L008_df.index = pd.to_datetime(L008_df.index, unit = 'ns')
L008_df_Monthly = L008_df.resample('M').mean()

L008_biascorrect = pd.merge(TP_South_Monthly,L008_df_Monthly, how='left', on='date')
L008_biascorrect_noNan = L008_biascorrect.dropna()
L008_biascorrect_noNan.to_csv('./L008_biascorrect_noNan.csv')
Obs_L008 = L008_biascorrect_noNan['L008_Obs_mg/m3'].values
Sim_S_L008 = L008_biascorrect_noNan['TP_Lake_S'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_S_L008, Obs_L008)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_L008_S:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_L008, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_L008 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_L008 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_L008 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_L008 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_S = Sim_S_L008.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_S
correction_factor_2 = average_below_50th / average_TP_Lake_S
correction_factor_3 = average_below_75th / average_TP_Lake_S
correction_factor_4 = average_above_75th / average_TP_Lake_S

sim_BC_Qnt_EntPer = np.zeros(len(Obs_L008))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_L008)):
    if Obs_L008[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L008[i] * correction_factor_1 
    elif Obs_L008[i] > obs_quantiles[0] and Obs_L008[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L008[i] * correction_factor_2 
    elif Obs_L008[i] > obs_quantiles[1] and Obs_L008[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L008[i] * correction_factor_3
    elif Obs_L008[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_L008[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_S_L008[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_L008)
r_squared = r_value**2
print("R2_L008_S_BC_Qnt:", r_squared)

# Station LZ40 South of the lake
#read observations at LZ40
LZ40_Obs = pd.read_csv('./water_quality_LZ40_PHOSPHATE, TOTAL AS P.csv')
LZ40_Obs['date'] = pd.to_datetime(LZ40_Obs['date'])
LZ40_df = pd.DataFrame()
LZ40_df['date'] = LZ40_Obs['date']
LZ40_df['LZ40_Obs_mg/m3'] = LZ40_Obs['LZ40_PHOSPHATE, TOTAL AS P_mg/L']*1000 #mg/m3
LZ40_df = LZ40_df.set_index(['date'])
LZ40_df.index = pd.to_datetime(LZ40_df.index, unit = 'ns')
LZ40_df_Monthly = LZ40_df.resample('M').mean()

LZ40_biascorrect = pd.merge(TP_South_Monthly,LZ40_df_Monthly, how='left', on='date')
LZ40_biascorrect_noNan = LZ40_biascorrect.dropna()
LZ40_biascorrect_noNan.to_csv('./LZ40_biascorrect_noNan.csv')
Obs_LZ40 = LZ40_biascorrect_noNan['LZ40_Obs_mg/m3'].values
Sim_S_LZ40 = LZ40_biascorrect_noNan['TP_Lake_S'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(Sim_S_LZ40, Obs_LZ40)
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_LZ40_S:", r_squared)

# Calculate the quantiles of the observed and simulated time series.
obs_quantiles = np.percentile(Obs_LZ40, [25, 50, 75])

# Calculate the average of values less than the 25th quantile
average_below_25th = np.mean([value for value in Obs_LZ40 if value <= obs_quantiles[0]])
# Calculate the average of values less than the 50th quantile (median)
average_below_50th = np.mean([value for value in Obs_LZ40 if obs_quantiles[0] < value <= obs_quantiles[1]])
# Calculate the average of values less than the 75th quantile
average_below_75th = np.mean([value for value in Obs_LZ40 if obs_quantiles[1] < value <= obs_quantiles[2]])
# Calculate the average of values more than the 75th quantile
average_above_75th = np.mean([value for value in Obs_LZ40 if value > obs_quantiles[2]])

#Average of Sim
average_TP_Lake_S = Sim_S_LZ40.mean()

# Calculate the correction factor
correction_factor_1 = average_below_25th / average_TP_Lake_S
correction_factor_2 = average_below_50th / average_TP_Lake_S
correction_factor_3 = average_below_75th / average_TP_Lake_S
correction_factor_4 = average_above_75th / average_TP_Lake_S

sim_BC_Qnt_EntPer = np.zeros(len(Obs_LZ40))
# Apply the correction factor to TP_Lake_N
for i in range(len(Obs_LZ40)):
    if Obs_LZ40[i] <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Sim_S_LZ40[i] * correction_factor_1 
    elif Obs_LZ40[i] > obs_quantiles[0] and Obs_LZ40[i] <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Sim_S_LZ40[i] * correction_factor_2 
    elif Obs_LZ40[i] > obs_quantiles[1] and Obs_LZ40[i] <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_LZ40[i] * correction_factor_3
    elif Obs_LZ40[i] > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Sim_S_LZ40[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Sim_S_LZ40[i]


# Calculate R-squared (coefficient of determination)        
slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, Obs_LZ40)
r_squared = r_value**2
print("R2_LZ40_S_BC_Qnt:", r_squared)
