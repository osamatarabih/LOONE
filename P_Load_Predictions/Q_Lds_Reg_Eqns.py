# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 12:13:12 2023

@author: osama
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import os
os.chdir('C:/Work/Caloosahatchee_Project/P_Load_Predictions')

#S65E TP Loads
S65E_Q = pd.read_csv('./S65E_S_FLOW_cmd.csv') #in cubic meters per day
S65E_P = pd.read_csv('water_quality_S65E_PHOSPHATE, TOTAL AS P.csv') #in mg/L
S65E_merge = pd.merge(S65E_Q,S65E_P,how='left',on='date')
S65E_merge = S65E_merge.loc[:,~S65E_merge.columns.str.startswith('Unnamed')]
S65E_merge.drop(columns=['days'], inplace=True)
S65E_merge = S65E_merge.dropna()
#drop zeros 
S65E_merge = S65E_merge[(S65E_merge != 0).all(1)]
# calculate P Loads (observations)
S65E_merge['TP_Loads_Calculated'] = S65E_merge['S65E_S_FLOW_cfs'] * S65E_merge['S65E_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
# LOG Q and LOG P
S65E_merge['S65E_Q_Log'] = np.log(S65E_merge['S65E_S_FLOW_cfs'])
S65E_merge['TP_Loads_Calculated_Log'] = np.log(S65E_merge['TP_Loads_Calculated'])
# Calculate P Loads as an exponential function of Flow (All in LOG)
S65E_merge['TP_Loads_Predicted_Log'] = 2.00040151533473* (S65E_merge['S65E_Q_Log']**0.837387838314323)
#Evaluate the Model
r_squared = r2_score(S65E_merge['TP_Loads_Calculated_Log'], S65E_merge['TP_Loads_Predicted_Log'])
print("R-squared_S65E(logged):", r_squared)
# Reverse the LOG
S65E_merge['TP_Loads_Predicted'] = np.exp(S65E_merge['TP_Loads_Predicted_Log'])

#S71 TP Loads
S71_Q = pd.read_csv('./S71_S_FLOW_cmd.csv')
S71_P = pd.read_csv('water_quality_S71_PHOSPHATE, TOTAL AS P.csv')
S71_merge = pd.merge(S71_Q,S71_P,how='left',on='date')
S71_merge = S71_merge.loc[:,~S71_merge.columns.str.startswith('Unnamed')]
S71_merge.drop(columns=['days'], inplace=True)
S71_merge = S71_merge.dropna()
S71_merge = S71_merge[(S71_merge != 0).all(1)]
S71_merge['TP_Loads_Calculated'] = S71_merge['S71_S_FLOW_cfs'] * S71_merge['S71_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S71_merge['S71_Q_Log'] = np.log(S71_merge['S71_S_FLOW_cfs'])
S71_merge['TP_Loads_Calculated_Log'] = np.log(S71_merge['TP_Loads_Calculated'])
S71_merge['TP_Loads_Predicted_Log'] = 2.55809777403484* (S71_merge['S71_Q_Log']**0.765894033054918)
#Evaluate the Model
r_squared = r2_score(S71_merge['TP_Loads_Calculated_Log'], S71_merge['TP_Loads_Predicted_Log'])
print("R-squared_S71(logged):", r_squared)
# Reverse the LOG
S71_merge['TP_Loads_Predicted'] = np.exp(S71_merge['TP_Loads_Predicted_Log'])

#S72 TP Loads
S72_Q = pd.read_csv('./S72_S_FLOW_cmd.csv')
S72_P = pd.read_csv('water_quality_S72_PHOSPHATE, TOTAL AS P.csv')
S72_merge = pd.merge(S72_Q,S72_P,how='left',on='date')
S72_merge = S72_merge.loc[:,~S72_merge.columns.str.startswith('Unnamed')]
S72_merge.drop(columns=['days'], inplace=True)
S72_merge = S72_merge.dropna()
S72_merge = S72_merge[(S72_merge != 0).all(1)]
S72_merge['TP_Loads_Calculated'] = S72_merge['S72_S_FLOW_cfs'] * S72_merge['S72_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S72_merge['S72_Q_Log'] = np.log(S72_merge['S72_S_FLOW_cfs'])
S72_merge['TP_Loads_Calculated_Log'] = np.log(S72_merge['TP_Loads_Calculated'])
S72_merge['TP_Loads_Predicted_Log'] = 2.85270576092534* (S72_merge['S72_Q_Log']**0.724935760736887)
#Evaluate the Model
r_squared = r2_score(S72_merge['TP_Loads_Calculated_Log'], S72_merge['TP_Loads_Predicted_Log'])
print("R-squared_S72(logged):", r_squared)
# Reverse the LOG
S72_merge['TP_Loads_Predicted'] = np.exp(S72_merge['TP_Loads_Predicted_Log'])

#S191 TP Loads
S191_Q = pd.read_csv('./S191_S_FLOW_cfs.csv')
S191_P = pd.read_csv('./water_quality_S191_PHOSPHATE, TOTAL AS P.csv')
S191_merge = pd.merge(S191_Q,S191_P,how='left',on='date')
S191_merge = S191_merge.loc[:,~S191_merge.columns.str.startswith('Unnamed')]
S191_merge.drop(columns=['days'], inplace=True)
S191_merge = S191_merge.dropna()
S191_merge = S191_merge[(S191_merge != 0).all(1)]
S191_merge['TP_Loads_Calculated'] = S191_merge['S191_S_FLOW_cfs'] * S191_merge['S191_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S191_merge['S191_Q_Log'] = np.log(S191_merge['S191_S_FLOW_cfs'])
S191_merge['TP_Loads_Calculated_Log'] = np.log(S191_merge['TP_Loads_Calculated'])
S191_merge['TP_Loads_Predicted_Log'] = 3.0257439276073* (S191_merge['S191_Q_Log']**0.721906661127014)
#Evaluate the Model
r_squared = r2_score(S191_merge['TP_Loads_Calculated_Log'], S191_merge['TP_Loads_Predicted_Log'])
print("R-squared_S191(logged):", r_squared)
# Reverse the LOG
S191_merge['TP_Loads_Predicted'] = np.exp(S191_merge['TP_Loads_Predicted_Log'])

#FECSR78 TP Loads
FECSR78_Q = pd.read_csv('./FISHP_FLOW_cmd.csv')
FECSR78_P = pd.read_csv('./water_quality_FECSR78_PHOSPHATE, TOTAL AS P.csv')
FECSR78_merge = pd.merge(FECSR78_Q,FECSR78_P,how='left',on='date')
FECSR78_merge = FECSR78_merge.loc[:,~FECSR78_merge.columns.str.startswith('Unnamed')]
FECSR78_merge.drop(columns=['days'], inplace=True)
FECSR78_merge = FECSR78_merge.dropna()
FECSR78_merge = FECSR78_merge[(FECSR78_merge != 0).all(1)]
FECSR78_merge['TP_Loads_Calculated'] = FECSR78_merge['FISHP_FLOW_cfs'] * FECSR78_merge['FECSR78_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
FECSR78_merge['FECSR78_Q_Log'] = np.log(FECSR78_merge['FISHP_FLOW_cfs'])
FECSR78_merge['TP_Loads_Calculated_Log'] = np.log(FECSR78_merge['TP_Loads_Calculated'])
FECSR78_merge['TP_Loads_Predicted_Log'] = 2.59223308404186* (FECSR78_merge['FECSR78_Q_Log']**0.756802713030507)
#Evaluate the Model
r_squared = r2_score(FECSR78_merge['TP_Loads_Calculated_Log'], FECSR78_merge['TP_Loads_Predicted_Log'])
print("R-squared_FECSR(logged):", r_squared)
# Reverse the LOG
FECSR78_merge['TP_Loads_Predicted'] = np.exp(FECSR78_merge['TP_Loads_Predicted_Log'])

#S4_P TP Loads
S4_P_Q = pd.read_csv('./S4_P_FLOW_cmd.csv')
S4_P_P = pd.read_csv('./water_quality_S4_PHOSPHATE, TOTAL AS P.csv')
S4_P_merge = pd.merge(S4_P_Q,S4_P_P,how='left',on='date')
S4_P_merge = S4_P_merge.loc[:,~S4_P_merge.columns.str.startswith('Unnamed')]
S4_P_merge.drop(columns=['days'], inplace=True)
S4_P_merge = S4_P_merge.dropna()
S4_P_merge = S4_P_merge[(S4_P_merge != 0).all(1)]
S4_P_merge['TP_Loads_Calculated'] = S4_P_merge['S4_P_FLOW_cfs'] * S4_P_merge['S4_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S4_P_merge['S4_P_Q_Log'] = np.log(S4_P_merge['S4_P_FLOW_cfs'])
S4_P_merge['TP_Loads_Calculated_Log'] = np.log(S4_P_merge['TP_Loads_Calculated'])
S4_P_merge['TP_Loads_Predicted_Log'] = 2.86495657296006* (S4_P_merge['S4_P_Q_Log']**0.72203267810211)
#Evaluate the Model
r_squared = r2_score(S4_P_merge['TP_Loads_Calculated_Log'], S4_P_merge['TP_Loads_Predicted_Log'])
print("R-squared_S4(logged):", r_squared)
# Reverse the LOG
S4_P_merge['TP_Loads_Predicted'] = np.exp(S4_P_merge['TP_Loads_Predicted_Log'])

#S84 TP Loads
S84_Q = pd.read_csv('./S84_S_FLOW_cmd.csv')
S84_P = pd.read_csv('./water_quality_S84_PHOSPHATE, TOTAL AS P.csv')
S84_merge = pd.merge(S84_Q,S84_P,how='left',on='date')
S84_merge = S84_merge.loc[:,~S84_merge.columns.str.startswith('Unnamed')]
S84_merge.drop(columns=['days'], inplace=True)
S84_merge = S84_merge.dropna()
S84_merge = S84_merge[(S84_merge != 0).all(1)]
S84_merge['TP_Loads_Calculated'] = S84_merge['S84_S_FLOW_cfs'] * S84_merge['S84_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S84_merge['S84_Q_Log'] = np.log(S84_merge['S84_S_FLOW_cfs'])
S84_merge['TP_Loads_Calculated_Log'] = np.log(S84_merge['TP_Loads_Calculated'])
S84_merge['TP_Loads_Predicted_Log'] = 2.53265243618408* (S84_merge['S84_Q_Log']**0.750938593484588)
#Evaluate the Model
r_squared = r2_score(S84_merge['TP_Loads_Calculated_Log'], S84_merge['TP_Loads_Predicted_Log'])
print("R-squared_S84(logged):", r_squared)
# Reverse the LOG
S84_merge['TP_Loads_Predicted'] = np.exp(S84_merge['TP_Loads_Predicted_Log'])

#S127_P TP Loads
S127_P_Q = pd.read_csv('./S127_P_FLOW_cmd.csv')
S127_P_P = pd.read_csv('./water_quality_S127_PHOSPHATE, TOTAL AS P.csv')
S127_P_merge = pd.merge(S127_P_Q,S127_P_P,how='left',on='date')
S127_P_merge = S127_P_merge.loc[:,~S127_P_merge.columns.str.startswith('Unnamed')]
S127_P_merge.drop(columns=['days'], inplace=True)
S127_P_merge = S127_P_merge.dropna()
S127_P_merge = S127_P_merge[(S127_P_merge != 0).all(1)]
S127_P_merge['TP_Loads_Calculated'] = S127_P_merge['S127_P_FLOW_cfs'] * S127_P_merge['S127_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S127_P_merge['S127_P_Q_Log'] = np.log(S127_P_merge['S127_P_FLOW_cfs'])
S127_P_merge['TP_Loads_Calculated_Log'] = np.log(S127_P_merge['TP_Loads_Calculated'])
S127_P_merge['TP_Loads_Predicted_Log'] = 2.34697955615531* (S127_P_merge['S127_P_Q_Log']**0.794046635942522)
#Evaluate the Model
r_squared = r2_score(S127_P_merge['TP_Loads_Calculated_Log'], S127_P_merge['TP_Loads_Predicted_Log'])
print("R-squared_S127P(logged):", r_squared)
# Reverse the LOG
S127_P_merge['TP_Loads_Predicted'] = np.exp(S127_P_merge['TP_Loads_Predicted_Log'])

#S127_C TP Loads
S127_C_Q = pd.read_csv('./S127_C_FLOW_cmd.csv')
S127_C_P = pd.read_csv('./water_quality_S127_PHOSPHATE, TOTAL AS P.csv')
S127_C_merge = pd.merge(S127_C_Q,S127_C_P,how='left',on='date')
S127_C_merge = S127_C_merge.loc[:,~S127_C_merge.columns.str.startswith('Unnamed')]
S127_C_merge.drop(columns=['days'], inplace=True)
S127_C_merge = S127_C_merge[S127_C_merge>0]
S127_C_merge = S127_C_merge.dropna()
S127_C_merge = S127_C_merge[(S127_C_merge != 0).all(1)]
S127_C_merge['TP_Loads_Calculated'] = S127_C_merge['S127_C_FLOW_cfs'] * S127_C_merge['S127_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S127_C_merge['S127_C_Q_Log'] = np.log(S127_C_merge['S127_C_FLOW_cfs'])
S127_C_merge['TP_Loads_Calculated_Log'] = np.log(S127_C_merge['TP_Loads_Calculated'])
S127_C_merge['TP_Loads_Predicted_Log'] = 2.73825064156312* (S127_C_merge['S127_C_Q_Log']**0.715023290260209)
#Evaluate the Model
r_squared = r2_score(S127_C_merge['TP_Loads_Calculated_Log'], S127_C_merge['TP_Loads_Predicted_Log'])
print("R-squared_S127C(logged):", r_squared)
# Reverse the LOG
S127_C_merge['TP_Loads_Predicted'] = np.exp(S127_C_merge['TP_Loads_Predicted_Log'])

#S133_P TP Loads
S133_P_Q = pd.read_csv('./S133_P_FLOW_cmd.csv')
S133_P_P = pd.read_csv('./water_quality_S133_PHOSPHATE, TOTAL AS P.csv')
S133_P_merge = pd.merge(S133_P_Q,S133_P_P,how='left',on='date')
S133_P_merge = S133_P_merge.loc[:,~S133_P_merge.columns.str.startswith('Unnamed')]
S133_P_merge.drop(columns=['days'], inplace=True)
S133_P_merge = S133_P_merge.dropna()
S133_P_merge = S133_P_merge[(S133_P_merge != 0).all(1)]
S133_P_merge['TP_Loads_Calculated'] = S133_P_merge['S133_P_FLOW_cfs'] * S133_P_merge['S133_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S133_P_merge['S133_P_Q_Log'] = np.log(S133_P_merge['S133_P_FLOW_cfs'])
S133_P_merge['TP_Loads_Calculated_Log'] = np.log(S133_P_merge['TP_Loads_Calculated'])
S133_P_merge['TP_Loads_Predicted_Log'] = 2.64107054734111* (S133_P_merge['S133_P_Q_Log']**0.756152588482486)
#Evaluate the Model
r_squared = r2_score(S133_P_merge['TP_Loads_Calculated_Log'], S133_P_merge['TP_Loads_Predicted_Log'])
print("R-squared_S133P(logged):", r_squared)
# Reverse the LOG
S133_P_merge['TP_Loads_Predicted'] = np.exp(S133_P_merge['TP_Loads_Predicted_Log'])

#S154_C TP Loads
S154_C_Q = pd.read_csv('./S154_C_FLOW_cmd.csv')
S154_C_P = pd.read_csv('./water_quality_S154_PHOSPHATE, TOTAL AS P.csv')
S154_C_merge = pd.merge(S154_C_Q,S154_C_P,how='left',on='date')
S154_C_merge = S154_C_merge.loc[:,~S154_C_merge.columns.str.startswith('Unnamed')]
S154_C_merge.drop(columns=['days'], inplace=True)
S154_C_merge = S154_C_merge[S154_C_merge>0]
S154_C_merge = S154_C_merge.dropna()
S154_C_merge = S154_C_merge[(S154_C_merge != 0).all(1)]
S154_C_merge['TP_Loads_Calculated'] = S154_C_merge['S154_C_FLOW_cfs'] * S154_C_merge['S154_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S154_C_merge['S154_C_Q_Log'] = np.log(S154_C_merge['S154_C_FLOW_cfs'])
S154_C_merge['TP_Loads_Calculated_Log'] = np.log(S154_C_merge['TP_Loads_Calculated'])
S154_C_merge['TP_Loads_Predicted_Log'] = 3.10305150879462* (S154_C_merge['S154_C_Q_Log']**0.7099895764193)
#Evaluate the Model
r_squared = r2_score(S154_C_merge['TP_Loads_Calculated_Log'], S154_C_merge['TP_Loads_Predicted_Log'])
print("R-squared_S154(logged):", r_squared)
# Reverse the LOG
S154_C_merge['TP_Loads_Predicted'] = np.exp(S154_C_merge['TP_Loads_Predicted_Log'])

#S135_P TP Loads
S135_P_Q = pd.read_csv('./S135_PMP_P_FLOW_cmd.csv')
S135_P_P = pd.read_csv('./water_quality_S135_PHOSPHATE, TOTAL AS P.csv')
S135_P_merge = pd.merge(S135_P_Q,S135_P_P,how='left',on='date')
S135_P_merge = S135_P_merge.loc[:,~S135_P_merge.columns.str.startswith('Unnamed')]
S135_P_merge.drop(columns=['days'], inplace=True)
S135_P_merge = S135_P_merge.dropna()
S135_P_merge = S135_P_merge[(S135_P_merge != 0).all(1)]
S135_P_merge['TP_Loads_Calculated'] = S135_P_merge['S135 PMP_P_FLOW_cfs'] * S135_P_merge['S135_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S135_P_merge['S135_P_Q_Log'] = np.log(S135_P_merge['S135 PMP_P_FLOW_cfs'])
S135_P_merge['TP_Loads_Calculated_Log'] = np.log(S135_P_merge['TP_Loads_Calculated'])
S135_P_merge['TP_Loads_Predicted_Log'] = 2.50975664040355* (S135_P_merge['S135_P_Q_Log']**0.760702496334553)
#Evaluate the Model
r_squared = r2_score(S135_P_merge['TP_Loads_Calculated_Log'], S135_P_merge['TP_Loads_Predicted_Log'])
print("R-squared_S135P(logged):", r_squared)
# Reverse the LOG
S135_P_merge['TP_Loads_Predicted'] = np.exp(S135_P_merge['TP_Loads_Predicted_Log'])

#S135_C TP Loads
S135_C_Q = pd.read_csv('./S135_C_FLOW_cmd.csv')
S135_C_P = pd.read_csv('./water_quality_S135_PHOSPHATE, TOTAL AS P.csv')
S135_C_merge = pd.merge(S135_C_Q,S135_C_P,how='left',on='date')
S135_C_merge = S135_C_merge.loc[:,~S135_C_merge.columns.str.startswith('Unnamed')]
S135_C_merge.drop(columns=['days'], inplace=True)
S135_C_merge = S135_C_merge[S135_C_merge>0]
S135_C_merge = S135_C_merge.dropna()
S135_C_merge = S135_C_merge[(S135_C_merge != 0).all(1)]
S135_C_merge['TP_Loads_Calculated'] = S135_C_merge['S135_C_FLOW_cfs'] * S135_C_merge['S135_PHOSPHATE, TOTAL AS P_mg/L'] * 1000
S135_C_merge['S135_C_Q_Log'] = np.log(S135_C_merge['S135_C_FLOW_cfs'])
S135_C_merge['TP_Loads_Calculated_Log'] = np.log(S135_C_merge['TP_Loads_Calculated'])
S135_C_merge['TP_Loads_Predicted_Log'] = 2.43076251736749* (S135_C_merge['S135_C_Q_Log']**0.759494593788417)
#Evaluate the Model
r_squared = r2_score(S135_C_merge['TP_Loads_Calculated_Log'], S135_C_merge['TP_Loads_Predicted_Log'])
print("R-squared_S135C(logged):", r_squared)
# Reverse the LOG
S135_C_merge['TP_Loads_Predicted'] = np.exp(S135_C_merge['TP_Loads_Predicted_Log'])