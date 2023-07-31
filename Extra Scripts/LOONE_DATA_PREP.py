# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 00:18:50 2023

@author: osama
"""

import pandas as pd
import numpy as np
from datetime import date, timedelta
import os
#Directory where the Script loactes
os.chdir('C:/Work/Research/Data Analysis/Tools/Python_Scripts')
from Data_Analyses_Fns import *
Working_dir = 'C:/Work/Research/LOONE' 
os.chdir('%s'%Working_dir) 
from Stg_Sto_Ar import Stg_Sto_Ar

M3_Yr = 2007
M3_M = 10
M3_D = 1
D2_Yr = 2007
D2_M = 12
D2_D = 30
St_Yr = 2008
St_M = 1
St_D = 1
En_Yr = 2022
En_M = 12
En_D = 31
# To create File (Average_LO_Storage_LORS20082023)
#Read LO Average Stage (ft)
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
LO_Stage = pd.read_csv('./LO_Stage_2023.csv')
# Create Column (EOD Stg(ft,NGVD)) in File (SFWMM_Daily_Outputs_LORS20082023)
LO_Stage = DF_Date_Range(LO_Stage, M3_Yr, M3_M, M3_D, En_Yr, En_M, En_D)
LO_Storage = Stg_Sto_Ar.stg2sto(LO_Stage['Average_Stage'], 0)
LO_SA = Stg_Sto_Ar.stg2ar(LO_Stage['Average_Stage'], 0)
LO_Stg_Sto_SA_df = pd.DataFrame(LO_Stage['date'],columns=['date'])
LO_Stg_Sto_SA_df['Stage_ft'] = LO_Stage['Average_Stage']
LO_Stg_Sto_SA_df['Stage_m'] = LO_Stg_Sto_SA_df['Stage_ft'].values * 0.3048 #ft to m
LO_Stg_Sto_SA_df['Storage_acft'] = LO_Storage
LO_Stg_Sto_SA_df['Storage_cmd'] = LO_Stg_Sto_SA_df['Storage_acft'] * 1233.48 #acft to m3/d
LO_Stg_Sto_SA_df['SA_acres'] = LO_SA #acres
 
#Read flow data cubic meters per day
Working_dir = 'C:/Work/Research/LOONE/Inflow_Data_2023' 
os.chdir('%s'%Working_dir) 
S65_total = pd.read_csv('./S65E_total.csv')
S65_total["S65E_tot_cmd"] = S65_total[["S65E_S_FLOW_cfs", "S65EX1_S_FLOW_cfs"]].sum(axis=1)
S71_S = pd.read_csv('./S71_S_FLOW_cmd.csv')
S72_S = pd.read_csv('./S72_S_FLOW_cmd.csv')
S84_S = pd.read_csv('./S84_S_FLOW_cmd.csv')
S127_C = pd.read_csv('./S127_C_FLOW_cmd.csv')
S127_P = pd.read_csv('./S127_P_FLOW_cmd.csv')
S129_C = pd.read_csv('./S129_C_FLOW_cmd.csv')
S129_P = pd.read_csv('./S129_PMP_P_FLOW_cmd.csv')
S133_P = pd.read_csv('./S133_P_FLOW_cmd.csv')
S135_C = pd.read_csv('./S135_C_FLOW_cmd.csv')
S135_P = pd.read_csv('./S135_PMP_P_FLOW_cmd.csv')
S154_C = pd.read_csv('./S154_C_FLOW_cmd.csv')
S191_S = pd.read_csv('./S191_S_FLOW_cmd.csv')
S308 = pd.read_csv('./S308.DS_FLOW_cmd.csv')
S351_S = pd.read_csv('./S351_S_FLOW_cmd.csv')
S352_S = pd.read_csv('./S352_S_FLOW_cmd.csv')
S354_S = pd.read_csv('./S354_S_FLOW_cmd.csv')
FISHP = pd.read_csv('./FISHP_FLOW_cmd.csv')
L8 = pd.read_csv('./L8.441_FLOW_cmd.csv')
S2_P = pd.read_csv('./S2_P_FLOW_cmd.csv')
S3_P = pd.read_csv('./S3_P_FLOW_cmd.csv')
S4_P = pd.read_csv('./S4_P_FLOW_cmd.csv')

Working_dir = 'C:/Work/Research/LOONE/Outflow_Data_2023' 
os.chdir('%s'%Working_dir) 
S77_S = pd.read_csv('./S77_S_FLOW_cmd.csv')
INDUST = pd.read_csv('./INDUST_FLOW_cmd.csv')

#Read Interpolated TP data
# Data_Interpolation Python Script is used to interpolate TP data for all inflow stations addressed below!
Working_dir = 'C:/Work/Research/LOONE/WQ_Jun23' 
os.chdir('%s'%Working_dir) 
S65_total_TP = pd.read_csv('./water_quality_S65E_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S71_TP = pd.read_csv('./water_quality_S71_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S72_TP = pd.read_csv('./water_quality_S72_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S84_TP = pd.read_csv('./water_quality_S84_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S127_TP = pd.read_csv('./water_quality_S127_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S133_TP = pd.read_csv('./water_quality_S133_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S135_TP = pd.read_csv('./water_quality_S135_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S154_TP = pd.read_csv('./water_quality_S154_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S191_TP = pd.read_csv('./water_quality_S191_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S308_TP = pd.read_csv('./water_quality_S308C_PHOSPHATE, TOTAL AS P_Interpolated.csv')
FISHP_TP = pd.read_csv('./water_quality_FECSR78_PHOSPHATE, TOTAL AS P_Interpolated.csv')
L8_TP = pd.read_csv('./water_quality_CULV10A_PHOSPHATE, TOTAL AS P_Interpolated.csv')
S4_TP = pd.read_csv('./water_quality_S4_PHOSPHATE, TOTAL AS P_Interpolated.csv')

#Set date range for S65 TP
S65_total_TP = DF_Date_Range(S65_total_TP, M3_Yr, M3_M, M3_D, En_Yr, En_M, En_D)

#Set Date Range 
Q_names = ['S65_Q','S71_Q', 'S72_Q','S84_Q','S127_C_Q','S127_P_Q','S129_C_Q','S129_P_Q','S133_P_Q','S135_C_Q','S135_P_Q','S154_Q','S191_Q',
            'S308_Q','S351_Q','S352_Q','S354_Q','FISHP_Q','L8_Q','S2_P_Q','S3_P_Q','S4_P_Q','S77_Q','INDUST_Q']
Q_list = {'S65_Q':S65_total,'S71_Q':S71_S,'S72_Q':S72_S,'S84_Q':S84_S,'S127_C_Q':S127_C,'S127_P_Q':S127_P,'S129_C_Q':S129_C,
           'S129_P_Q':S129_P,'S133_P_Q':S133_P,'S135_C_Q':S135_C,'S135_P_Q':S135_P,'S154_Q':S154_C,'S191_Q':S191_S,'S308_Q':S308,
           'S351_Q':S351_S,'S352_Q':S352_S,'S354_Q':S354_S,'FISHP_Q':FISHP,'L8_Q':L8,'S2_P_Q':S2_P,'S3_P_Q':S3_P,'S4_P_Q':S4_P,
           'S77_Q':S77_S,'INDUST_Q':INDUST}
# Identify date range
date = pd.date_range(start = '%s/%s/%s'%(M3_M,M3_D,M3_Yr),end = '%s/%s/%s'%(En_M,En_D,En_Yr),freq = 'D')
# Create Flow Dataframe
Flow_df = pd.DataFrame(date, columns=['date'])
for i in range(len(Q_names)):
    x = DF_Date_Range(Q_list[Q_names[i]], M3_Yr, M3_M, M3_D, En_Yr, En_M, En_D)
    Flow_df['%s'%Q_names[i]] = x.iloc[:,-1:].values
Flow_df['S127_C_Q'] = Flow_df['S127_C_Q'][Flow_df['S127_C_Q']>=0]
Flow_df['S127_C_Q'] = Flow_df['S127_C_Q'].fillna(0)
Flow_df['S127_In'] = Flow_df[["S127_C_Q", "S127_P_Q"]].sum(axis=1)
Flow_df['S129_C_Q'] = Flow_df['S129_C_Q'][Flow_df['S129_C_Q']>=0]
Flow_df['S129_C_Q'] = Flow_df['S129_C_Q'].fillna(0)
Flow_df['S129_In'] = Flow_df[["S129_C_Q", "S129_P_Q"]].sum(axis=1)
Flow_df['S135_C_Q'] = Flow_df['S135_C_Q'][Flow_df['S135_C_Q']>=0]
Flow_df['S135_C_Q'] = Flow_df['S135_C_Q'].fillna(0)
Flow_df['S135_In'] = Flow_df[["S135_C_Q", "S135_P_Q"]].sum(axis=1)
Flow_df['S308_In'] = Flow_df['S308_Q'][Flow_df['S308_Q']<0]
Flow_df['S308_In'] = Flow_df['S308_In'] * -1
Flow_df['S308_In'] = Flow_df['S308_In'].fillna(0)
Flow_df['S77_In'] = Flow_df['S77_Q'][Flow_df['S77_Q']<0]
Flow_df['S77_In'] = Flow_df['S77_In'] * -1
Flow_df['S77_In'] = Flow_df['S77_In'].fillna(0)
Flow_df['S351_In'] = Flow_df['S351_Q'][Flow_df['S351_Q']<0]
Flow_df['S351_In'] = Flow_df['S351_In'] * -1
Flow_df['S351_In'] = Flow_df['S351_In'].fillna(0)
Flow_df['S352_In'] = Flow_df['S352_Q'][Flow_df['S352_Q']<0]
Flow_df['S352_In'] = Flow_df['S352_In'] * -1
Flow_df['S352_In'] = Flow_df['S352_In'].fillna(0)
Flow_df['S354_In'] = Flow_df['S354_Q'][Flow_df['S354_Q']<0]
Flow_df['S354_In'] = Flow_df['S354_In'] * -1
Flow_df['S354_In'] = Flow_df['S354_In'].fillna(0)
Flow_df['L8_In'] = Flow_df['L8_Q'][Flow_df['L8_Q']<0]
Flow_df['L8_In'] = Flow_df['L8_In'] * -1
Flow_df['L8_In'] = Flow_df['L8_In'].fillna(0)
Flow_df['S308_Out'] = Flow_df['S308_Q'][Flow_df['S308_Q']>=0]
Flow_df['S308_Out'] = Flow_df['S308_Out'].fillna(0)
Flow_df['S77_Out'] = Flow_df['S77_Q'][Flow_df['S77_Q']>=0]
Flow_df['S77_Out'] = Flow_df['S77_Out'].fillna(0)
Flow_df['INDUST_Out'] = Flow_df['INDUST_Q'][Flow_df['INDUST_Q']>=0]
Flow_df['INDUST_Out'] = Flow_df['INDUST_Out'].fillna(0)
Flow_df['S351_Out'] = Flow_df['S351_Q'][Flow_df['S351_Q']>=0]
Flow_df['S351_Out'] = Flow_df['S351_Out'].fillna(0)
Flow_df['S352_Out'] = Flow_df['S352_Q'][Flow_df['S352_Q']>=0]
Flow_df['S352_Out'] = Flow_df['S352_Out'].fillna(0)
Flow_df['S354_Out'] = Flow_df['S354_Q'][Flow_df['S354_Q']>=0]
Flow_df['S354_Out'] = Flow_df['S354_Out'].fillna(0)
Flow_df['L8_Out'] = Flow_df['L8_Q'][Flow_df['L8_Q']>=0]
Flow_df['L8_Out'] = Flow_df['L8_Out'].fillna(0)
Flow_df['Inflows'] = Flow_df[["S65_Q", "S71_Q",'S72_Q','S84_Q','S127_In','S129_In','S133_P_Q','S135_In',
                              'S154_Q','S191_Q','S308_In','S77_In','S351_In','S352_In','S354_In','L8_In','FISHP_Q',
                              'S2_P_Q','S3_P_Q','S4_P_Q']].sum(axis=1)
Flow_df['Netflows'] = Flow_df['Inflows'] - Flow_df['INDUST_Out']
Flow_df['Outflows'] = Flow_df[["S308_Out", "S77_Out",'S351_Out','S352_Out','S354_Out','INDUST_Out','L8_Out']].sum(axis=1)
TP_names = ['S65_TP','S71_TP','S72_TP','S84_TP','S127_TP','S133_TP','S135_TP','S154_TP','S191_TP','S308_TP','FISHP_TP','L8_TP','S4_TP']
TP_list = {'S65_TP':S65_total_TP,'S71_TP':S71_TP,'S72_TP':S72_TP,'S84_TP':S84_TP,'S127_TP':S127_TP,'S133_TP':S133_TP,'S135_TP':S135_TP,
           'S154_TP':S154_TP,'S191_TP':S191_TP,'S308_TP':S308_TP,'FISHP_TP':FISHP_TP,'L8_TP':L8_TP,'S4_TP':S4_TP}
#Create TP Concentrations Dataframe
TP_df = pd.DataFrame(date, columns=['date'])
for i in range(len(TP_names)):
    y = DF_Date_Range(TP_list[TP_names[i]], M3_Yr, M3_M, M3_D, En_Yr, En_M, En_D)
    TP_df['%s'%TP_names[i]] = y.iloc[:,-1:].values

#Determine TP Loads (mg)
TP_Loads_In = pd.DataFrame(date, columns=['date'])
TP_Loads_In['S65_P_Ld'] = Flow_df['S65_Q'] * TP_df['S65_TP'] * 1000 #(m3/d * mg/L * 1000 = mg/d)
TP_Loads_In['S71_P_Ld'] = Flow_df['S71_Q'] * TP_df['S71_TP'] * 1000
TP_Loads_In['S72_P_Ld'] = Flow_df['S72_Q'] * TP_df['S72_TP'] * 1000
TP_Loads_In['S84_P_Ld'] = Flow_df['S84_Q'] * TP_df['S84_TP'] * 1000
TP_Loads_In['S127_P_Ld'] = Flow_df['S127_In'] * TP_df['S127_TP'] * 1000
TP_Loads_In['S133_P_Ld'] = Flow_df['S133_P_Q'] * TP_df['S133_TP'] * 1000
TP_Loads_In['S135_P_Ld'] = Flow_df['S135_In'] * TP_df['S135_TP'] * 1000
TP_Loads_In['S154_P_Ld'] = Flow_df['S154_Q'] * TP_df['S154_TP'] * 1000
TP_Loads_In['S191_P_Ld'] = Flow_df['S191_Q'] * TP_df['S191_TP'] * 1000
TP_Loads_In['S308_P_Ld'] = Flow_df['S308_In'] * TP_df['S308_TP'] * 1000
TP_Loads_In['FISHP_P_Ld'] = Flow_df['FISHP_Q'] * TP_df['FISHP_TP'] * 1000
TP_Loads_In['L8_P_Ld'] = Flow_df['L8_In'] * TP_df['L8_TP'] * 1000
TP_Loads_In['S4_P_Ld'] = Flow_df['S4_P_Q'] * TP_df['S4_TP'] * 1000
#Calculate the total External Loads to Lake Okeechobee
TP_Loads_In['External_P_Ld_mg'] = TP_Loads_In.sum(axis=1)

# Create File (LO_External_Loadings_3MLag_LORS20082023)
TP_Loads_In_3MLag = DF_Date_Range(TP_Loads_In, M3_Yr, M3_M, M3_D, En_Yr, En_M, En_D)
TP_Loads_In_3MLag_df = pd.DataFrame(TP_Loads_In_3MLag['date'],columns=['date'])
TP_Loads_In_3MLag_df['External_Loads'] = TP_Loads_In_3MLag['External_P_Ld_mg']

#Create File (LO_Inflows_BK_LORS20082023)
LO_Inflows_BK = pd.DataFrame(Flow_df['date'],columns=['date'])
LO_Inflows_BK['Inflows_cmd'] = Flow_df['Inflows']
LO_Inflows_BK = DF_Date_Range(LO_Inflows_BK, St_Yr, St_M, St_D, En_Yr, En_M, En_D)

#Create File (Outflows_consd_20082023)
Outflows_consd= pd.DataFrame(Flow_df['date'],columns=['date'])
Outflows_consd['Outflows_acft'] = Flow_df['Outflows']/1233.48 #acft 
Outflows_consd = DF_Date_Range(Outflows_consd, St_Yr, St_M, St_D, En_Yr, En_M, En_D)

# Create File (INDUST_Outflow_20082023) 
INDUST_Outflows = pd.DataFrame(Flow_df['date'],columns=['date'])
INDUST_Outflows['INDUST'] = Flow_df['INDUST_Out']

#Create File (Netflows_acft_LORS20082023)
# This is also Column (Net Inflow) in File (SFWMM_Daily_Outputs_LORS20082023)
Netflows = pd.DataFrame(Flow_df['date'],columns=['date'])
Netflows['Netflows_acft'] = Flow_df['Netflows']/1233.48 #acft 
Netflows = DF_Date_Range(Netflows, D2_Yr, D2_M, D2_D, En_Yr, En_M, En_D)


# Create File (TotalQWCA_Obs_LORS20082023)
# This is also Column (RegWCA) in File (SFWMM_Daily_Outputs_LORS20082023)
TotalQWCA = pd.DataFrame(Flow_df['date'],columns=['date'])
TotalQWCA['S351_Out'] = Flow_df['S351_Out'] * (35.3147/86400) #cmd to cfs
TotalQWCA['S354_Out'] = Flow_df['S354_Out'] * (35.3147/86400) 
TotalQWCA['RegWCA_cfs'] = TotalQWCA.sum(axis=1) #cfs
TotalQWCA['RegWCA_acft'] = TotalQWCA['RegWCA_cfs'] *1.9835 #acft
TotalQWCA = DF_Date_Range(TotalQWCA, D2_Yr, D2_M, D2_D, En_Yr, En_M, En_D)

# Create Column (RegL8C51) in the File (SFWMM_Daily_Outputs_LORS20082023)
L8C51 = pd.DataFrame(Flow_df['date'],columns=['date'])
L8C51['S352_Out'] = Flow_df['S352_Out'].values * (35.3147/86400) #cmd to cfs
L8C51['L8_O_cfs'] = Flow_df['L8_Out'].values * (35.3147/86400) #cmd to cfs
L8C51['L8C51_cfs'] = L8C51.sum(axis=1) #cfs
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
L8C51.to_csv('./L8C51.csv')

#C43 RO C44 RO
#Create Files (C43RO_LORS20082023, C43RO_Monthly_LORS20082023,C44RO_LORS20082023, C44RO_Monthly_LORS20082023)
#As well as Columns C43Runoff and C44Runoff in File (SFWMM_Daily_Outputs_LORS20082023)
Working_dir = 'C:/Work/Research/LOONE/Outflow_Data_2023' 
os.chdir('%s'%Working_dir) 
S79 = pd.read_csv('./S79.csv')
S79 = S79.fillna(0)
S80 = pd.read_csv('./S80.csv')
S80 = S80.fillna(0)
S79['Q_cmd'] = S79['S79_TOT_FLOW_cfs'] * 0.0283168466 * 86400
S80['Q_cmd'] = S80['S80_S_FLOW_cfs'] * 0.0283168466 * 86400
S79 = DF_Date_Range(S79, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
S80 = DF_Date_Range(S80, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
C43RO_df = pd.DataFrame(S79['date'], columns=['date'])
C44RO_df = pd.DataFrame(S79['date'], columns=['date'])
C43RO = np.zeros(len(C43RO_df.index))
C44RO = np.zeros(len(C44RO_df.index))
for i in range(len(C44RO_df.index)):
    if S79['Q_cmd'].iloc[i] - Flow_df['S77_Out'].iloc[i] + Flow_df['S77_In'].iloc[i] < 0:
        C43RO[i] = 0
    else:
        C43RO[i] = S79['Q_cmd'].iloc[i] - Flow_df['S77_Out'].iloc[i] + Flow_df['S77_In'].iloc[i]
for i in range(len(C44RO_df.index)):
    if S80['Q_cmd'].iloc[i] - Flow_df['S308_Out'].iloc[i] + Flow_df['S308_In'].iloc[i] < 0:
        C44RO[i] = 0
    else:
        C44RO[i] = S80['Q_cmd'].iloc[i] - Flow_df['S308_Out'].iloc[i] + Flow_df['S308_In'].iloc[i]
C43RO_df['C43RO_cmd'] = C43RO
C44RO_df['C44RO_cmd'] = C44RO
C43RO_df['C43RO_cfs'] = C43RO_df['C43RO_cmd']/(0.0283168466 * 86400)
C44RO_df['C44RO_cfs'] = C44RO_df['C44RO_cmd']/(0.0283168466 * 86400)
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
C43RO_df.to_csv('./C43RO.csv')
C44RO_df.to_csv('./C44RO.csv')
C43RO_df = C43RO_df.set_index(C43RO_df['date'])
C44RO_df = C44RO_df.set_index(C44RO_df['date'])
C43RO_df.index = pd.to_datetime(C43RO_df.index, unit = 'ns')
C44RO_df.index = pd.to_datetime(C44RO_df.index, unit = 'ns')
C43Mon = C43RO_df.resample('M').mean()
C44Mon = C44RO_df.resample('M').mean()
C43Mon.to_csv('./C43RO_Monthly_LORS20082023.csv')
C44Mon.to_csv('./C44RO_Monthly_LORS20082023.csv')

#SLTRIB
#Create File (SLTRIB_Monthly_LORS20082023)
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
S48_S = pd.read_csv('./S48_S.csv')
S49_S = pd.read_csv('./S49_S.csv')
S48_S = DF_Date_Range(S48_S, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
S49_S = DF_Date_Range(S49_S, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
SLTRIB = pd.DataFrame(S48_S['date'],columns=['date'])
SLTRIB['SLTRIB_cmd'] = S48_S['S48_S_FLOW_cmd'] + S49_S['S49_S_FLOW_cmd']  
SLTRIB['SLTRIB_cfs'] = SLTRIB['SLTRIB_cmd']/(0.0283168466 * 86400)
SLTRIB = SLTRIB.set_index(SLTRIB['date'])
SLTRIB.index = pd.to_datetime(SLTRIB.index, unit = 'ns')
SLTRIBMon = SLTRIB.resample('M').mean()
SLTRIB = SLTRIB.reset_index()
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
SLTRIB.to_csv('./SLTRIB.csv')
SLTRIBMon.to_csv('./SLTRIB_Monthly_LORS20082023.csv')
Basin_RO = pd.DataFrame(SLTRIBMon.index,columns=['date'])
Basin_RO['SLTRIB'] = SLTRIBMon['SLTRIB_cfs'].values * 1.9835 #cfs to acft
Basin_RO['C44RO'] = C44Mon['C44RO_cfs'].values * 86400
Basin_RO['C43RO'] = C43Mon['C43RO_cfs'].values * 86400
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
Basin_RO.to_csv('./Basin_RO_inputs_LORS20082023.csv')

#EAA MIA RUNOFF
#Create File (EAA_MIA_RUNOFF_Inputs_LORS20082023)
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
S3_Miami_data = pd.read_csv('./S3_Miami.csv')
S3_Miami_data = DF_Date_Range(S3_Miami_data, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
S3_Miami = S3_Miami_data['S3_FLOW_cfs']
S2_NNR_data = pd.read_csv('./S2_NNR.csv')
S2_NNR_data = DF_Date_Range(S2_NNR_data, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
S2_NNR = S2_NNR_data['S2 NNR_FLOW_cfs']
EAA_MIA_RO = pd.DataFrame(date,columns=['date'])
EAA_MIA_RO['MIA'] = S3_Miami.values
EAA_MIA_RO['NNR'] = S2_NNR.values
EAA_MIA_RO['WPB'] = Flow_df['S352_Out']/(0.0283168466 * 86400)
EAA_MIA_RO['S2PMP'] = Flow_df['S2_P_Q']/(0.0283168466 * 86400)
EAA_MIA_RO['S3PMP'] = Flow_df['S3_P_Q']/(0.0283168466 * 86400)
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
EAA_MIA_RO.to_csv('./EAA_MIA_RUNOFF_Inputs_LORS20082023.csv')

#Weekly Tributary Conditions
#Create File (Trib_cond_wkly_data_LORS20082023)
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
#Net RF Inch
RF_data = pd.read_csv('./LAKE_RAINFALL_DATA_2008-2023.csv')
ET_data = pd.read_csv('./LOONE_AVERAGE_ETPI_DATA_2008-2023.csv')
Net_RF = pd.DataFrame(RF_data['date'], columns=['date'])
Net_RF['NetRF_In'] = RF_data['Avg_RF_In'] - ET_data['Avg_ET_In']
Net_RF = Net_RF.set_index(['date'])
Net_RF.index = pd.to_datetime(Net_RF.index, unit = 'ns')
Net_RF_Weekly =  Net_RF.resample('W-FRI').sum()
#Net Inflows cfs
Net_Inflows = pd.DataFrame(Flow_df['date'], columns=['date'])
Net_Inflows = DF_Date_Range(Net_Inflows, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
Net_Inflows['Net_Inflows'] = Flow_df['Netflows']/(0.0283168466 * 86400) #cmd to cfs
Net_Inflows = Net_Inflows.set_index(['date'])
Net_Inflows.index = pd.to_datetime(Net_Inflows.index, unit = 'ns')
Net_Inflow_Weekly = Net_Inflows.resample('W-FRI').mean()
# S65 cfs
S65E = pd.DataFrame(Flow_df['date'], columns=['date'])
S65E = DF_Date_Range(S65E, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
S65E['S65E'] = Flow_df['S65_Q']/(0.0283168466 * 86400) #cmd to cfs
S65E = S65E.set_index(['date'])
S65E.index = pd.to_datetime(S65E.index, unit = 'ns')
S65E_Weekly = S65E.resample('W-FRI').mean()
# PI 
#TODO
#This is prepared manually 
#Weekly data is downloaded from https://www.ncei.noaa.gov/access/monitoring/weekly-palmers/time-series/0804
#State:Florida Division:4.South Central
PI = pd.DataFrame(S65E_Weekly.index,columns=['date'])
PI_data = pd.read_csv('./PI_2008-2023.csv')
PI['PI'] = PI_data['PI']

Trib_Cond_Wkly = pd.DataFrame(S65E_Weekly.index,columns=['date'])
Trib_Cond_Wkly['NetRF'] = Net_RF_Weekly['NetRF_In'].values
Trib_Cond_Wkly['NetInf'] = Net_Inflow_Weekly['Net_Inflows'].values
Trib_Cond_Wkly['S65E'] = S65E_Weekly['S65E'].values
Trib_Cond_Wkly['Palmer'] = PI['PI'].values
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
Trib_Cond_Wkly.to_csv('./Trib_cond_wkly_data_LORS20082023.csv')


#Wind Speed
#Create File (LOWS)
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
L001WS = pd.read_csv('./L001_WNDS_MPH.csv')
L005WS = pd.read_csv('./L005_WNDS_MPH.csv')
L006WS = pd.read_csv('./L006_WNDS_MPH.csv')
LZ40WS = pd.read_csv('./LZ40_WNDS_MPH.csv')
L001WS = DF_Date_Range(L001WS, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
L005WS = DF_Date_Range(L005WS, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
L006WS = DF_Date_Range(L006WS, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
LZ40WS = DF_Date_Range(LZ40WS, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
LOWS = pd.DataFrame(L001WS['date'], columns=['date'])
LOWS['L001WS'] = L001WS['L001_WNDS_MPH']
LOWS['L005WS'] = L005WS['L005_WNDS_MPH']
LOWS['L006WS'] = L006WS['L006_WNDS_MPH']
LOWS['LZ40WS'] = LZ40WS['LZ40_WNDS_MPH']
LOWS['LO_Avg_WS_MPH'] = LOWS.mean(axis=1)
LOWS.to_csv('./LOWS.csv')

#RFVol acft
#Create File (RF_Volume_LORS20082023)
RFVol = pd.DataFrame(RF_data['date'],columns=['date'])
RFVol['RFVol_acft'] = (RF_data['Avg_RF_In'].values/12) * LO_Stg_Sto_SA_df['SA_acres'].values
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
RFVol.to_csv('./RFVol_LORS_20082023.csv')

#ETVol acft
#Create File (ETVol_LORS20082023)
ETVol = pd.DataFrame(ET_data['date'],columns=['date'])
ETVol['ETVol_acft'] = (ET_data['Avg_ET_In'].values/12) * LO_Stg_Sto_SA_df['SA_acres'].values
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
ETVol.to_csv('./ETVol_LORS_20082023.csv')


#WCA Stages
#Create File (WCA_Stages_Inputs_LORS20082023)
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
Stg_3ANW = pd.read_csv('./Stg_3ANW.csv')
Stg_3ANW = DF_Date_Range(Stg_3ANW, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
Stg_2A17 = pd.read_csv('./Stg_2A17.csv')
Stg_2A17 = DF_Date_Range(Stg_2A17, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
Stg_3A3 = pd.read_csv('./Stg_3A3.csv')
Stg_3A3 = DF_Date_Range(Stg_3A3, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
Stg_3A4 = pd.read_csv('./Stg_3A4.csv')
Stg_3A4 = DF_Date_Range(Stg_3A4, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
Stg_3A28 = pd.read_csv('./Stg_3A28.csv')
Stg_3A28 = DF_Date_Range(Stg_3A28, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
WCA_Stg = pd.DataFrame(Stg_3A28['date'],columns=['date'])
WCA_Stg['3A-NW'] = Stg_3ANW['3A-NW_STG_ft NGVD29'].values
WCA_Stg['2A-17'] = Stg_2A17['2-17_GAGHT_feet'].values
WCA_Stg['3A-3'] = Stg_3A3['3-63_GAGHT_feet'].values
WCA_Stg['3A-4'] = Stg_3A4['3-64_GAGHT_feet'].values
WCA_Stg['3A-28'] = Stg_3A28['3-65_GAGHT_feet'].values
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
WCA_Stg.to_csv('./WCA_Stages_Inputs_LORS20082023.csv')


# Predict Water Temp Function of Air Temp
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
Water_Temp_data = pd.read_csv('./Temp_Avg.csv')
L001_AirT = pd.read_csv('./L001_AirT.csv')
L001_AirT = DF_Date_Range(L001_AirT, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
L005_AirT = pd.read_csv('./L005_AirT.csv')
L005_AirT = DF_Date_Range(L005_AirT, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
L006_AirT = pd.read_csv('./L006_AirT.csv')
L006_AirT = DF_Date_Range(L006_AirT, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
LZ40_AirT = pd.read_csv('./LZ40_AirT.csv')
LZ40_AirT = DF_Date_Range(LZ40_AirT, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 

Water_Temp_data['L001_WaterT'] = Water_Temp_data[['L001_H2OT_C_1', 'L001_H2OT_C_2', 'L001_H2OT_C_3']].mean(axis=1)
Water_Temp_data['L005_WaterT'] = Water_Temp_data[['L005_H2OT_C_1', 'L005_H2OT_C_2', 'L005_H2OT_C_3']].mean(axis=1)
Water_Temp_data['L006_WaterT'] = Water_Temp_data[['L006_H2OT_C_1', 'L006_H2OT_C_2', 'L006_H2OT_C_3']].mean(axis=1)
Water_Temp_data['LZ40_WaterT'] = Water_Temp_data[['LZ40_H2OT_C_1', 'LZ40_H2OT_C_2', 'LZ40_H2OT_C_3']].mean(axis=1)

Water_Temp_data = DF_Date_Range(Water_Temp_data, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 

WaterT_pred_df = pd.DataFrame(L001_AirT['date'], columns=['date'])
WaterT_pred_df['L001_WaterT_pred'] = 1.862667 + 0.936899 * L001_AirT['L001_AIRT_Degrees Celsius'].values
WaterT_pred_df['L005_WaterT_pred'] = 1.330211 + 0.909713 * L005_AirT['L005_AIRT_Degrees Celsius'].values
WaterT_pred_df['L006_WaterT_pred'] = -0.88564 + 1.01585 * L006_AirT['L006_AIRT_Degrees Celsius'].values
WaterT_pred_df['LZ40_WaterT_pred'] = 0.388231 + 0.980154 * LZ40_AirT['LZ40_AIRT_Degrees Celsius'].values
WaterT_pred_df['WaterT_pred_Mean'] = WaterT_pred_df[['L001_WaterT_pred','L005_WaterT_pred','L006_WaterT_pred','LZ40_WaterT_pred']].mean(axis=1)
WaterT_pred_df_1 = DF_Date_Range(WaterT_pred_df, St_Yr, St_M, St_D, 2020, 8, 25) 
WaterT_pred_df_2 = DF_Date_Range(WaterT_pred_df, 2020, 8, 26, En_Yr, En_M, En_D) 
Filled_WaterT_1 = np.zeros(len(WaterT_pred_df_1.index))
Filled_WaterT_2 = np.zeros(len(WaterT_pred_df_2.index))
for i in range(len(Water_Temp_data.index)):
    if isnan(Water_Temp_data['WaterT_Mean'].iloc[i]):
        Filled_WaterT_1[i] = WaterT_pred_df_1['WaterT_pred_Mean'].iloc[i]
    else:
        Filled_WaterT_1[i] = Water_Temp_data['WaterT_Mean'].iloc[i]
        
Filled_WaterT_2 = WaterT_pred_df_2['WaterT_pred_Mean']
Filled_WaterT_1df = pd.DataFrame(WaterT_pred_df_1['date'],columns=['date'])
Filled_WaterT_2df = pd.DataFrame(WaterT_pred_df_2['date'],columns=['date'])
Filled_WaterT_1df['Water_T'] = Filled_WaterT_1
Filled_WaterT_2df['Water_T'] = Filled_WaterT_2
Filled_WaterT = pd.concat([Filled_WaterT_1df, Filled_WaterT_2df]).reset_index(drop= True)
Filled_WaterT.to_csv('./Filled_WaterT_20082023.csv')      

# TP Observations in Lake 
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_TP = pd.read_csv('./water_quality_L001_PHOSPHATE, TOTAL AS P.csv')
L004_TP = pd.read_csv('./water_quality_L004_PHOSPHATE, TOTAL AS P.csv')
L005_TP = pd.read_csv('./water_quality_L005_PHOSPHATE, TOTAL AS P.csv')
L006_TP = pd.read_csv('./water_quality_L006_PHOSPHATE, TOTAL AS P.csv')
L007_TP = pd.read_csv('./water_quality_L007_PHOSPHATE, TOTAL AS P.csv')
L008_TP = pd.read_csv('./water_quality_L008_PHOSPHATE, TOTAL AS P.csv')
LZ40_TP = pd.read_csv('./water_quality_LZ40_PHOSPHATE, TOTAL AS P.csv')

LO_TP_data = pd.merge(L001_TP,L004_TP, how = 'left', on = 'date')
LO_TP_data = pd.merge(LO_TP_data,L005_TP, how = 'left', on = 'date')
LO_TP_data = pd.merge(LO_TP_data,L006_TP, how = 'left', on = 'date')
LO_TP_data = pd.merge(LO_TP_data,L007_TP, how = 'left', on = 'date')
LO_TP_data = pd.merge(LO_TP_data,L008_TP, how = 'left', on = 'date')
LO_TP_data = pd.merge(LO_TP_data,LZ40_TP, how = 'left', on = 'date')
LO_TP_data = LO_TP_data.loc[:,~LO_TP_data.columns.str.startswith('Unnamed')]
LO_TP_data['Mean_TP'] = LO_TP_data.mean(axis=1)
LO_TP_data = LO_TP_data.set_index(['date'])
LO_TP_data.index = pd.to_datetime(LO_TP_data.index, unit = 'ns')
LO_TP_Monthly = LO_TP_data.resample('M').mean()
LO_TP_Monthly.to_csv('./LO_TP_Monthly.csv')

# Interpolated TP Observations in Lake 
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_TP_Inter = pd.read_csv('./water_quality_L001_PHOSPHATE, TOTAL AS P_Interpolated.csv')
L004_TP_Inter = pd.read_csv('./water_quality_L004_PHOSPHATE, TOTAL AS P_Interpolated.csv')
L005_TP_Inter = pd.read_csv('./water_quality_L005_PHOSPHATE, TOTAL AS P_Interpolated.csv')
L006_TP_Inter = pd.read_csv('./water_quality_L006_PHOSPHATE, TOTAL AS P_Interpolated.csv')
L007_TP_Inter = pd.read_csv('./water_quality_L007_PHOSPHATE, TOTAL AS P_Interpolated.csv')
L008_TP_Inter = pd.read_csv('./water_quality_L008_PHOSPHATE, TOTAL AS P_Interpolated.csv')
LZ40_TP_Inter = pd.read_csv('./water_quality_LZ40_PHOSPHATE, TOTAL AS P_Interpolated.csv')

LO_TP_data_Inter = pd.merge(L001_TP_Inter,L004_TP_Inter, how = 'left', on = 'date')
LO_TP_data_Inter = pd.merge(LO_TP_data_Inter,L005_TP_Inter, how = 'left', on = 'date')
LO_TP_data_Inter = pd.merge(LO_TP_data_Inter,L006_TP_Inter, how = 'left', on = 'date')
LO_TP_data_Inter = pd.merge(LO_TP_data_Inter,L007_TP_Inter, how = 'left', on = 'date')
LO_TP_data_Inter = pd.merge(LO_TP_data_Inter,L008_TP_Inter, how = 'left', on = 'date')
LO_TP_data_Inter = pd.merge(LO_TP_data_Inter,LZ40_TP_Inter, how = 'left', on = 'date')
LO_TP_data_Inter = LO_TP_data_Inter.loc[:,~LO_TP_data_Inter.columns.str.startswith('Unnamed')]
LO_TP_data_Inter['Mean_TP'] = LO_TP_data_Inter.mean(axis=1)
LO_TP_data_Inter = LO_TP_data_Inter.set_index(['date'])
LO_TP_data_Inter.index = pd.to_datetime(LO_TP_data_Inter.index, unit = 'ns')
LO_TP_Monthly_Inter = LO_TP_data_Inter.resample('M').mean()
Max = LO_TP_Monthly_Inter.max(axis=1)
Min = LO_TP_Monthly_Inter.min(axis=1)
LO_TP_Monthly_Inter['Max'] = Max.values
LO_TP_Monthly_Inter['Min'] = Min.values
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
LO_TP_Monthly_Inter.to_csv('./LO_TP_Monthly.csv')

# Interpolated OP Observations in Lake 
# Create File (LO_Avg_OP_2008-2022)
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_OP_Inter = pd.read_csv('./water_quality_L001_PHOSPHATE, ORTHO AS P_Interpolated.csv')
L004_OP_Inter = pd.read_csv('./water_quality_L004_PHOSPHATE, ORTHO AS P_Interpolated.csv')
L005_OP_Inter = pd.read_csv('./water_quality_L005_PHOSPHATE, ORTHO AS P_Interpolated.csv')
L006_OP_Inter = pd.read_csv('./water_quality_L006_PHOSPHATE, ORTHO AS P_Interpolated.csv')
L007_OP_Inter = pd.read_csv('./water_quality_L007_PHOSPHATE, ORTHO AS P_Interpolated.csv')
L008_OP_Inter = pd.read_csv('./water_quality_L008_PHOSPHATE, ORTHO AS P_Interpolated.csv')
LZ40_OP_Inter = pd.read_csv('./water_quality_LZ40_PHOSPHATE, ORTHO AS P_Interpolated.csv')

LO_OP_data_Inter = pd.merge(L001_OP_Inter,L004_OP_Inter, how = 'left', on = 'date')
LO_OP_data_Inter = pd.merge(LO_OP_data_Inter,L005_OP_Inter, how = 'left', on = 'date')
LO_OP_data_Inter = pd.merge(LO_OP_data_Inter,L006_OP_Inter, how = 'left', on = 'date')
LO_OP_data_Inter = pd.merge(LO_OP_data_Inter,L007_OP_Inter, how = 'left', on = 'date')
LO_OP_data_Inter = pd.merge(LO_OP_data_Inter,L008_OP_Inter, how = 'left', on = 'date')
LO_OP_data_Inter = pd.merge(LO_OP_data_Inter,LZ40_OP_Inter, how = 'left', on = 'date')
LO_OP_data_Inter = LO_OP_data_Inter.loc[:,~LO_OP_data_Inter.columns.str.startswith('Unnamed')]
LO_OP_data_Inter['Mean_OP'] = LO_OP_data_Inter.mean(axis=1)
LO_OP_data_Inter = DF_Date_Range(LO_OP_data_Inter, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
LO_OP_data_Inter.to_csv('./LO_OP.csv')


# Interpolated NH4 Observations in Lake
#Create File (LO_Avg_NH4_2008-2022) 
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_NH4_Inter = pd.read_csv('./water_quality_L001_AMMONIA-N_Interpolated.csv')
L004_NH4_Inter = pd.read_csv('./water_quality_L004_AMMONIA-N_Interpolated.csv')
L005_NH4_Inter = pd.read_csv('./water_quality_L005_AMMONIA-N_Interpolated.csv')
L006_NH4_Inter = pd.read_csv('./water_quality_L006_AMMONIA-N_Interpolated.csv')
L007_NH4_Inter = pd.read_csv('./water_quality_L007_AMMONIA-N_Interpolated.csv')
L008_NH4_Inter = pd.read_csv('./water_quality_L008_AMMONIA-N_Interpolated.csv')
LZ40_NH4_Inter = pd.read_csv('./water_quality_LZ40_AMMONIA-N_Interpolated.csv')

LO_NH4_data_Inter = pd.merge(L001_NH4_Inter,L004_NH4_Inter, how = 'left', on = 'date')
LO_NH4_data_Inter = pd.merge(LO_NH4_data_Inter,L005_NH4_Inter, how = 'left', on = 'date')
LO_NH4_data_Inter = pd.merge(LO_NH4_data_Inter,L006_NH4_Inter, how = 'left', on = 'date')
LO_NH4_data_Inter = pd.merge(LO_NH4_data_Inter,L007_NH4_Inter, how = 'left', on = 'date')
LO_NH4_data_Inter = pd.merge(LO_NH4_data_Inter,L008_NH4_Inter, how = 'left', on = 'date')
LO_NH4_data_Inter = pd.merge(LO_NH4_data_Inter,LZ40_NH4_Inter, how = 'left', on = 'date')
LO_NH4_data_Inter.to_csv('./LO_NH4_Inter.csv')
#Read clean LO_NH4 data
LO_NH4_Clean_Inter = pd.read_csv('./LO_NH4_Inter.csv')
LO_NH4_Clean_Inter['Mean_NH4'] = LO_NH4_Clean_Inter.mean(axis=1)
LO_NH4_Clean_Inter.to_csv('./LO_NH4_Clean_daily.csv')
LO_NH4_Clean_Inter = LO_NH4_Clean_Inter.set_index(['date'])
LO_NH4_Clean_Inter.index = pd.to_datetime(LO_NH4_Clean_Inter.index, unit = 'ns')
LO_NH4_Monthly_Inter = LO_NH4_Clean_Inter.resample('M').mean()
LO_NH4_Monthly_Inter.to_csv('./LO_NH4_Monthly_Inter.csv')

# Interpolated NO Observations in Lake
#Create File (LO_Avg_NO_2008-2022) and (LO_NO_Obs20082022)
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_NO_Inter = pd.read_csv('./water_quality_L001_NITRATE+NITRITE-N_Interpolated.csv')
L004_NO_Inter = pd.read_csv('./water_quality_L004_NITRATE+NITRITE-N_Interpolated.csv')
L005_NO_Inter = pd.read_csv('./water_quality_L005_NITRATE+NITRITE-N_Interpolated.csv')
L006_NO_Inter = pd.read_csv('./water_quality_L006_NITRATE+NITRITE-N_Interpolated.csv')
L007_NO_Inter = pd.read_csv('./water_quality_L007_NITRATE+NITRITE-N_Interpolated.csv')
L008_NO_Inter = pd.read_csv('./water_quality_L008_NITRATE+NITRITE-N_Interpolated.csv')
LZ40_NO_Inter = pd.read_csv('./water_quality_LZ40_NITRATE+NITRITE-N_Interpolated.csv')

LO_NO_data_Inter = pd.merge(L001_NO_Inter,L004_NO_Inter, how = 'left', on = 'date')
LO_NO_data_Inter = pd.merge(LO_NO_data_Inter,L005_NO_Inter, how = 'left', on = 'date')
LO_NO_data_Inter = pd.merge(LO_NO_data_Inter,L006_NO_Inter, how = 'left', on = 'date')
LO_NO_data_Inter = pd.merge(LO_NO_data_Inter,L007_NO_Inter, how = 'left', on = 'date')
LO_NO_data_Inter = pd.merge(LO_NO_data_Inter,L008_NO_Inter, how = 'left', on = 'date')
LO_NO_data_Inter = pd.merge(LO_NO_data_Inter,LZ40_NO_Inter, how = 'left', on = 'date')
LO_NO_data_Inter = LO_NO_data_Inter.loc[:,~LO_NO_data_Inter.columns.str.startswith('Unnamed')]
LO_NO_data_Inter['Mean_NO'] = LO_NO_data_Inter.mean(axis=1)
# LO_NO_data_Inter.to_csv('./LO_NO_Clean_daily.csv')
LO_NO_data_Inter = LO_NO_data_Inter.set_index(['date'])
LO_NO_data_Inter.index = pd.to_datetime(LO_NO_data_Inter.index, unit = 'ns')
LO_NO_Monthly_Inter = LO_NO_data_Inter.resample('M').mean()
NO_Max = LO_NO_Monthly_Inter.max(axis=1)
NO_Min = LO_NO_Monthly_Inter.min(axis=1)
LO_NO_Monthly_Inter['Max'] = NO_Max.values
LO_NO_Monthly_Inter['Min'] = NO_Min.values
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
LO_NO_Monthly_Inter.to_csv('./LO_NO_Monthly_Inter.csv')

#Create File (LO_DIN_2008-2022)
date_DIN = pd.date_range(start = '%s/%s/%s'%(St_M,St_D,St_Yr),end = '%s/%s/%s'%(En_M,En_D,En_Yr),freq = 'D')
LO_DIN  = pd.DataFrame(date_DIN,columns=['date'])
LO_NH4_Clean_Inter = DF_Date_Range(LO_NH4_Clean_Inter, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
LO_NO_Clean_Inter = DF_Date_Range(LO_NO_Clean_Inter, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
LO_DIN['NH4'] = LO_NH4_Clean_Inter['Mean_NH4'].values 
LO_DIN['NO'] = LO_NO_Clean_Inter['Mean_NO'].values 
LO_DIN['DIN_mg/m3'] = LO_DIN[['NH4','NO']].sum(axis=1)*1000
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
LO_DIN.to_csv('./LO_DIN.csv')

# Interpolated DO Observations in Lake 
#Create File (LO_Avg_DO_2008-2022)
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_DO_Inter = pd.read_csv('./water_quality_L001_DISSOLVED OXYGEN_Interpolated.csv')
L004_DO_Inter = pd.read_csv('./water_quality_L004_DISSOLVED OXYGEN_Interpolated.csv')
L005_DO_Inter = pd.read_csv('./water_quality_L005_DISSOLVED OXYGEN_Interpolated.csv')
L006_DO_Inter = pd.read_csv('./water_quality_L006_DISSOLVED OXYGEN_Interpolated.csv')
L007_DO_Inter = pd.read_csv('./water_quality_L007_DISSOLVED OXYGEN_Interpolated.csv')
L008_DO_Inter = pd.read_csv('./water_quality_L008_DISSOLVED OXYGEN_Interpolated.csv')
LZ40_DO_Inter = pd.read_csv('./water_quality_LZ40_DISSOLVED OXYGEN_Interpolated.csv')

LO_DO_data_Inter = pd.merge(L001_DO_Inter,L004_DO_Inter, how = 'left', on = 'date')
LO_DO_data_Inter = pd.merge(LO_DO_data_Inter,L005_DO_Inter, how = 'left', on = 'date')
LO_DO_data_Inter = pd.merge(LO_DO_data_Inter,L006_DO_Inter, how = 'left', on = 'date')
LO_DO_data_Inter = pd.merge(LO_DO_data_Inter,L007_DO_Inter, how = 'left', on = 'date')
LO_DO_data_Inter = pd.merge(LO_DO_data_Inter,L008_DO_Inter, how = 'left', on = 'date')
LO_DO_data_Inter = pd.merge(LO_DO_data_Inter,LZ40_DO_Inter, how = 'left', on = 'date')
LO_DO_data_Inter = LO_DO_data_Inter.loc[:,~LO_DO_data_Inter.columns.str.startswith('Unnamed')]
#Read clean LO_DO data
LO_DO_data_Inter['Mean_DO'] = LO_DO_data_Inter.mean(axis=1)
LO_DO_data_Inter = DF_Date_Range(LO_DO_data_Inter, St_Yr, St_M, St_D, En_Yr, En_M, En_D)
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
LO_DO_data_Inter.to_csv('./LO_DO_Clean_daily.csv')
LO_DO_data_Inter = LO_DO_data_Inter.set_index(['date'])
LO_DO_data_Inter.index = pd.to_datetime(LO_DO_data_Inter.index, unit = 'ns')
LO_DO_Monthly_Inter = LO_DO_data_Inter.resample('M').mean()
LO_DO_Monthly_Inter.to_csv('./LO_DO_Monthly_Inter.csv')


#RADT Data in Lake Okeechobee
#Create File (LO_RADT_2008-2022)
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_Weather_Data/Weather_data_2023'
os.chdir('%s'%Working_dir) 
L001_RADT = pd.read_csv('./L001_RADT.csv')
L005_RADT = pd.read_csv('./L005_RADT.csv')
L006_RADT = pd.read_csv('./L006_RADT.csv')
LZ40_RADT = pd.read_csv('./LZ40_RADT.csv')
LO_RADT_data = pd.merge(L006_RADT,L001_RADT, how = 'left', on = 'date')
LO_RADT_data = pd.merge(LO_RADT_data,L005_RADT, how = 'left', on = 'date')
LO_RADT_data = pd.merge(LO_RADT_data,LZ40_RADT, how = 'left', on = 'date')
LO_RADT_data = LO_RADT_data.loc[:,~LO_RADT_data.columns.str.startswith('Unnamed')]
LO_RADT_data['Mean_RADT'] = LO_RADT_data.mean(axis=1)
LO_RADT_data = DF_Date_Range(LO_RADT_data, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
LO_RADT_data.to_csv('./LO_RADT_data_20082022.csv')

#RADP Data in Lake Okeechobee
#Create File (LO_RADP_2008-2022)
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_Weather_Data/Weather_data_2023'
os.chdir('%s'%Working_dir) 
L001_RADP = pd.read_csv('./L001_RADP.csv')
L005_RADP = pd.read_csv('./L005_RADP.csv')
L006_RADP = pd.read_csv('./L006_RADP.csv')
LZ40_RADP = pd.read_csv('./LZ40_RADP.csv')
LO_RADP_data = pd.merge(L006_RADP,L001_RADP, how = 'left', on = 'date')
LO_RADP_data = pd.merge(LO_RADP_data,L005_RADP, how = 'left', on = 'date')
LO_RADP_data = pd.merge(LO_RADP_data,LZ40_RADP, how = 'left', on = 'date')
LO_RADP_data = LO_RADP_data.loc[:,~LO_RADP_data.columns.str.startswith('Unnamed')]
LO_RADP_data['Mean_RADP'] = LO_RADP_data.mean(axis=1)
LO_RADP_data = DF_Date_Range(LO_RADP_data, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
LO_RADP_data.to_csv('./LO_RADP_data_20082022.csv')

# Interpolated Chla Corrected Observations in Lake 
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_Chla_Inter = pd.read_csv('./water_quality_L001_CHLOROPHYLL-A, CORRECTED_Interpolated.csv')
L004_Chla_Inter = pd.read_csv('./water_quality_L004_CHLOROPHYLL-A, CORRECTED_Interpolated.csv')
L005_Chla_Inter = pd.read_csv('./water_quality_L005_CHLOROPHYLL-A, CORRECTED_Interpolated.csv')
L006_Chla_Inter = pd.read_csv('./water_quality_L006_CHLOROPHYLL-A, CORRECTED_Interpolated.csv')
L007_Chla_Inter = pd.read_csv('./water_quality_L007_CHLOROPHYLL-A, CORRECTED_Interpolated.csv')
L008_Chla_Inter = pd.read_csv('./water_quality_L008_CHLOROPHYLL-A, CORRECTED_Interpolated.csv')
LZ40_Chla_Inter = pd.read_csv('./water_quality_LZ40_CHLOROPHYLL-A, CORRECTED_Interpolated.csv')

LO_Chla_data_Inter = pd.merge(L001_Chla_Inter,L004_Chla_Inter, how = 'left', on = 'date')
LO_Chla_data_Inter = pd.merge(LO_Chla_data_Inter,L005_Chla_Inter, how = 'left', on = 'date')
LO_Chla_data_Inter = pd.merge(LO_Chla_data_Inter,L006_Chla_Inter, how = 'left', on = 'date')
LO_Chla_data_Inter = pd.merge(LO_Chla_data_Inter,L007_Chla_Inter, how = 'left', on = 'date')
LO_Chla_data_Inter = pd.merge(LO_Chla_data_Inter,L008_Chla_Inter, how = 'left', on = 'date')
LO_Chla_data_Inter = pd.merge(LO_Chla_data_Inter,LZ40_Chla_Inter, how = 'left', on = 'date')
LO_Chla_data_Inter = LO_Chla_data_Inter.loc[:,~LO_Chla_data_Inter.columns.str.startswith('Unnamed')]
#Read clean LO_Chla data
LO_Chla_data_Inter['Mean_Chla'] = LO_Chla_data_Inter.mean(axis=1)
LO_Chla_data_Inter.to_csv('./LO_Chla_Clean_daily.csv')
#Monthly
LO_Chla_data_Inter = LO_Chla_data_Inter.set_index(['date'])
LO_Chla_data_Inter.index = pd.to_datetime(LO_Chla_data_Inter.index, unit = 'ns')
LO_Chla_Monthly_Inter = LO_Chla_data_Inter.resample('M').mean()
LO_Chla_Monthly_Inter.to_csv('./LO_Chla_Monthly_Inter.csv')

# Interpolated Chla LC Observations in Lake 
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/LO_WQ_June2023'
os.chdir('%s'%Working_dir) 
L001_Chla_LC_Inter = pd.read_csv('./water_quality_L001_CHLOROPHYLL-A(LC)_Interpolated.csv')
L004_Chla_LC_Inter = pd.read_csv('./water_quality_L004_CHLOROPHYLL-A(LC)_Interpolated.csv')
L005_Chla_LC_Inter = pd.read_csv('./water_quality_L005_CHLOROPHYLL-A(LC)_Interpolated.csv')
L006_Chla_LC_Inter = pd.read_csv('./water_quality_L006_CHLOROPHYLL-A(LC)_Interpolated.csv')
L007_Chla_LC_Inter = pd.read_csv('./water_quality_L007_CHLOROPHYLL-A(LC)_Interpolated.csv')
L008_Chla_LC_Inter = pd.read_csv('./water_quality_L008_CHLOROPHYLL-A(LC)_Interpolated.csv')
LZ40_Chla_LC_Inter = pd.read_csv('./water_quality_LZ40_CHLOROPHYLL-A(LC)_Interpolated.csv')

LO_Chla_LC_data_Inter = pd.merge(L001_Chla_LC_Inter,L004_Chla_LC_Inter, how = 'left', on = 'date')
LO_Chla_LC_data_Inter = pd.merge(LO_Chla_LC_data_Inter,L005_Chla_LC_Inter, how = 'left', on = 'date')
LO_Chla_LC_data_Inter = pd.merge(LO_Chla_LC_data_Inter,L006_Chla_LC_Inter, how = 'left', on = 'date')
LO_Chla_LC_data_Inter = pd.merge(LO_Chla_LC_data_Inter,L007_Chla_LC_Inter, how = 'left', on = 'date')
LO_Chla_LC_data_Inter = pd.merge(LO_Chla_LC_data_Inter,L008_Chla_LC_Inter, how = 'left', on = 'date')
LO_Chla_LC_data_Inter = pd.merge(LO_Chla_LC_data_Inter,LZ40_Chla_LC_Inter, how = 'left', on = 'date')
LO_Chla_LC_data_Inter = LO_Chla_LC_data_Inter.loc[:,~LO_Chla_LC_data_Inter.columns.str.startswith('Unnamed')]
#Read clean LO_Chla_LC data
LO_Chla_LC_data_Inter['Mean_Chla_LC'] = LO_Chla_LC_data_Inter.mean(axis=1)
LO_Chla_LC_data_Inter.to_csv('./LO_Chla_LC_Clean_daily.csv')
#Monthly
LO_Chla_LC_data_Inter = LO_Chla_LC_data_Inter.set_index(['date'])
LO_Chla_LC_data_Inter.index = pd.to_datetime(LO_Chla_LC_data_Inter.index, unit = 'ns')
LO_Chla_LC_Monthly_Inter = LO_Chla_LC_data_Inter.resample('M').mean()
LO_Chla_LC_Monthly_Inter.to_csv('./LO_Chla_LC_Monthly_Inter.csv')

#Merge the Chla Data
#Create Files LO_Avg_Chla_2008-2022 and Obs_Chla_LO_2008-2022
# Chla_date = pd.date_range(start = LO_Chla_data_Inter['date'].iloc[0],end =LO_Chla_LC_data_Inter['date'].iloc[-1],freq = 'D')
LO_Chla_data_Inter = DF_Date_Range(LO_Chla_data_Inter, St_Yr, St_M, St_D, 2010, 10, 19) 
LO_Chla_df = pd.DataFrame(LO_Chla_data_Inter['date'],columns=['date'])
LO_Chla_df['Chla'] = LO_Chla_data_Inter['Mean_Chla']
LO_Chla_LC_df = pd.DataFrame(LO_Chla_LC_data_Inter['date'],columns=['date'])
LO_Chla_LC_df['Chla'] = LO_Chla_LC_data_Inter['Mean_Chla_LC']

LO_Chla_Merge = pd.concat([LO_Chla_df, LO_Chla_LC_df]).reset_index(drop= True)
LO_Chla_Merge.to_csv('./LO_Merged_Chla.csv')

LO_Chla_Merge = LO_Chla_Merge.set_index(['date'])
LO_Chla_Merge.index = pd.to_datetime(LO_Chla_Merge.index, unit = 'ns')
LO_Chla_Merge_Monthly_Inter = LO_Chla_Merge.resample('M').mean()
LO_Chla_Merge_Monthly_Inter.to_csv('./LO_Chla_Monthly_Inter.csv')


#NO Loads 
#Create File (Daily_NOx_External_Loads_2008-2022)
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/WQ_Data_May2023' 
os.chdir('%s'%Working_dir) 
S65_NO = pd.read_csv('./S65E_NO_Interpolated.csv')
S71_NO = pd.read_csv('./S71_NO_Interpolated.csv')
S72_NO = pd.read_csv('./S72_NO_Interpolated.csv')
S84_NO = pd.read_csv('./S84_NO_Interpolated.csv')
S127_NO = pd.read_csv('./S127_NO_Interpolated.csv')
S133_NO = pd.read_csv('./S133_NO_Interpolated.csv')
S135_NO = pd.read_csv('./S135_NO_Interpolated.csv')
S154_NO = pd.read_csv('./S154_NO_Interpolated.csv')
S191_NO = pd.read_csv('./S191_NO_Interpolated.csv')
S308_NO = pd.read_csv('./S308C_NO_Interpolated.csv')
FISHP_NO = pd.read_csv('./FECSR78_NO_Interpolated.csv')
L8_NO = pd.read_csv('./CULV10A_NO_Interpolated.csv')
S4_NO = pd.read_csv('./S4_NO_Interpolated.csv')

NO_names = ['S65_NO','S71_NO','S72_NO','S84_NO','S127_NO','S133_NO','S135_NO','S154_NO','S191_NO','S308_NO','FISHP_NO','L8_NO','S4_NO']
NO_list = {'S65_NO':S65_NO,'S71_NO':S71_NO,'S72_NO':S72_NO,'S84_NO':S84_NO,'S127_NO':S127_NO,'S133_NO':S133_NO,'S135_NO':S135_NO,
           'S154_NO':S154_NO,'S191_NO':S191_NO,'S308_NO':S308_NO,'FISHP_NO':FISHP_NO,'L8_NO':L8_NO,'S4_NO':S4_NO}
date_NO = pd.date_range(start = '1/1/2008',end ='12/31/2022',freq = 'D')

NO_df = pd.DataFrame(date_NO, columns=['date'])
for i in range(len(NO_names)):
    y = DF_Date_Range(NO_list[NO_names[i]], St_Yr, St_M, St_D, En_Yr, En_M, En_D)
    NO_df['%s'%NO_names[i]] = y.iloc[:,-1:].values

Flow_df = DF_Date_Range(Flow_df, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 

#Determine NO Loads
NO_Loads_In = pd.DataFrame(date_NO, columns=['date'])
NO_Loads_In['S65_NO_Ld'] = Flow_df['S65_Q'].values * NO_df['S65_NO'].values * 1000
NO_Loads_In['S71_NO_Ld'] = Flow_df['S71_Q'].values * NO_df['S71_NO'].values * 1000
NO_Loads_In['S72_NO_Ld'] = Flow_df['S72_Q'].values * NO_df['S72_NO'].values * 1000
NO_Loads_In['S84_NO_Ld'] = Flow_df['S84_Q'].values * NO_df['S84_NO'].values * 1000
NO_Loads_In['S127_NO_Ld'] = Flow_df['S127_In'].values * NO_df['S127_NO'].values * 1000
NO_Loads_In['S133_NO_Ld'] = Flow_df['S133_P_Q'].values * NO_df['S133_NO'].values * 1000
NO_Loads_In['S135_NO_Ld'] = Flow_df['S135_In'].values * NO_df['S135_NO'].values * 1000
NO_Loads_In['S154_NO_Ld'] = Flow_df['S154_Q'].values * NO_df['S154_NO'].values * 1000
NO_Loads_In['S191_NO_Ld'] = Flow_df['S191_Q'].values * NO_df['S191_NO'].values * 1000
NO_Loads_In['S308_NO_Ld'] = Flow_df['S308_In'].values * NO_df['S308_NO'].values * 1000
NO_Loads_In['FISHP_NO_Ld'] = Flow_df['FISHP_Q'].values * NO_df['FISHP_NO'].values * 1000
NO_Loads_In['L8_NO_Ld'] = Flow_df['L8_In'].values * NO_df['L8_NO'].values * 1000
NO_Loads_In['S4_NO_Ld'] = Flow_df['S4_P_Q'].values * NO_df['S4_NO'].values * 1000
#Calculate the total External Loads to Lake Okeechobee
NO_Loads_In['External_NO_Ld_mg'] = NO_Loads_In.sum(axis=1)
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
NO_Loads_In.to_csv('./LO_External_Loadings_NO.csv')

#Determine Chla Loads
#Create File (Chla_Loads_In_2008-2022)
Working_dir = 'C:/Work/Research/Data Analysis/Lake_O_water_Q_data/WQ_Data_May2023' 
os.chdir('%s'%Working_dir) 
S65E_Chla = pd.read_csv('./S65E_Chla_Merged.csv')
S65E_Chla = DF_Date_Range(S65E_Chla, St_Yr, St_M, St_D, En_Yr, En_M, En_D) 
Chla_Loads_In = pd.DataFrame(date_NO, columns=['date'])
Chla_Loads_In['Chla_Loads'] = Flow_df['Inflows'].values * S65E_Chla['Data'].values
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
Chla_Loads_In.to_csv('./Chla_Loads_In.csv')



#Write Data into csv files
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
#write Avg Stage (ft,m) Storage (acft, m3) SA (acres) to csv
LO_Stg_Sto_SA_df.to_csv('./Average_LO_Storage_3MLagLORS20082023.csv')
#Write S65 TP concentrations (mg/L)
S65_total_TP.to_csv('./S65_TP_3MLag.csv')
# TP External Loads 3 Months Lag (mg)
TP_Loads_In_3MLag_df.to_csv('./LO_External_Loadings_3MLag_LORS20082023.csv')
# Flow dataframe including Inflows, NetFlows, and Outflows (all in m3/day)
Flow_df.to_csv('./Flow_df_3MLag.csv')
#Inflows (cmd)
LO_Inflows_BK.to_csv('./LO_Inflows_BK_LORS20082023.csv')
#Outflows (cmd)
Outflows_consd.to_csv('./Outflows_consd_LORS20082023.csv')
# NetFlows (cmd)
Netflows.to_csv('./Netflows_acft_LORS20082023.csv')
#Total flows to WCAs (acft)
TotalQWCA.to_csv('./TotalQWCA_Obs_LORS20082023.csv')
# INDUST Outflows (cmd)
INDUST_Outflows.to_csv('./INDUST_Outflows.csv')
