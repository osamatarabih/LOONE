# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 18:44:37 2021

@author: osama
"""

import os
import pandas as pd
from datetime import datetime
import numpy as np
from calendar import monthrange  
from Model_Config import Model_Config 
Working_Path = Model_Config.Working_Path
os.chdir('%s'%Working_Path) 
from Pre_defined_Variables import Pre_defined_Variables 
from Model_variables import M_var
from LO_FNs import LO_FNs
from Stg_Sto_Ar import Stg_Sto_Ar 
from LONINO_FNs import LONINO_FNs
from Dec_Tree_FNs import Dec_Tree_FNs
from WCA_Stages_Cls import WCA_Stages_Cls
from Additional_Fncs import Add_Fn
from THC_Class import THC_Class
from df_WSMs import WSMs
from Data import Data
from Trib_HC import Trib_HC
from scipy.integrate import odeint
os.chdir('%s/Code/'%Working_Path) 
from LOONE_NCHLA_FNS import *
os.chdir('%s'%Working_Path) 


def LOONE_HydNO(): 
    print("LOONE Q Module is Running!")
    # Based on the defined Start and End year, month, and day on the Pre_defined_Variables File, Startdate and enddate are defined. 
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    startdate = datetime(year, month, day).date() 
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    begdateCS = datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.enddate_entry)
    enddate = datetime(year, month, day).date()
        
#############################################################################
    Results_data = pd.read_csv('./Outputs/Dec_Var_StLQ400_Temp25_EcoEnv_RegSouth_Hbound.csv')
    Results = Results_data['Value']
    P_1 = Results[0]
    P_2 = Results[1]
    NO_1 = Results[2]
    NO_2 = Results[3]
    Chla_1 = Results[4]
    Chla_2 = Results[5]
    S77_DV = Results[6:18]
    S77_DV = S77_DV.reset_index(drop=True)
    S308_DV = Results[18:30]
    S308_DV = S308_DV.reset_index(drop=True)
    # RegSouth_DV = Results[26:38]
    # RegSouth_DV = RegSouth_DV.reset_index(drop=True)

    if Model_Config.Sim_type == 0 or Model_Config.Sim_type == 1:
        WSMs()
        df_WSMs = pd.read_csv('./Data/df_WSMs.csv')
    else:
        df_WSMs = pd.read_csv('./Data/df_WSMs.csv')
    #The Following Code interpolates daily LOSA demand from weekly data for 6 differnet datasets where the user defines the LOSA demand that will be used based on a Code (1:6).
    #Set time frame for model run
    date_rng_2 = pd.date_range(start=startdate, end = enddate, freq= 'D')
    #Create a data frame with a date column
    Water_dmd = pd.DataFrame(date_rng_2, columns =['date'])
    
    N = []
    Wk = []
    #Generate a count list
    for i in Water_dmd['date']:
        if i.month == 1 and i.day == 1:
            n = 0
        else:
            n = n + 1
        N.append(n)
    Water_dmd['count'] = N
    #Calculate the week number for all rows in the data frame
    for i in Water_dmd['count']:
        if i > 363:
            J = 52
        else:
            J = int(i/7)+1
        Wk.append(J)
    Water_dmd['Week_num'] = Wk
    dd = [] #daily demand
    #Calculate daily water demand
    for i in Water_dmd['Week_num']:
        D = ((Data.Weekly_dmd['C%s'%Pre_defined_Variables.Code].iloc[i-1])/7)*(Pre_defined_Variables.Multiplier/100)
        dd.append(D)
    Water_dmd['Daily_demand'] = dd
    ##############################################################################################
    #Determine Tributary Hydrologic Conditions
    TC_LONINO_df = Trib_HC()
    #Determine WCA Stages
    WCA_Stages_df = WCA_Stages_Cls(TC_LONINO_df)
    #A dataframe to determine eachday's season (Months 11,12,1,2 are Season 1, Months 3,4,5 are season 2, Months 6,7 are season 3, Months 8,9,10 are season 4 )
    date_rng_5 = pd.date_range(start = startdate, end = enddate, freq ='D')
    Seasons = pd.DataFrame(date_rng_5, columns =['date'])
    Seas_Count = len(Seasons.index)
    for i in range(Seas_Count):
        if Seasons['date'].iloc[i].month > 2 and Seasons['date'].iloc[i].month < 6:
            S = 2
        elif Seasons['date'].iloc[i].month > 5 and Seasons['date'].iloc[i].month < 8:
            S = 3
        elif Seasons['date'].iloc[i].month > 7 and Seasons['date'].iloc[i].month < 11:
            S = 4
        else:
            S = 1
        M_var.Daily_Seasons[i] = S
        M_var.Mon[i] = Seasons['date'].iloc[i].month
    Seasons['Season'] = M_var.Daily_Seasons 
    Seasons['Month'] = M_var.Mon
##################################################################################################################   
    #This following Script runs the main model daily simulations.
    date_rng_6 = pd.date_range(start='12/30/%d'%(Pre_defined_Variables.startyear-1), end = enddate, freq= 'D')
    LO_Model = pd.DataFrame(date_rng_6, columns =['date'])
    LO_Model['Net_Inflow'] = Data.NetInf_Input['Netflows_acft']
    # Geoglows = pd.read_csv('C:/LOONE_WQ/Data/LORS2020/geoglows_flow_df_ens_07_predicted/geoglows_flow_df_ens_01_predicted.csv')
    # last_15_values = Geoglows['Netflows'].tail(15).values*50/1233.48
    # Replace the last 15 values in the target column in df1 with the values from df2
    # LO_Model['Net_Inflow'].iloc[-15:] = last_15_values

    n_rows = len(LO_Model.index)
    LO_Model['LOSA_dmd_SFWMM'] = Data.SFWMM_W_dmd['LOSA_dmd'] * (Pre_defined_Variables.Mult_LOSA/100)
    LO_Model['C44RO'] = Data.C44_Runoff['C44RO']
    ##################################
    DecTree_df = pd.DataFrame(date_rng_5, columns = ['Date'])
    DecTree_df['Zone_B_MetFcast'] = TC_LONINO_df['LONINO_Seasonal_Classes']
    #Create a dataframe that includes Monthly Mean Basin Runoff & BaseFlow-Runoff & Runoff-Baseflow (cfs)
    date_rng_11 = pd.date_range(start=startdate, end = enddate, freq= 'MS')
    date_rng_11d = pd.date_range(start=startdate, end = enddate, freq= 'D')
    date_rng_11d.name = 'Date'
    Basin_RO = pd.DataFrame(date_rng_11, columns =['date'])
    #Baseflows
    Outlet1_baseflow = Data.S77_RegRelRates['Zone_D0'].iloc[0]   
    Outlet2_baseflow = Data.S80_RegRelRates['Zone_D0'].iloc[0]
    #Calculta number of months in the timeseries data.
    num_B_R = len(Basin_RO.index)
    BS_C43RO = np.zeros(num_B_R)
    BS_C44RO = np.zeros(num_B_R)
    C44RO_SLTRIB = np.zeros(num_B_R)
    C44RO_BS = np.zeros(num_B_R)
    Num_days = np.zeros(num_B_R)
    for i in range(num_B_R) :
        Num_days[i] = monthrange(Basin_RO['date'].iloc[i].year, Basin_RO['date'].iloc[i].month)[1] #no. of days in each time step month.
        BS_C43RO[i] = max(0, (Outlet1_baseflow - Data.C43RO['C43RO'].iloc[i]))
        BS_C44RO[i] = max(0, (Outlet2_baseflow - Data.C44RO['C44RO'].iloc[i]))
        C44RO_SLTRIB[i] = BS_C44RO[i] + Data.SLTRIB['SLTRIB_cfs'].iloc[i]
        C44RO_BS[i] = max(0, Data.C44RO['C44RO'].iloc[i] - Outlet2_baseflow)*Num_days[i]
    Basin_RO['Ndays'] = Num_days
    Basin_RO['C43RO'] = Data.C43RO['C43RO']
    Basin_RO['BS-C43RO'] = BS_C43RO
    Basin_RO['C44RO'] = Data.C44RO['C44RO']
    Basin_RO['BS-C44RO'] = BS_C44RO
    Basin_RO['SLTRIB'] = Data.SLTRIB['SLTRIB_cfs']
    Basin_RO['C44RO_SLTRIB'] = C44RO_SLTRIB
    Basin_RO['C44RO-BS'] = C44RO_BS
    LO_Model['C43RO'] = Data.C43RO_Daily['C43RO']
    S80avgL1 = Data.Pulses['S-80_L1_%s'%Pre_defined_Variables.Schedule].mean()
    S80avgL2 = Data.Pulses['S-80_L2_%s'%Pre_defined_Variables.Schedule].mean()
    S80avgL3 = Data.Pulses['S-80_L3_%s'%Pre_defined_Variables.Schedule].mean()
    S77avgL1 = Data.Pulses['S-77_L1_%s'%Pre_defined_Variables.Schedule].mean() #LORS
    S77avgL2 = Data.Pulses['S-77_L2_%s'%Pre_defined_Variables.Schedule].mean() #LORS
    S77avgL3 = Data.Pulses['S-77_L3_%s'%Pre_defined_Variables.Schedule].mean()
    Basin_RO = Basin_RO.set_index(['date'])
    Basin_RO.index = pd.to_datetime(Basin_RO.index)
    Basin_RO_Daily = Basin_RO.reindex(date_rng_11d, method='ffill')
    Basin_RO = Basin_RO.reset_index()
    VLOOKUP1 = Basin_RO_Daily['BS-C44RO']
    VLOOKUP1_c = [x for x in VLOOKUP1 if ~np.isnan(x)]
    ##################################################################################################################
    #This following script contains the logic and calculations for the proposed Lake Okeechobee Adaptive Protocol.
    AdapProt_df = pd.DataFrame(date_rng_5, columns = ['date'])
    #Calculate Late Dry Season (Apr-May) logic.
    Late_Dry_Season = []
    for i in AdapProt_df['date']:
        if i.month > 3 and i.month < 6:
            L = True
        else:
            L= False
        Late_Dry_Season.append(L)
    AdapProt_df['Late_Dry_Season'] = Late_Dry_Season
    AdapProt_df['Tributary Hydrologic Condition'] = TC_LONINO_df['Tributary_Condition']       
    #Define "Low Chance" 6/1 stg<11'
    if Pre_defined_Variables.Opt_Date_Targ_Stg ==1:
        Targ_Stg = Data.Targ_Stg_June_1st
    else:
        Targ_Stg = Data.Targ_Stg_May_1st
    
    Targ_Stg_df = pd.DataFrame(date_rng_5, columns = ['dates'])
    for i in range(len(Targ_Stg_df)):
        M_var.V10per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,10,Targ_Stg)   
        M_var.V20per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,20,Targ_Stg)
        M_var.V25per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,25,Targ_Stg)
        M_var.V30per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,30,Targ_Stg)   
        M_var.V40per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,40,Targ_Stg) 
        M_var.V45per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,45,Targ_Stg)
        M_var.V50per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,50,Targ_Stg)   
        M_var.V60per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,60,Targ_Stg)   

    V10per_c = [x for x in M_var.V10per if ~np.isnan(x)]
    V20per_c = [x for x in M_var.V20per if ~np.isnan(x)]
    V25per_c = [x for x in M_var.V25per if ~np.isnan(x)]
    V30per_c = [x for x in M_var.V30per if ~np.isnan(x)]
    V40per_c = [x for x in M_var.V40per if ~np.isnan(x)]
    V45per_c = [x for x in M_var.V45per if ~np.isnan(x)]
    V50per_c = [x for x in M_var.V50per if ~np.isnan(x)]
    V60per_c = [x for x in M_var.V60per if ~np.isnan(x)]
    Targ_Stg_df['10%'] = V10per_c
    Targ_Stg_df['20%'] = V20per_c
    Targ_Stg_df['25%'] = V25per_c
    Targ_Stg_df['30%'] = V30per_c
    Targ_Stg_df['40%'] = V40per_c
    Targ_Stg_df['45%'] = V45per_c
    Targ_Stg_df['50%'] = V50per_c
    Targ_Stg_df['60%'] = V60per_c
    
    # Outlet1_baseflow = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    Outlet1_baseflow = 450 #cfs
    VLOOKUP2 = Basin_RO_Daily['BS-C43RO']
    VLOOKUP2_c = [x for x in VLOOKUP2 if ~np.isnan(x)]
####################################################################################################################
    # Q_in = pd.read_csv('./Model_Data_Filled_20082023/LO_Inflows_20082023.csv')
    NOx_In = pd.read_csv('./Model_Data_Filled_20082023/Daily_NOx_External_Loads_20082023.csv') #mg
    S65E_NO_data = pd.read_csv('./Model_Data_Filled_20082023/S65E_NO_Interpolated.csv') #mg/m3
    Load_ext = pd.read_csv('./Data/%s/ts_data/LO_External_Loadings_%s_3M.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    Q_in = pd.read_csv('./Data/%s/ts_data/LO_Inflows_BK_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))

########################################################################################################################
    M_var.Lake_Stage[0] = Pre_defined_Variables.begstageCS
    M_var.Lake_Stage[1] = Pre_defined_Variables.begstageCS
    M_var.DecTree_Relslevel[0] = np.nan
    M_var.DecTree_Relslevel[1] = np.nan
    if startdate.month == LO_Model['date'].iloc[2].month and startdate.day == LO_Model['date'].iloc[2].day:
        X1 = 'SimDay1'
    elif begdateCS.year == LO_Model['date'].iloc[2].year and begdateCS.month == LO_Model['date'].iloc[2].month and begdateCS.day == LO_Model['date'].iloc[2].day:
        X1 = 'CS start date'
    else:
        X1 = LO_Model['date'].iloc[2]
    M_var.DayFlags[2] = X1     
    StartStorage = Stg_Sto_Ar.stg2sto(Pre_defined_Variables.startstage,0)
    M_var.Storage[0] = StartStorage
    M_var.Storage[1] = StartStorage
    # Flood = np.zeros(n_rows, dtype = object)
    ##Here, I will insert the Storage Deviaiton Values as Input!
    Storage_dev = Data.Stroage_dev_df['DS_dev'] 
    #Create a Choose Function for AP Post Baseflow
    # if Pre_defined_Variables.Opt_AdapProt == 0:
    #     C = 450
    # elif Pre_defined_Variables.Opt_AdapProt == 1:
    #     C = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    # Choose_1 = C
    Choose_1 = 450 #cfs
##############################################################################################################    
    #Read Required Data
    Q_in = pd.read_csv('./Model_Data_Filled_20082023/LO_Inflows_20082023.csv')
    # Q_out = pd.read_csv('./Model_Data_Filled_20082023/LO_Outflows_20082023.csv')
    Temp_data = pd.read_csv('./Model_Data_Filled_20082023/Temp_Avg_20082023.csv')
    DO_data = pd.read_csv('./Model_Data_Filled_20082023/LO_Avg_DO_20082023.csv')
    RAD_data = pd.read_csv('./Model_Data_Filled_20082023/LO_RADT_20082023.csv')
    Storage = pd.read_csv('./Model_Data_Filled_20082023/Average_LO_Storage_20082023.csv')
    Storage_deviation = pd.read_csv('./Model_Data_Filled_20082023/Storage_Dev_20082023.csv')
    Chla_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_Merged_Chla.csv') #microgram/L
    Chla_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_Merged_Chla.csv') #microgram/L
    NOx_In = pd.read_csv('./Model_Data_Filled_20082023/Daily_NOx_External_Loads_20082023.csv') #mg
    S65E_NO_data = pd.read_csv('./Model_Data_Filled_20082023/S65E_NO_Interpolated.csv') #mg/m3
    Chla_In = pd.read_csv('./Model_Data_Filled_20082023/Chla_Loads_In_20082023.csv') #mg
    S65E_Chla_data = pd.read_csv('./Model_Data_Filled_20082023/S65E_Chla_20082023.csv') #mg/m3

    Photoperiod = pd.read_csv('./Model_Data_Filled_20082023/PhotoPeriod_20082023.csv') 
    Photoperiod['Date'] = pd.to_datetime(Photoperiod['Date'])
    LO_DIP_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_OP.csv') #mg/m3
    LO_DIP_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_OP.csv') #mg/m3
    LO_DIN_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_DIN.csv') #mg/m3
    LO_DIN_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_DIN.csv') #mg/m3

    RAD = RAD_data['RAD'].astype(float) * 4.6*1000 

    Temp = Temp_data['Mean_T'].astype(float)
    DO = DO_data['DO'].astype(float)
    External_NO = NOx_In['NOx_Load_mg'].astype(float) #mg
    S65E_NO = S65E_NO_data['NO'].astype(float) #mg/m3
    
    External_Chla = Chla_In['Load_mg'].astype(float)*3 #mg
    S65E_Chla = S65E_Chla_data['Chla'].astype(float)

    # N-S Procedure
    N_Per = 0.43
    S_Per = 0.57
    
    volume  = Storage['Storage_m3'].astype(float) #m3
    volume_N = volume*N_Per
    volume_S = volume*S_Per
    Storage_dev = Storage_deviation['DS_dev'].astype(float) #acft
    Stage = Storage['Stage_ft'].astype(float) * 0.3048 #m
    Q_I = Q_in['Inflows_cmd'].astype(float) #m3
    # Q_O = Q_out['LO_Outflows_cmd'].astype(float) #m3

   
    # Simulated Q
    S77_Q = LOONE_Q_Outputs['S77_Q'] * 0.0283168 * 86400 #cfs to cubic meters per day
    S308_Q = LOONE_Q_Outputs['S308_Q'] * 0.0283168 * 86400 #cfs to cubic meters per day
    TotRegSo = LOONE_Q_Outputs['TotRegSo'] #acft/day
    TotRegEW = LOONE_Q_Outputs['TotRegEW'] #acft/day
   
    # # Observed S77 S308 South
    # Outflows_Obs = pd.read_csv('./Data/%s/ts_data/Flow_df_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    # S77_Q = Outflows_Obs['S77_Out']
    # S308_Q = Outflows_Obs['S308_Out']
    # TotRegSo = Outflows_Obs[['S351_Out','S354_Out','S352_Out','L8_Out']].sum(axis=1)/1233.48    #m3/day to acft
    
    Q_O = S77_Q + S308_Q + TotRegSo * 1233.48  #cmd

    NH4_N_Obs = LO_DIN_N_data['NH4'].astype(float) #mg/m3
    NH4_S_Obs = LO_DIN_S_data['NH4'].astype(float) #mg/m3
    NH4_Obs = (NH4_N_Obs+NH4_S_Obs)/2
    NOx_N_Obs = LO_DIN_N_data['NO'].astype(float) #mg/m3
    NOx_S_Obs = LO_DIN_S_data['NO'].astype(float) #mg/m3
    NOx_Obs = (NOx_N_Obs+NOx_S_Obs)/2
    Chla_N = Chla_N_data['Chla'].astype(float) #microgram/L = mg/m3
    Chla_S = Chla_S_data['Chla'].astype(float) #microgram/L = mg/m3
    Chla = (Chla_N+Chla_S)/2 
    DIN_N_Obs = LO_DIN_N_data['DIN'].astype(float)
    DIN_S_Obs = LO_DIN_S_data['DIN'].astype(float)
    DIN_Obs = (DIN_N_Obs+DIN_S_Obs)/2
    DIP_N_Obs = LO_DIP_N_data['OP'].astype(float)
    DIP_S_Obs = LO_DIP_S_data['OP'].astype(float)
    DIP_Obs = (DIP_N_Obs+DIP_S_Obs)/2
    
    Nit_Denit_Opt = 'Opt1'
    
    #NO Atmospheric Deposition
    Atm_Deposition = (714.6*1000*1000) #mg/day
    Atm_Deposition_N = N_Per * Atm_Deposition
    Atm_Deposition_S = S_Per * Atm_Deposition
    
    #Calibration parameters 
    Cal_Par_Opt = pd.read_csv('./Model_Data_Filled_20082023/Cal_Par.csv')
    Cal_Par_NO = Cal_Par_Opt['Par_NO']
    Cal_Par_Chla = Cal_Par_Opt['Par_Chla']
 
       
    #nitrification
    #Either a simple method where Nit_R_T = nitrification rate in water column (Tot_nit_R) * Temp Coeff for nitrification * (T-20)
    #Or a more detailed method R = k0 + k * fam * fox
    # two options for fam and fox 
    #Option1
    # fam = (Cam/Ks+Cam) and fox = (Cox/Ksox+Cox)
    #Option 2
    #Fam = Cam and fox = (1-foxmin) * (Cox-Coxc/Coxo-Coxc) ^10a + foxmin
    #The nitrification rate for the simple method
    Tot_nit_R = 0.3
    # nitrification parameters for the detailed method
    kni0_NO = Cal_Par_NO[0]
    kni0_Chla = Cal_Par_Chla[0]
    kni_NO = Cal_Par_NO[1]
    kni_Chla = Cal_Par_Chla[1]

    Theta_ni = 1.06
    
    #Denitrification
    #Either a simple method where Denit_R_T = denitrification rate in water column (Tot_denit_R) * Temp Coeff for denitrification * (T-20)
    #Or a more detailed method R = k0 + k * fni * fox
    # two options for fam and fox 
    #Option1
    # fni = (Cni/Ks+Cni) and fox = 1- (Cox/Ksox+Cox)
    #Option 2
    #Fni = Cni and fox = (Coxc-Cox/Coxc-Coxo) 
    #The denitrification rate for the simple method
    Tot_denit_R = 0.1
    #The denitrification rate for the simple method
    kden0_NO =  Cal_Par_NO[2]
    kden0_Chla =  Cal_Par_Chla[2]

    kden_NO =  Cal_Par_NO[3]
    kden_Chla =  Cal_Par_Chla[3]

    Theta_deni = 1.06
    
    #Light Attenuation due to water Kw and algae Kc
    Kw_NO = Cal_Par_NO[4] #1.7
    Kw_Chla = Cal_Par_Chla[4] #1.7
    Kc_NO = Cal_Par_NO[5] #0.029
    Kc_Chla = Cal_Par_Chla[5] #0.029
    #Light limitation and inhibition coefficients K1 and K2
    K1_NO = Cal_Par_NO[6] #100
    K1_Chla = Cal_Par_Chla[6] #100
    K2_NO = Cal_Par_NO[7] #700
    K2_Chla = Cal_Par_Chla[7] #700
    # Maximum Phytoplankton growth rate
    G_max_NO = Cal_Par_NO[8] #1.5
    G_max_Chla = Cal_Par_Chla[8] #1.5

    #Half Saturation Coefficients
    K_NH_NO = Cal_Par_NO[9] #0.1
    K_NH_Chla = Cal_Par_Chla[9] #0.1
    K_Nitr_NO = Cal_Par_NO[10]
    K_Nitr_Chla = Cal_Par_Chla[10]
    K_DIN_NO = K_NH_NO + K_Nitr_NO
    K_DIN_Chla = K_NH_Chla + K_Nitr_Chla
    KP_NO = Cal_Par_NO[11]
    KP_Chla = Cal_Par_Chla[11]
    K_TN_NO = Cal_Par_NO[12] #0.1
    K_TN_Chla = Cal_Par_Chla[12] #0.1
    YNOChla_NO = Cal_Par_NO[13] #0.1
    YNOChla_Chla = Cal_Par_Chla[13] #0.1

    
    #Temperatures
    T_opt_NO = Cal_Par_NO[14] #15 #C
    T_min_NO = Cal_Par_NO[15] #10 #C
    T_max_NO = Cal_Par_NO[16] #30 #C
    
    T_opt_Chla = Cal_Par_Chla[14] #15 #C
    T_min_Chla = Cal_Par_Chla[15] #10 #C
    T_max_Chla = Cal_Par_Chla[16] #30 #C

    # Dissolved Oxygen    
    KDO_NO = Cal_Par_NO[17]
    KDO_Chla = Cal_Par_Chla[17]

    #Sediment Release
    S_NO_NO = Cal_Par_NO[18] #mg/m2/d
    S_NO_Chla = Cal_Par_Chla[18] #mg/m2/d

    Theta_NO = 1.06
    NO_Temp_Adj = Theta_NO**(Temp-20)  
    Area = 1730 *1E6 #m2
    Area_N = Area * N_Per
    Area_S = Area * S_Per
    
    X = len(Q_in.index)
    
    #Nitrogen and Phosphorus Limiting
    f_P_NO = DIP_Obs/(KP_NO+DIP_Obs)
    f_P_Chla = DIP_Obs/(KP_Chla+DIP_Obs)

    f_N_NO = DIN_Obs/(K_DIN_NO+DIN_Obs)
    f_N_Chla = DIN_Obs/(K_DIN_Chla+DIN_Obs)

    f_P_N_NO = DIP_N_Obs/(KP_NO+DIP_N_Obs)
    f_P_N_Chla = DIP_N_Obs/(KP_Chla+DIP_N_Obs)

    f_P_S_NO = DIP_S_Obs/(KP_NO+DIP_S_Obs)
    f_P_S_Chla = DIP_S_Obs/(KP_Chla+DIP_S_Obs)

    f_N_N_NO = DIN_N_Obs/(K_DIN_NO+DIN_N_Obs)
    f_N_N_Chla = DIN_N_Obs/(K_DIN_Chla+DIN_N_Obs)

    f_N_S_NO = DIN_S_Obs/(K_DIN_NO+DIN_S_Obs)
    f_N_S_Chla = DIN_S_Obs/(K_DIN_Chla+DIN_S_Obs)

    Q_I_M = np.zeros(X,dtype = object)
    Q_O_M = np.zeros(X,dtype = object)
    Q_N2S = np.zeros(X,dtype = object)
    External_NO_M = np.zeros(X,dtype = object)
    External_Chla_M = np.zeros(X,dtype = object)

    Nit_R_N_NO = np.zeros(X,dtype = object)
    Nit_R_N_T_NO = np.zeros(X,dtype = object)
    Nit_R_S_NO = np.zeros(X,dtype = object)
    Nit_R_S_T_NO = np.zeros(X,dtype = object)
    Denit_R_N_NO = np.zeros(X,dtype = object)
    Denit_R_N_T_NO = np.zeros(X,dtype = object)
    Denit_R_S_NO = np.zeros(X,dtype = object)
    Denit_R_S_T_NO = np.zeros(X,dtype = object)
    

    fT_NO = np.zeros(X,dtype = object)
    fL_NO = np.zeros(X,dtype = object)
    
    fT_Chla = np.zeros(X,dtype = object)
    fL_Chla = np.zeros(X,dtype = object)

    NO_N = np.zeros(X,dtype = object)
    NO_S = np.zeros(X,dtype = object)
    NO_MEAN = np.zeros(X,dtype = object)
    
    NO_Load_Cal = np.zeros(X,dtype = object)
    NO_Load_StL = np.zeros(X,dtype = object)
    NO_Load_South = np.zeros(X,dtype = object)

    Sim_Chla_N = np.zeros(X,dtype = object)
    Sim_Chla_S = np.zeros(X,dtype = object)
    Sim_Chla = np.zeros(X,dtype = object)
    
    Chla_Load_Cal = np.zeros(X,dtype = object)
    Chla_Load_StL = np.zeros(X,dtype = object)
    Chla_Load_South = np.zeros(X,dtype = object)

    
    #Water depth could be stage - B.L. (variable) or a constant average depth (e.g., 2.5 m)

    z = np.zeros(X,dtype = object)
    # B_L = 1 #m
    # z = Stage - B_L #m

    for i in range(X):
        z[i] = 2.5
    
    fox_min = Cal_Par_NO[19]
    DO_Cr_n = Cal_Par_NO[20]
    DO_Opt_n = Cal_Par_NO[21]
    a = Cal_Par_NO[22]
    DO_Cr_d = Cal_Par_NO[23]
    DO_Opt_d = Cal_Par_NO[24]
   
    vc = Cal_Par_Chla[28] # phytoplankton settling velocity (m/d).

    K_r = Cal_Par_Chla[25]
    Theta_r = 1.06
    K_r_T = K_r*Theta_r**(Temp-20) 
    # K_r_T = K_r

    YNHChla = Cal_Par_Chla[26] #0.05
    #Chla mortality
    k_m = Cal_Par_Chla[27]
    Theta_m = 1.06
    K_m_T = k_m*Theta_m**(Temp-20) 
    # K_m_T = k_m
    # Grazing = 0.01
    # Theta_G = 1.06
    # Grazing_T = Grazing*Theta_G**(Temp-20)
    
    Im = Cal_Par_Chla[31]

    
    Grazing = 0
    Theta_G = 1.06
    Grazing_T = Grazing*Theta_G**(Temp-20)
    NO_N[0] = NOx_N_Obs[0]
    NO_S[0] = NOx_S_Obs[0]
    
    NO_MEAN[0] = (NO_N[0]+NO_S[0])*0.5
    
    Sim_Chla_N[0] = Chla_N[0]
    Sim_Chla_S[0] = Chla_S[0]
    Sim_Chla[0] = Chla[0]

    # Map the monthly values to the DataFrame
    Nitro_Model_Output = pd.DataFrame(Q_in['date'],columns=['date'])
    Nitro_Model_Output['date'] = pd.to_datetime(Nitro_Model_Output['date'])
    print("LOONE Nitrogen Module is Running!")
    S65_P_Conc = pd.read_csv('./Data/%s/ts_data/S65E_TP_Conc_3M.csv'%Pre_defined_Variables.Schedule)
    S65_P = S65_P_Conc['TP_mg/m3'].astype(float)
    L_ext = Load_ext['TP_Loads_In_mg'] #mg
    Atm_Dep_N = TP_Variables.N_Per * Load_ext['Atm_Loading_mg']
    Atm_Dep_S = TP_Variables.S_Per * Load_ext['Atm_Loading_mg']
    Wind_ShearStr = pd.read_csv('./Data/%s/ts_data/WindShearStress_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    W_SS = Wind_ShearStr['ShearStress'] #Dyne/cm2
    Current_Stress =  pd.read_csv('./Data/%s/ts_data/Current_ShearStress.csv'%Pre_defined_Variables.Schedule)
    Current_SS = Current_Stress['Current_Stress']
    Bottom_Stress = W_SS + Current_SS
    nu_ts = pd.read_csv('./Data/%s/ts_data/nu_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    LO_BL = 0.5 # m (Bed Elevation of LO)
    # LO_WD = pd.to_numeric(Stage_Storage['Stage_m'])-LO_BL
    g = 9.8 #m/s2 gravitational acceleration
    Cal_Res = pd.read_csv('./Data/%s/nondominated_Sol_var.csv'%Pre_defined_Variables.Schedule)
    Par = Cal_Res['Par']
    d_c = Par[20] # m (particle diameter 10 microm /1E6 to convert to m) clay
    d_s = Par[21] # m sand
    nu_d = nu_ts['nu']
    R = 1.65 #submerged specific gravity (1.65 for quartz in water)
    C_1_c = Par[16]
    C_2_c = Par[17]
    C_1_s = Par[18]
    C_2_s = Par[19]
    #Parameters associated with sediment resuspension
    Crtcl_ShStr = Par[22] #0.32 #Dyne/cm2
    E_0 = Par[23]
    E_1 = Par[24]
    E_2 = Par[25]
    Td = Par[26] #days
    Stage2ar = np.zeros(n_rows,dtype = object)
    LO_WD = np.zeros(n_rows,dtype = object)
    Lake_O_Storage_N = np.zeros(n_rows,dtype = object)
    Lake_O_Storage_S = np.zeros(n_rows,dtype = object)
    Lake_O_A_N = np.zeros(n_rows,dtype = object)
    Lake_O_A_S = np.zeros(n_rows,dtype = object)
    Lake_O_A_M_N = np.zeros(n_rows,dtype = object)
    Lake_O_A_S_N = np.zeros(n_rows,dtype = object)
    Lake_O_A_R_N = np.zeros(n_rows,dtype = object)
    Lake_O_A_P_N = np.zeros(n_rows,dtype = object)
    Lake_O_A_M_S = np.zeros(n_rows,dtype = object)
    Lake_O_A_S_S = np.zeros(n_rows,dtype = object)
    Lake_O_A_R_S = np.zeros(n_rows,dtype = object)
    Lake_O_A_P_S = np.zeros(n_rows,dtype = object)
    DIP_Lake_N = np.zeros(n_rows,dtype = object)
    DIP_Lake_S = np.zeros(n_rows,dtype = object)
    TP_Lake_Mean = np.zeros(n_rows,dtype = object)
    J_des_M_N = np.zeros(n_rows,dtype = object)
    J_des_S_N = np.zeros(n_rows,dtype = object)
    J_des_R_N = np.zeros(n_rows,dtype = object)
    J_des_P_N = np.zeros(n_rows,dtype = object)
    J_des_M_S = np.zeros(n_rows,dtype = object)
    J_des_S_S = np.zeros(n_rows,dtype = object)
    J_des_R_S = np.zeros(n_rows,dtype = object)
    J_des_P_S = np.zeros(n_rows,dtype = object)
    J_ads_M_N = np.zeros(n_rows,dtype = object)
    J_ads_S_N = np.zeros(n_rows,dtype = object)
    J_ads_R_N = np.zeros(n_rows,dtype = object)
    J_ads_P_N = np.zeros(n_rows,dtype = object)
    J_ads_M_S = np.zeros(n_rows,dtype = object)
    J_ads_S_S = np.zeros(n_rows,dtype = object)
    J_ads_R_S = np.zeros(n_rows,dtype = object)
    J_ads_P_S = np.zeros(n_rows,dtype = object)
    P_sed_M_N = np.zeros(n_rows,dtype = object)
    P_sed_S_N = np.zeros(n_rows,dtype = object)
    P_sed_R_N = np.zeros(n_rows,dtype = object)
    P_sed_P_N = np.zeros(n_rows,dtype = object)
    P_sed_M_S = np.zeros(n_rows,dtype = object)
    P_sed_S_S = np.zeros(n_rows,dtype = object)
    P_sed_R_S = np.zeros(n_rows,dtype = object)
    P_sed_P_S = np.zeros(n_rows,dtype = object)
    J_sedburial_M_N = np.zeros(n_rows,dtype = object)
    J_sedburial_S_N = np.zeros(n_rows,dtype = object)
    J_sedburial_R_N = np.zeros(n_rows,dtype = object)
    J_sedburial_P_N = np.zeros(n_rows,dtype = object)
    J_sedburial_M_S = np.zeros(n_rows,dtype = object)
    J_sedburial_S_S = np.zeros(n_rows,dtype = object)
    J_sedburial_R_S = np.zeros(n_rows,dtype = object)
    J_sedburial_P_S = np.zeros(n_rows,dtype = object)
    J_Γburial_M_N = np.zeros(n_rows,dtype = object)
    J_Γburial_S_N = np.zeros(n_rows,dtype = object)
    J_Γburial_R_N = np.zeros(n_rows,dtype = object)
    J_Γburial_P_N = np.zeros(n_rows,dtype = object)
    J_Γburial_M_S = np.zeros(n_rows,dtype = object)
    J_Γburial_S_S = np.zeros(n_rows,dtype = object)
    J_Γburial_R_S = np.zeros(n_rows,dtype = object)
    J_Γburial_P_S = np.zeros(n_rows,dtype = object)
    Γ_M_N = np.zeros(n_rows,dtype = object)
    Γ_S_N = np.zeros(n_rows,dtype = object)
    Γ_R_N = np.zeros(n_rows,dtype = object)
    Γ_P_N= np.zeros(n_rows,dtype = object)
    Γ_M_S = np.zeros(n_rows,dtype = object)
    Γ_S_S = np.zeros(n_rows,dtype = object)
    Γ_R_S = np.zeros(n_rows,dtype = object)
    Γ_P_S = np.zeros(n_rows,dtype = object)
    DIP_pore_M_N = np.zeros(n_rows,dtype = object)
    DIP_pore_S_N = np.zeros(n_rows,dtype = object)
    DIP_pore_R_N = np.zeros(n_rows,dtype = object)
    DIP_pore_P_N = np.zeros(n_rows,dtype = object)
    DIP_pore_M_S = np.zeros(n_rows,dtype = object)
    DIP_pore_S_S = np.zeros(n_rows,dtype = object)
    DIP_pore_R_S = np.zeros(n_rows,dtype = object)
    DIP_pore_P_S = np.zeros(n_rows,dtype = object)
    TP_Lake_N = np.zeros(n_rows,dtype = object)
    TP_Lake_S = np.zeros(n_rows,dtype = object)
    Sed_Resusp_M_N = np.zeros(n_rows,dtype = object)
    Sed_Resusp_S_N = np.zeros(n_rows,dtype = object)
    Sed_Resusp_R_N = np.zeros(n_rows,dtype = object)
    Sed_Resusp_P_N = np.zeros(n_rows,dtype = object)
    Sed_Resusp_M_S = np.zeros(n_rows,dtype = object)
    Sed_Resusp_S_S= np.zeros(n_rows,dtype = object)
    Sed_Resusp_R_S = np.zeros(n_rows,dtype = object)
    Sed_Resusp_P_S = np.zeros(n_rows,dtype = object)
    
    Suspended_Sed_N = np.zeros(n_rows,dtype = object)
    Suspended_Sed_S = np.zeros(n_rows,dtype = object)

    J_decomp_M_N = np.zeros(n_rows,dtype = object)
    J_decomp_S_N = np.zeros(n_rows,dtype = object)
    J_decomp_R_N = np.zeros(n_rows,dtype = object)
    J_decomp_P_N = np.zeros(n_rows,dtype = object)
    J_decomp_M_S = np.zeros(n_rows,dtype = object)
    J_decomp_S_S = np.zeros(n_rows,dtype = object)
    J_decomp_R_S = np.zeros(n_rows,dtype = object)
    J_decomp_P_S = np.zeros(n_rows,dtype = object)
     
    Settling_P_N = np.zeros(n_rows,dtype = object)
    Settling_P_S = np.zeros(n_rows,dtype = object)
    
    P_diff_M_N = np.zeros(n_rows,dtype = object)
    P_diff_S_N = np.zeros(n_rows,dtype = object)
    P_diff_R_N = np.zeros(n_rows,dtype = object)
    P_diff_P_N = np.zeros(n_rows,dtype = object)
    P_diff_M_S = np.zeros(n_rows,dtype = object)
    P_diff_S_S = np.zeros(n_rows,dtype = object)
    P_diff_R_S = np.zeros(n_rows,dtype = object)
    P_diff_P_S = np.zeros(n_rows,dtype = object)
    
    L_ext_M = np.zeros(n_rows,dtype = object)

    Q_I = Q_in['Flow_cmd'].astype(float) 
    Q_I_M = np.zeros(n_rows,dtype = object)
    Q_O = np.zeros(n_rows,dtype = object)
    Q_O_M = np.zeros(n_rows,dtype = object)
    Q_N2S = np.zeros(n_rows,dtype = object)
    Storage = np.zeros(n_rows,dtype = object)
    P_Load_Cal = np.zeros(n_rows,dtype = object)
    P_Load_StL = np.zeros(n_rows,dtype = object)
    P_Load_South = np.zeros(n_rows,dtype = object)
    
    v_settle_N_c = np.zeros(n_rows,dtype = object)
    v_settle_N_s = np.zeros(n_rows,dtype = object)
    v_settle_N = np.zeros(n_rows,dtype = object)
    v_settle_S_c = np.zeros(n_rows,dtype = object)
    v_settle_S_s = np.zeros(n_rows,dtype = object)
    v_settle_S = np.zeros(n_rows,dtype = object)
    Ext_Load_Rate_N = np.zeros(n_rows,dtype = object)
    Ext_Load_Rate_S = np.zeros(n_rows,dtype = object)
    Load_Out_N = np.zeros(n_rows,dtype = object)
    Load_Out_S = np.zeros(n_rows,dtype = object)
    ##Initial Values##
    #S.A. is calculated based on the Lake's previous time step Stage, but for the S.A. at i=0 I used same time step Stage!
    StartStorage = Stg_Sto_Ar.stg2sto(Pre_defined_Variables.startstage,0)
    Stage2ar[0] = Stg_Sto_Ar.stg2ar(M_var.Lake_Stage[0],0)
    Stage2ar[1] = Stg_Sto_Ar.stg2ar(M_var.Lake_Stage[1],0)
    Storage[0] = StartStorage #ac-ft
    Storage[1] = Stg_Sto_Ar.stg2sto(M_var.Lake_Stage[1],0) #ac-ft
    #TP_MassBalanceModel Initial Values.
    TP_Lake_N[0] = 235 #mg/m3
    TP_Lake_S[0] = 255 #mg/m3
    TP_Lake_Mean[0] = (TP_Lake_N[0] + TP_Lake_S[0])/2
    Γ_M_N[0] = 25 #mg/kg 
    Γ_S_N[0] = 25 #mg/kg 
    Γ_R_N[0] = 25 #mg/kg 
    Γ_P_N[0] = 25 #mg/kg 
    Γ_M_S[0] = 25 #mg/kg
    Γ_S_S[0] = 25 #mg/kg 
    Γ_R_S[0] = 25 #mg/kg 
    Γ_P_S[0] = 25 #mg/kg  
    
    DIP_pore_M_N[0] = 700#760 #mg/m3 
    DIP_pore_S_N[0] = 240#205 #mg/m3 
    DIP_pore_R_N[0] = 240#205 #mg/m3 
    DIP_pore_P_N[0] = 160#160 #mg/m3 
    DIP_pore_M_S[0] = 700#760 #mg/m3 
    DIP_pore_S_S[0] = 240#205 #mg/m3 
    DIP_pore_R_S[0] = 240#205 #mg/m3
    DIP_pore_P_S[0] = 160#160 #mg/m3 
    P_sed_M_N[0] = 1100 #mg/kg 
    P_sed_S_N[0] = 300 #mg/kg  
    P_sed_R_N[0] = 300 #mg/kg 
    P_sed_P_N[0] = 200 #mg/kg 
    P_sed_M_S[0] = 1100 #mg/kg 
    P_sed_S_S[0] = 300 #mg/kg 
    P_sed_R_S[0] = 300 #mg/kg 
    P_sed_P_S[0] = 200 #mg/kg 
    
    
    Θ_M = 1-((TP_Variables.Bulk_density_M/TP_Variables.Particle_density_M)*((100-TP_Variables.Per_H2O_M)/100))
    Θ_S = 1-((TP_Variables.Bulk_density_S/TP_Variables.Particle_density_S)*((100-TP_Variables.Per_H2O_S)/100))
    Θ_R = 1-((TP_Variables.Bulk_density_R/TP_Variables.Particle_density_R)*((100-TP_Variables.Per_H2O_R)/100))
    Θ_P = 1-((TP_Variables.Bulk_density_P/TP_Variables.Particle_density_P)*((100-TP_Variables.Per_H2O_P)/100))
    #Mass of sediment in surfacial mix Mud layer in the North Region(kg)
    Mass_sed_M_N = TP_Variables.A_Mud_N * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_M)/100) * TP_Variables.Bulk_density_M * 1000
    #Mass of sediment in surfacial mix Sand layer in the North Region(kg)
    Mass_sed_S_N = TP_Variables.A_Sand_N * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_S)/100) * TP_Variables.Bulk_density_S * 1000
    #Mass of sediment in surfacial mix Rock layer in the North Region(kg)
    Mass_sed_R_N = TP_Variables.A_Rock_N * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_R)/100) * TP_Variables.Bulk_density_R * 1000
    #Mass of sediment in surfacial mix Peat layer in the North Region(kg)
    Mass_sed_P_N = TP_Variables.A_Peat_N * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_P)/100) * TP_Variables.Bulk_density_P * 1000
    #Mass of sediment in surfacial mix Mud layer in the South Region(kg)
    Mass_sed_M_S = TP_Variables.A_Mud_S * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_M)/100) * TP_Variables.Bulk_density_M * 1000
    #Mass of sediment in surfacial mix Sand layer in the South Region(kg)
    Mass_sed_S_S = TP_Variables.A_Sand_S * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_S)/100) * TP_Variables.Bulk_density_S * 1000
    #Mass of sediment in surfacial mix Rock layer in the South Region(kg)
    Mass_sed_R_S = TP_Variables.A_Rock_S * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_R)/100) * TP_Variables.Bulk_density_R * 1000
    #Mass of sediment in surfacial mix Peat layer in the South Region(kg)
    Mass_sed_P_S = TP_Variables.A_Peat_S * TP_Variables.Z_sed * ((100-TP_Variables.Per_H2O_P)/100) * TP_Variables.Bulk_density_P * 1000

    #####################################################################################################################
    M_var.Zone_Code[0] = LO_FNs.Zone_Code(M_var.Lake_Stage[0],df_WSMs['A'].iloc[0],df_WSMs['B'].iloc[0],df_WSMs['C'].iloc[0],df_WSMs['D3'].iloc[0],df_WSMs['D2'].iloc[0],df_WSMs['D1'].iloc[0],df_WSMs['D0'].iloc[0],df_WSMs['WSM1'].iloc[0])
    M_var.LO_Zone[0] = LO_FNs.LO_Zone(M_var.Zone_Code[0])
    LO_EcoRecEnv = pd.read_csv('./Data/LO_EcoRecEnv_Interpolated_ts.csv')

    for i in range(n_rows-2):
        M_var.WSM_Zone[i+2] = LO_FNs.WSM_Zone(M_var.Lake_Stage[i+1],df_WSMs.at[i+1, 'WSM4'],df_WSMs.at[i+1, 'WSM3'],df_WSMs.at[i+1, 'WSM2'],df_WSMs.at[i+1, 'WSM1'])
    #Calculate Daily Maximum Water Supply
    # Note that in LOSA_dmd we used (i) because this file starts from 1/1/2008 so i at this point =0.
    #Cutbacks are determined based on the WSM Zone. 
        M_var.Max_Supply[i+2] = LO_FNs.Max_Supply(M_var.WSM_Zone[i+2],Water_dmd.at[i, 'Daily_demand'],Pre_defined_Variables.Z1_cutback,Pre_defined_Variables.Z2_cutback,Pre_defined_Variables.Z3_cutback,Pre_defined_Variables.Z4_cutback)
    #Actual Daily Water supply
        M_var.LOSA_Supply[i+2] = LO_FNs.LOSA_Supply(M_var.WSM_Zone[i+2],LO_Model.at[i+2, 'LOSA_dmd_SFWMM'],M_var.Max_Supply[i+2],Pre_defined_Variables.Opt_LOSAws)
    # NetInflow - LOSA Supply
        M_var.NI_Supply[i+2] = LO_Model.at[i+2, 'Net_Inflow'] - M_var.LOSA_Supply[i+2]
    #TODO Note: for the pass statement, We will read the Daily Water supply from the SFWMM as an input.
    #Calculate the cutback where Cutback = Demand - Supply
        ctbk = LO_Model.at[i+2, 'LOSA_dmd_SFWMM'] - M_var.LOSA_Supply[i+2]        
        M_var.Cut_back[i+2] = ctbk
    #Calculate percentage of the demand that is not supplied for each day
        if LO_Model.at[i+2, 'LOSA_dmd_SFWMM'] == 0:
            DNS = 0
        else:
            DNS = (M_var.Cut_back[i+2] / LO_Model.at[i+2, 'LOSA_dmd_SFWMM'])*100
        M_var.Dem_N_Sup[i+2] = DNS
    # Calculate the Zone Code
    #Note that to calculate the Zone Code in Dec 31 2020 we needed the WSM and breakpoint zones in 1/1/2021!
    #Note Also that i = 0 in Stage indicates Dec 30 1964 while i = 0 in df_WSMs indicates Dec 31 1964!
        M_var.Zone_Code[i+1] = LO_FNs.Zone_Code(M_var.Lake_Stage[i+1],df_WSMs.at[i+1, 'A'],df_WSMs.at[i+1, 'B'],df_WSMs.at[i+1, 'C'],df_WSMs.at[i+1, 'D3'],df_WSMs.at[i+1, 'D2'],df_WSMs.at[i+1, 'D1'],df_WSMs.at[i+1, 'D0'],df_WSMs.at[i+1, 'WSM1'])
    #Generate the Zone Column based on the corresponding Zone Code.
        M_var.LO_Zone[i+1] = LO_FNs.LO_Zone(M_var.Zone_Code[i+1])
        M_var.Zone_D_Trib[i] = Dec_Tree_FNs.Zone_D_Trib(TC_LONINO_df.at[i, 'Tributary_Condition'],Pre_defined_Variables.Opt_NewTree)
        M_var.Zone_D_stage[i] = Dec_Tree_FNs.Zone_D_stage(M_var.Lake_Stage[i+1],df_WSMs.at[i, 'C-b'])
        M_var.Zone_D_Seas[i] = Dec_Tree_FNs.Zone_D_Seas(TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],M_var.Zone_D_Trib[i],Pre_defined_Variables.Opt_NewTree)
        M_var.Zone_D_MSeas[i] = Dec_Tree_FNs.Zone_D_MSeas(TC_LONINO_df.at[i, 'LONINO_MultiSeasonal_Classes'])
        M_var.Zone_D_Branch_Code[i] = M_var.Zone_D_Trib[i]*1000 + M_var.Zone_D_stage[i]*100 + M_var.Zone_D_Seas[i]*10 + M_var.Zone_D_MSeas[i]*1
        M_var.Zone_D_Rel_Code[i] = Dec_Tree_FNs.Zone_D_Rel_Code(M_var.Zone_D_Branch_Code[i],Pre_defined_Variables.Opt_DecTree)  
        M_var.Zone_C_Trib[i] = Dec_Tree_FNs.Zone_C_Trib(TC_LONINO_df.at[i, 'Tributary_Condition'],Pre_defined_Variables.Opt_NewTree)
        M_var.Zone_C_Seas[i] = Dec_Tree_FNs.Zone_C_Seas(TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],Pre_defined_Variables.Opt_NewTree)  
        M_var.Zone_C_MSeas[i] = Dec_Tree_FNs.Zone_C_MSeas(TC_LONINO_df.at[i, 'LONINO_MultiSeasonal_Classes'])
        M_var.Zone_C_MetFcast[i] = Dec_Tree_FNs.Zone_C_MetFcast(M_var.Zone_C_Seas[i],TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],Pre_defined_Variables.Zone_C_MetFcast_Indicator)
        M_var.Zone_C_Branch_Code[i] = M_var.Zone_C_Trib[i]*1000 + M_var.Zone_C_MetFcast[i]*100 + M_var.Zone_C_Seas[i]*10 + M_var.Zone_C_MSeas[i]*1
        M_var.Zone_C_Rel_Code[i] = Dec_Tree_FNs.Zone_C_Rel_Code(M_var.Zone_C_Branch_Code[i],Pre_defined_Variables.Opt_DecTree)   
        M_var.Zone_B_Trib[i] = Dec_Tree_FNs.Zone_B_Trib(TC_LONINO_df.at[i, 'Tributary_Condition'],Pre_defined_Variables.Opt_NewTree)
        M_var.Zone_B_Stage[i] = Dec_Tree_FNs.Zone_B_Stage(M_var.Lake_Stage[i+1],Seasons.at[i, 'Season'])
        M_var.Zone_B_Seas[i] = Dec_Tree_FNs.Zone_B_Seas(TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'])
        M_var.Zone_B_Branch_Code[i] = M_var.Zone_B_Trib[i]*1000 + M_var.Zone_B_Stage[i]*100 + DecTree_df.at[i, 'Zone_B_MetFcast']*10 + M_var.Zone_B_Seas[i]*1
        M_var.Zone_B_Rel_Code[i] = Dec_Tree_FNs.Zone_B_Rel_Code(M_var.Zone_B_Branch_Code[i],Pre_defined_Variables.Opt_DecTree)
        M_var.DecTree_Relslevel[i+2] = LO_FNs.DecTree_Relslevel(M_var.Zone_Code[i+1],M_var.Zone_D_Rel_Code[i],M_var.Zone_C_Rel_Code[i],M_var.Zone_B_Rel_Code[i])
        if i >= 3:
            if startdate.month == LO_Model.at[i, 'date'].month and startdate.day == LO_Model.at[i, 'date'].day and (Pre_defined_Variables.CSflag == 0 or startdate.year == LO_Model.at[i, 'date'].year):
                X2 = 'SimDay1'
            else:
                X2 = LO_Model.at[i, 'date'].date()
            M_var.DayFlags[i] = X2
        M_var.PlsDay[i+2] = LO_FNs.PlsDay(M_var.DayFlags[i+2],M_var.DecTree_Relslevel[i+2],Pre_defined_Variables.PlsDay_Switch)
        M_var.Release_Level[i+2] = LO_FNs.Release_Level(M_var.Release_Level[i+1],M_var.Lake_Stage[i+1],TC_LONINO_df.at[i, 'Tributary_Condition'],M_var.PlsDay[i+2],M_var.Zone_Code[i+1],M_var.DecTree_Relslevel[i+2],Pre_defined_Variables.MaxQstgTrigger)
        if i >= 6:
            dh = M_var.Lake_Stage[i+1] - M_var.Lake_Stage[i-6]
            M_var.dh_7days[i+1] = dh
        M_var.ZoneCodeminus1Code[i+1] = LO_FNs.ZoneCodeminus1Code(M_var.Zone_Code[i+1],df_WSMs.at[i+1, 'WSM1'],df_WSMs.at[i+1, 'D0'],df_WSMs.at[i+1, 'D1'],df_WSMs.at[i+1, 'D2'],df_WSMs.at[i+1, 'D3'],df_WSMs.at[i+1, 'C'],df_WSMs.at[i+1, 'B'],df_WSMs.at[i+1, 'A'])
        M_var.ZoneCodeCode[i+1] = LO_FNs.ZoneCodeCode(M_var.Zone_Code[i+1],df_WSMs.at[i+1, 'WSM1'],df_WSMs.at[i+1, 'D0'],df_WSMs.at[i+1, 'D1'],df_WSMs.at[i+1, 'D2'],df_WSMs.at[i+1, 'D3'],df_WSMs.at[i+1, 'C'],df_WSMs.at[i+1, 'B'],df_WSMs.at[i+1, 'A'])
        M_var.Fraction_of_Zone_height[i+1] = LO_FNs.Fraction_of_Zone_height(M_var.Zone_Code[i+1],M_var.Lake_Stage[i+1],M_var.ZoneCodeminus1Code[i+1],M_var.ZoneCodeCode[i+1])
        M_var.ReLevelCode_1[i+2] = LO_FNs.ReLevelCode_1(M_var.Release_Level[i+2],Pre_defined_Variables.dstar_D1,Pre_defined_Variables.dstar_D2,Pre_defined_Variables.dstar_D3,Pre_defined_Variables.dstar_C,Pre_defined_Variables.dstar_B)
        M_var.ReLevelCode_2[i+2] = LO_FNs.ReLevelCode_2(M_var.Release_Level[i+2],Pre_defined_Variables.astar_D1,Pre_defined_Variables.astar_D2,Pre_defined_Variables.astar_D3,Pre_defined_Variables.astar_C,Pre_defined_Variables.astar_B)
        M_var.ReLevelCode_3_S80[i+2] = LO_FNs.ReLevelCode_3_S80(M_var.Release_Level[i+2],Pre_defined_Variables.bstar_S80_D1,Pre_defined_Variables.bstar_S80_D2,Pre_defined_Variables.bstar_S80_D3,Pre_defined_Variables.bstar_S80_C,Pre_defined_Variables.bstar_S80_B)
        M_var.Outlet2DS_Mult[i+2] = LO_FNs.Outlet2DS_Mult(Seasons.at[i, 'Season'],Seasons.at[i, 'Month'],M_var.dh_7days[i+1],M_var.ReLevelCode_1[i+2],M_var.Fraction_of_Zone_height[i+1],M_var.ReLevelCode_2[i+2],M_var.ReLevelCode_3_S80[i+2],Pre_defined_Variables.Opt_QregMult)
        M_var.Outlet2DS_Mult_2[i+2] = LO_FNs.Outlet2DS_Mult_2(LO_Model.at[i+2, 'date'].month,LO_Model.at[i+2, 'date'].day,M_var.PlsDay[i+2],M_var.Outlet2DS_Mult[i+2-M_var.PlsDay[i+2]],M_var.Outlet2DS_Mult[i+2],Pre_defined_Variables.Opt_QregMult)        
        M_var.Outlet2DSRS[i+2] = LO_FNs.Outlet2DSRS(M_var.Release_Level[i+2],Data.S80_RegRelRates.at[0, 'Zone_D1'],S80avgL1,Data.Pulses.at[M_var.PlsDay[i+2]-1 if M_var.PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-80_L1_%s'%Pre_defined_Variables.Schedule],M_var.Outlet2DS_Mult_2[i+2],Data.CE_SLE_turns.at[LO_Model.at[i+2, 'date'].year-Pre_defined_Variables.startyear, 'SLEturn'],Data.S80_RegRelRates.at[0, 'Zone_D2'],S80avgL2,Data.Pulses.at[M_var.PlsDay[i+2]-1 if M_var.PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-80_L2_%s'%Pre_defined_Variables.Schedule],Data.S80_RegRelRates.at[0, 'Zone_D3'],S80avgL3,Data.Pulses.at[M_var.PlsDay[i+2]-1 if M_var.PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-80_L3_%s'%Pre_defined_Variables.Schedule],Data.S80_RegRelRates.at[0, 'Zone_C'],Data.S80_RegRelRates.at[0, 'Zone_B'],Data.S80_RegRelRates.at[0, 'Zone_A'])
        M_var.Outlet2USRG1[i+2] = max(0,M_var.Outlet2DSRS[i+2]-LO_Model.at[i+2, 'C44RO'])
        M_var.Sum_Outlet2USRG1[i+2] = LO_FNs.Sum_Outlet2USRG1(LO_Model.at[i+2, 'date'].day,M_var.Outlet2USRG1[i+2])
        M_var.Outlet2DSBS[i+2] = LO_FNs.Outlet2DSBS(M_var.Release_Level[i+2],M_var.Sum_Outlet2USRG1[i+2],VLOOKUP1_c[i],Outlet2_baseflow,Pre_defined_Variables.Option_S80Baseflow)
        M_var.Outlet2USBK[i+2] = LO_FNs.Outlet2USBK(M_var.Lake_Stage[i+1],df_WSMs.at[i+1, 'D1'],M_var.Outlet2USRG[i+1],LO_Model.at[i+2, 'C44RO'],Data.SFWMM_Daily_Outputs.at[i+2, 'S308BK'],Pre_defined_Variables.Opt_S308,Pre_defined_Variables.S308BK_Const,Pre_defined_Variables.S308_BK_Thr)
        M_var.ROeast[i+2] = LO_Model.at[i+2, 'C44RO'] - M_var.Outlet2USBK[i+2]
        M_var.Outlet2USBS[i+2] = LO_FNs.Outlet2USBS(M_var.Outlet2DSBS[i+2],M_var.Outlet2USRG1[i+2],M_var.ROeast[i+2],Pre_defined_Variables.Option_S80Baseflow)
        M_var.Sum_Outlet2USBK[i+2] = LO_FNs.Sum_Outlet2USBK(LO_Model.at[i+2, 'date'].day,M_var.Outlet2USBK[i+2])
        M_var.Outlet2USRG_Code[i+2] = LO_FNs.Outlet2USRG_Code(M_var.Outlet2USRG1[i+2],M_var.Outlet2USBS[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S308RG'],Data.SFWMM_Daily_Outputs.at[i+2, 'STEST'],Pre_defined_Variables.Option_RegS77S308)
        if Model_Config.Sim_type == 0:
            M_var.Outlet2USRG[i+2] = LO_FNs.Outlet2USRG(M_var.Outlet2USRG_Code[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S308RG'],Data.SFWMM_Daily_Outputs.at[i+2, 'STEST'],Pre_defined_Variables.Opt_S308,Pre_defined_Variables.S308RG_Const)
        else:
            if M_var.Lake_Stage[i+1] >= 18:
                M_var.Outlet2USRG[i+2] = 7200
            elif M_var.Lake_Stage[i+1] <= 8:
                M_var.Outlet2USRG[i+2] = 0
            elif (Sim_Chla_S[i] <= Chla_1) and (date_rng_6[i+2].month in [1,2,3,4,11,12]):
                M_var.Outlet2USRG[i+2] = S308_DV[(date_rng_6[i+2].month)-1]
            elif (Sim_Chla_S[i] <= Chla_2) and (date_rng_6[i+2].month in [5,6,7,8,9,10]):
                M_var.Outlet2USRG[i+2] = S308_DV[(date_rng_6[i+2].month)-1]
            # elif M_var.Lake_Stage[i+1] > LO_EcoRecEnv['Upper bound'].iloc[i+1] :  # Soft constraint
            #     M_var.Outlet2USRG[i+2] = S308_DV[(date_rng_6[i+2].month)-1]  
            else:
                M_var.Outlet2USRG[i+2] = 0
        M_var.Outlet2DS[i+2] = LO_FNs.S80(M_var.ROeast[i+2],M_var.Outlet2USRG[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S80'],Pre_defined_Variables.S80_Const)
        M_var.ReLevelCode_3_S77[i+2] = LO_FNs.ReLevelCode_3_S77(M_var.Release_Level[i+2],Pre_defined_Variables.bstar_S77_D1,Pre_defined_Variables.bstar_S77_D2,Pre_defined_Variables.bstar_S77_D3,Pre_defined_Variables.bstar_S77_C,Pre_defined_Variables.bstar_S77_B)
        M_var.Outlet1US_Mult[i+2] = LO_FNs.Outlet1US_Mult(Seasons.at[i, 'Season'],Seasons.at[i, 'Month'],M_var.dh_7days[i+1],M_var.ReLevelCode_1[i+2],M_var.Fraction_of_Zone_height[i+1],M_var.ReLevelCode_2[i+2],M_var.ReLevelCode_3_S77[i+2],Pre_defined_Variables.Opt_QregMult)
        M_var.Outlet1US_Mult_2[i+2] = LO_FNs.Outlet1US_Mult_2(LO_Model.at[i+2, 'date'].month,LO_Model.at[i+2, 'date'].day,M_var.PlsDay[i+2],M_var.Outlet1US_Mult[i+2-M_var.PlsDay[i+2]],M_var.Outlet1US_Mult[i+2],Pre_defined_Variables.Opt_QregMult)
        M_var.Outlet1USRS[i+2] = LO_FNs.Outlet1USRS(M_var.Release_Level[i+2],Data.S77_RegRelRates.at[0, 'Zone_D1'],S77avgL1,Data.Pulses.at[M_var.PlsDay[i+2]-1 if M_var.PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-77_L1_%s'%Pre_defined_Variables.Schedule],M_var.Outlet1US_Mult_2[i+2],LO_Model.at[i+2, 'C43RO'],Data.CE_SLE_turns.at[LO_Model.at[i+2, 'date'].year-Pre_defined_Variables.startyear, 'CEturn'],Data.S77_RegRelRates.at[0, 'Zone_D2'],S77avgL2,Data.Pulses.at[M_var.PlsDay[i+2]-1 if M_var.PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-77_L2_%s'%Pre_defined_Variables.Schedule],M_var.Zone_Code[i+1],Data.S77_RegRelRates.at[0, 'Zone_D3'],S77avgL3,Data.Pulses.at[M_var.PlsDay[i+2]-1 if M_var.PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-77_L3_%s'%Pre_defined_Variables.Schedule],Data.S77_RegRelRates.at[0, 'Zone_C'],Data.S77_RegRelRates.at[0, 'Zone_B'],Data.S77_RegRelRates.at[0, 'Zone_A'],Pre_defined_Variables.Opt_Outlet1DSRG)
        M_var.Sum_Outlet1USRS[i+2] = LO_FNs.Sum_Outlet1USRS(LO_Model.at[i+2, 'date'].day,M_var.Outlet1USRS[i+2])
        M_var.Outlet1USBK[i+2] = LO_FNs.Outlet1USBK(M_var.Lake_Stage[i+1],M_var.Outlet1USRS[i+2],M_var.Outlet1USBSAP[i+1],M_var.Outlet1USEWS[i+1],LO_Model.at[i+2, 'C43RO'],Data.SFWMM_Daily_Outputs.at[i+2, 'S77BK'],Pre_defined_Variables.Outlet1USBK_Switch,Pre_defined_Variables.Outlet1USBK_Threshold)
        M_var.ROwest[i+2] = LO_Model.at[i+2, 'C43RO'] - M_var.Outlet1USBK[i+2]
        M_var.Outlet1DSBS[i+2] = LO_FNs.Outlet1DSBS(M_var.Release_Level[i+2],M_var.Sum_Outlet1USRS[i+2],VLOOKUP2_c[i],Outlet1_baseflow,Pre_defined_Variables.Option_S77Baseflow)
        M_var.Outlet1USBS[i+2] = LO_FNs.Outlet1USBS(M_var.Outlet1DSBS[i+2],M_var.Outlet1USRS[i+2],M_var.ROwest[i+2],Pre_defined_Variables.Option_S77Baseflow)
        #Define THC Class Normal or above
        if i < (n_rows-2):
            M_var.Post_Ap_Baseflow[i] = THC_Class(i,M_var.THC_Class_normal_or_above,M_var.Lake_O_Stage_AP,M_var.Lake_O_Schedule_Zone,M_var.LStgCorres,M_var.LowChance_Check,M_var.Outlet1USRS_AP,M_var.Outlet1USBS_AP,
              M_var.Outlet1USRS_Pre_AP_S77_Baseflow,M_var.Forecast_D_Sal,M_var.n30d_mavg,M_var.n30davgForecast,M_var.LORS08_bf_rel,M_var.LDS_LC6_1,M_var.S_O,M_var.All_4,
              M_var.Sabf,M_var.Swbf,M_var.Swbu,M_var.All_4andStage,M_var.All_4andStagein,M_var.P_AP_BF_Stg,M_var.Logic_test_1,M_var.Post_Ap_Baseflow,M_var.Outlet1USRSplusPreAPS77bsf,
              M_var.AndEstNeedsLakeWater,M_var.AndLowChance61stagelessth11,M_var.ATHCnora,M_var.Choose_PAPEWS_1,M_var.Choose_PAPEWS_2,M_var.Post_AP_EWS,
              M_var.Post_AP_Baseflow_EWS_cfs,AdapProt_df,M_var.Lake_Stage,M_var.Zone_Code,df_WSMs,Targ_Stg_df,M_var.Outlet1USRS,M_var.Outlet1USBS,Data.Estuary_needs_water,
              Choose_1,M_var.WSM_Zone)['Post_Ap_Baseflow']
            M_var.Post_AP_EWS[i] = THC_Class(i,M_var.THC_Class_normal_or_above,M_var.Lake_O_Stage_AP,M_var.Lake_O_Schedule_Zone,M_var.LStgCorres,M_var.LowChance_Check,M_var.Outlet1USRS_AP,M_var.Outlet1USBS_AP,
              M_var.Outlet1USRS_Pre_AP_S77_Baseflow,M_var.Forecast_D_Sal,M_var.n30d_mavg,M_var.n30davgForecast,M_var.LORS08_bf_rel,M_var.LDS_LC6_1,M_var.S_O,M_var.All_4,
              M_var.Sabf,M_var.Swbf,M_var.Swbu,M_var.All_4andStage,M_var.All_4andStagein,M_var.P_AP_BF_Stg,M_var.Logic_test_1,M_var.Post_Ap_Baseflow,M_var.Outlet1USRSplusPreAPS77bsf,
              M_var.AndEstNeedsLakeWater,M_var.AndLowChance61stagelessth11,M_var.ATHCnora,M_var.Choose_PAPEWS_1,M_var.Choose_PAPEWS_2,M_var.Post_AP_EWS,
              M_var.Post_AP_Baseflow_EWS_cfs,AdapProt_df,M_var.Lake_Stage,M_var.Zone_Code,df_WSMs,Targ_Stg_df,M_var.Outlet1USRS,M_var.Outlet1USBS,Data.Estuary_needs_water,
              Choose_1,M_var.WSM_Zone)['Post_AP_EWS']
        M_var.Outlet1USBSAP[i+2] = LO_FNs.Outlet1USBSAP(M_var.Outlet1USBS[i+2],M_var.Post_Ap_Baseflow[i],Pre_defined_Variables.Opt_AdapProt)
        M_var.Outlet1USEWS[i+2] = LO_FNs.Outlet1USEWS(M_var.Post_AP_EWS[i],Data.SFWMM_Daily_Outputs.at[i+2, 'CAEST'],Pre_defined_Variables.Outlet1USEWS_Switch,Pre_defined_Variables.Opt_AdapProt)
        if Model_Config.Sim_type == 0:
            M_var.Outlet1USREG[i+2] = LO_FNs.Outlet1USREG(M_var.Outlet1USRS[i+2],M_var.Outlet1USBSAP[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S77RG'],Pre_defined_Variables.Outlet1USREG_Switch,Pre_defined_Variables.Option_RegS77S308)
        else:
            if M_var.Lake_Stage[i+1] >= 18:
                M_var.Outlet1USREG[i+2] = 7800
            elif M_var.Lake_Stage[i+1] <= 8:
                M_var.Outlet1USREG[i+2] = 0
            elif (Sim_Chla_S[i] <= Chla_1) and (date_rng_6[i+2].month in [1,2,3,4,11,12]):
                M_var.Outlet1USREG[i+2] = S77_DV[(date_rng_6[i+2].month)-1]
            elif (Sim_Chla_S[i] <= Chla_2) and (date_rng_6[i+2].month in [5,6,7,8,9,10]):
                M_var.Outlet1USREG[i+2] = S77_DV[(date_rng_6[i+2].month)-1]
            # elif M_var.Lake_Stage[i+1] > LO_EcoRecEnv['Upper bound'].iloc[i+1] :  # Soft constraint
            #     M_var.Outlet1USREG[i+2] = S77_DV[(date_rng_6[i+2].month)-1]  
            else:
                M_var.Outlet1USREG[i+2] = 0
        M_var.Outlet1DS[i+2] = LO_FNs.Outlet1DS(M_var.Outlet1USREG[i+2],M_var.Outlet1USEWS[i+2],M_var.ROwest[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S79'],Pre_defined_Variables.Outlet1DS_Switch)
        M_var.TotRegEW[i+2] = (M_var.Outlet1USREG[i+2] + M_var.Outlet2USRG[i+2])*1.9835
        M_var.Choose_WCA[i+2] = LO_FNs.Choose_WCA(Data.SFWMM_Daily_Outputs.at[i+2, 'RegWCA'],Pre_defined_Variables.Option_RegWCA,Pre_defined_Variables.Constant_RegWCA)
        M_var.RegWCA[i+2] = min(Pre_defined_Variables.MaxCap_RegWCA , Pre_defined_Variables.Multiplier_RegWCA*M_var.Choose_WCA[i+2])
        M_var.Choose_L8C51[i+2] = LO_FNs.Choose_L8C51(Data.SFWMM_Daily_Outputs.at[i+2, 'RegL8C51'],Pre_defined_Variables.Option_RegL8C51,Pre_defined_Variables.Constant_RegL8C51)
        M_var.RegL8C51[i+2] = min(Pre_defined_Variables.MaxCap_RegL8C51 , Pre_defined_Variables.Multiplier_RegL8C51*M_var.Choose_L8C51[i+2])
        
        # Total Reg South
        # if Model_Config.Sim_type == 0:
        #     M_var.TotRegSo[i+2] = (M_var.RegWCA[i+2] + M_var.RegL8C51[i+2]) * 1.9835
        # else:
        #     if M_var.Lake_Stage[i+1] < LO_EcoRecEnv['Lower bound'].iloc[i+1]:
        #         M_var.TotRegSo[i+2] = 0
        #     else:
        #         M_var.TotRegSo[i+2] = RegSouth_DV[(date_rng_6[i+2].month)-1]* 1.9835
        
        M_var.TotRegSo[i+2] = (M_var.RegWCA[i+2] + M_var.RegL8C51[i+2]) * 1.9835

        
        M_var.Stage2ar[i+2] = Stg_Sto_Ar.stg2ar(M_var.Lake_Stage[i+1],0)
        M_var.Stage2marsh[i+2] = Stg_Sto_Ar.stg2mar(M_var.Lake_Stage[i+1],0)
        M_var.RF[i+2] = Data.RF_Vol.at[i+2, 'RF_acft']
        M_var.ET[i+2] = LO_FNs.ET(Data.SFWMM_Daily_Outputs.at[i+2, 'et_dry'],M_var.Stage2ar[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'et_litoral'],M_var.Stage2marsh[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'et_open'],Data.ET_Vol.at[i+2, 'ETVol_acft'],Pre_defined_Variables.ET_Switch)
        M_var.Choose_WSA_1[i+2] = LO_FNs.Choose_WSA_1(df_WSMs.at[i+2, 'WSM1'],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSAtrig2,Pre_defined_Variables.WSAoff2)   
        M_var.Choose_WSA_2[i+2] = LO_FNs.Choose_WSA_2(df_WSMs.at[i+2, 'WSM1'],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSAtrig1,Pre_defined_Variables.WSAoff1)    
        M_var.WSA_MIA[i+2] = LO_FNs.WSA_MIA(WCA_Stages_df.at[i, 'Are WCA stages too low?'],TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],M_var.Lake_Stage[i+1],M_var.Choose_WSA_1[i+2],Data.EAA_MIA_RUNOFF.at[i, 'MIA'],Data.EAA_MIA_RUNOFF.at[i, 'S3PMP'],M_var.Choose_WSA_2[i+2],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSA_THC,Pre_defined_Variables.MIAcap2,Pre_defined_Variables.MIAcap1)
        M_var.WSA_NNR[i+2] = LO_FNs.WSA_NNR(WCA_Stages_df.at[i, 'Are WCA stages too low?'],TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],M_var.Lake_Stage[i+1],M_var.Choose_WSA_1[i+2],Data.EAA_MIA_RUNOFF.at[i, 'NNR'],Data.EAA_MIA_RUNOFF.at[i, 'S2PMP'],M_var.Choose_WSA_2[i+2],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSA_THC,Pre_defined_Variables.NNRcap2,Pre_defined_Variables.NNRcap1)
        M_var.DSto[i+2] = M_var.NI_Supply[i+2] + M_var.RF[i+2] - M_var.ET[i+2] + 1.9835*(M_var.Outlet2USBK[i+2]\
                     + M_var.Outlet1USBK[i+2] + M_var.WSA_MIA[i+2] + M_var.WSA_NNR[i+2]\
                     - M_var.Outlet1USEWS[i+2]) - M_var.TotRegEW[i+2] - M_var.TotRegSo[i+2] + Storage_dev[i+2]
        M_var.Storage[i+2] = LO_FNs.Storage(M_var.DayFlags[i+2],M_var.Storage[i],StartStorage,M_var.Storage[i+1],M_var.DSto[i+2])
        M_var.Lake_Stage[i+2] = LO_FNs.Lake_Stage(Stg_Sto_Ar.stg2sto(M_var.Storage[i+2],1),Data.SFWMM_Daily_Outputs.at[i+2, 'EOD Stg(ft,NGVD)'],Pre_defined_Variables.Option_Stage)
        # if M_var.Lake_Stage[i+2] >= 18:
        #     Flood[i+2] = 1
        # else:
        #     Flood[i+2] = 0
        #######################################################################################################################
        Q_O[i+2] = (M_var.Outlet1USEWS[i+2] *0.028316847 + ((M_var.TotRegEW[i+2] + M_var.TotRegSo[i+2])/70.0456)) * 3600 * 24
        
        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48 #m3/d
            Q_O_M[i] = Q_O[i]
            External_NO_M[i] = External_NO[i] + Q_I_M[i] * S65E_NO[i]            
            External_Chla_M[i] = External_Chla[i] + Q_I_M[i] * S65E_Chla[i]    
            L_ext_M[i] = L_ext[i] + Q_I_M[i] * S65_P[i]            

        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48 #m3/d
            Q_I_M[i] = Q_I[i]
            External_NO_M[i] = External_NO[i]         
            External_Chla_M[i] = External_Chla[i]
            L_ext_M[i] = L_ext[i]

        Q_N2S[i] = (Q_I_M[i]*1 + Q_O_M[i]*0)   
    
        t = np.linspace(1,2,num = 2)
        Nit_R_N_NO[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0_NO,kni_NO,NH4_N_Obs[i],K_NH_NO,DO[i],KDO_NO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_N_T_NO[i] = Nit_R_N_NO[i]*Theta_ni**(Temp[i]-20) 
        Nit_R_S_NO[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0_NO,kni_NO,NH4_S_Obs[i],K_NH_NO,DO[i],KDO_NO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_S_T_NO[i] = Nit_R_S_NO[i]*Theta_ni**(Temp[i]-20) 
    
        Denit_R_N_NO[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0_NO,kden_NO,NO_N[i],K_Nitr_NO,DO[i],KDO_NO,DO_Cr_d,DO_Opt_d)
        Denit_R_N_T_NO[i] = Denit_R_N_NO[i]*Theta_deni**(Temp[i]-20) 
        Denit_R_S_NO[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0_NO,kden_NO,NO_S[i],K_Nitr_NO,DO[i],KDO_NO,DO_Cr_d,DO_Opt_d)
        Denit_R_S_T_NO[i] = Denit_R_S_NO[i]*Theta_deni**(Temp[i]-20) 
    
    
        fT_NO[i] = f_T_alt1(Temp[i],T_opt_NO,T_min_NO,T_max_NO)
        fL_NO[i] = f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw_NO,Kc_NO,Sim_Chla[i],z[i],K1_NO,K2_NO)

        DfEq_Res_N = odeint(NOx_N_DiffEq,NO_N[i],t,args=(External_NO[i],Atm_Deposition_N,Q_N2S[i],Nit_R_N_T_NO[i],Denit_R_N_T_NO[i],NH4_N_Obs[i],volume_N[i],G_max_NO,fT_NO[i],fL_NO[i],f_P_N_NO[i],f_N_N_NO[i],K_NH_NO,K_TN_NO,YNOChla_NO,Sim_Chla_N[i],S_NO_NO,Area_N,NO_Temp_Adj[i],))
        NO_N[i+1] = DfEq_Res_N[:,0][1]
        DfEq_Res_S = odeint(NOx_S_DiffEq,NO_S[i],t,args=(Atm_Deposition_S,Q_N2S[i],Q_O_M[i],NO_N[i],Nit_R_S_T_NO[i],Denit_R_S_T_NO[i],NH4_S_Obs[i],volume_S[i],G_max_NO,fT_NO[i],fL_NO[i],f_P_S_NO[i],f_N_S_NO[i],K_NH_NO,K_TN_NO,YNOChla_NO,Sim_Chla_S[i],S_NO_NO,Area_S,NO_Temp_Adj[i],))
        NO_S[i+1] = DfEq_Res_S[:,0][1]
        NO_MEAN[i+1] = (NO_N[i+1] + NO_S[i+1])*0.5
          
        
        NO_Load_Cal[i] = M_var.Outlet1USREG[i]*0.028316847*3600*24*NO_S[i] #mg
        NO_Load_StL[i] = M_var.Outlet2USRG[i]*0.028316847*3600*24*NO_S[i] #mg
        NO_Load_South[i] = M_var.TotRegSo[i]*1233.48 *NO_S[i] #mg
    
        
        ##### Chla 
        fT_Chla[i] = f_T__Chla_alt1(Temp[i],T_opt_Chla,T_min_Chla,T_max_Chla,Nitro_Model_Output['date'].iloc[i].month)
        
        fL_Chla[i] = f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw_Chla,Kc_Chla,Sim_Chla[i],z[i],K1_Chla,K2_Chla)*1.2 if Nitro_Model_Output['date'].iloc[i].month in (6,7,8,9,10) else f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw_Chla,Kc_Chla,Sim_Chla[i],z[i],K1_Chla,K2_Chla)*1

        Sim_Chla_N[i+1] = Chla_N_alt1(External_Chla_M[i], Q_N2S[i], Sim_Chla_N[i], vc, z[i], K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla[i], fL_Chla[i], f_P_N_Chla[i], f_N_N_Chla[i], volume_N[i])
        Sim_Chla_S[i+1] = Chla_S_alt1(Q_N2S[i], Q_O_M[i], Sim_Chla_N[i], Sim_Chla_S[i], vc, z[i], K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla[i], fL_Chla[i], f_P_S_Chla[i], f_N_S_Chla[i], volume_S[i])
        Sim_Chla[i+1] = (Sim_Chla_N[i+1] + Sim_Chla_S[i+1])/2
    
        Chla_Load_Cal[i] = M_var.Outlet1USREG[i]*0.028316847*3600*24*Sim_Chla_S[i] #mg
        Chla_Load_StL[i] = M_var.Outlet2USRG[i]*0.028316847*3600*24*Sim_Chla_S[i] #mg
        Chla_Load_South[i] = M_var.TotRegSo[i]*1233.48 *Sim_Chla_S[i] #mg


        Stage2ar[i+2] = Stg_Sto_Ar.stg2ar(M_var.Lake_Stage[i+2],0)
        LO_WD[i] = M_var.Lake_Stage[i]*0.3048 - LO_BL
        Lake_O_Storage_N[i] = M_var.Storage[i] * TP_Variables.N_Per * 4046.85642 * 0.305 #m3
        Lake_O_Storage_S[i] = M_var.Storage[i] * TP_Variables.S_Per * 4046.85642 * 0.305 #m3
        Lake_O_A_N[i] = Stage2ar[i] * TP_Variables.N_Per * 4046.85642 #m2
        Lake_O_A_S[i] = Stage2ar[i] * TP_Variables.S_Per * 4046.85642 #m2
        Lake_O_A_M_N[i] = Lake_O_A_N[i] * TP_Variables.A_Mud_N/(TP_Variables.A_Mud_N+TP_Variables.A_Sand_N+TP_Variables.A_Rock_N+TP_Variables.A_Peat_N)
        Lake_O_A_S_N[i] = Lake_O_A_N[i] * TP_Variables.A_Sand_N/(TP_Variables.A_Mud_N+TP_Variables.A_Sand_N+TP_Variables.A_Rock_N+TP_Variables.A_Peat_N)
        Lake_O_A_R_N[i] = Lake_O_A_N[i] * TP_Variables.A_Rock_N/(TP_Variables.A_Mud_N+TP_Variables.A_Sand_N+TP_Variables.A_Rock_N+TP_Variables.A_Peat_N)
        Lake_O_A_P_N[i] = Lake_O_A_N[i] * TP_Variables.A_Peat_N/(TP_Variables.A_Mud_N+TP_Variables.A_Sand_N+TP_Variables.A_Rock_N+TP_Variables.A_Peat_N)
        Lake_O_A_M_S[i] = Lake_O_A_S[i] * TP_Variables.A_Mud_S/(TP_Variables.A_Mud_S+TP_Variables.A_Sand_S+TP_Variables.A_Rock_S+TP_Variables.A_Peat_S)
        Lake_O_A_S_S[i] = Lake_O_A_S[i] * TP_Variables.A_Sand_S/(TP_Variables.A_Mud_S+TP_Variables.A_Sand_S+TP_Variables.A_Rock_S+TP_Variables.A_Peat_S)
        Lake_O_A_R_S[i] = Lake_O_A_S[i] * TP_Variables.A_Rock_S/(TP_Variables.A_Mud_S+TP_Variables.A_Sand_S+TP_Variables.A_Rock_S+TP_Variables.A_Peat_S)
        Lake_O_A_P_S[i] = Lake_O_A_S[i] * TP_Variables.A_Peat_S/(TP_Variables.A_Mud_S+TP_Variables.A_Sand_S+TP_Variables.A_Rock_S+TP_Variables.A_Peat_S)
        
        DIP_Lake_N[i] = TP_MBFR.DIP_Lake(TP_Lake_N[i])
        DIP_Lake_S[i] = TP_MBFR.DIP_Lake(TP_Lake_S[i])
        
        
        v_settle_N_c[i] = (R*g*d_c**2)/(C_1_c*nu_d[i]+(0.75*C_2_c*R*g*d_c**3)**0.5)
        v_settle_N_s[i] = (R*g*d_s**2)/(C_1_s*nu_d[i]+(0.75*C_2_s*R*g*d_s**3)**0.5)
        v_settle_N[i] = v_settle_N_c[i]*((TP_Variables.A_Mud_N+TP_Variables.A_Peat_N)/TP_Variables.A_N) + v_settle_N_s[i]*((TP_Variables.A_Sand_N + TP_Variables.A_Rock_N)/TP_Variables.A_N)
        
        v_settle_S_c[i] = (R*g*d_c**2)/(C_1_c*nu_d[i]+(0.75*C_2_c*R*g*d_c**3)**0.5)
        v_settle_S_s[i] = (R*g*d_s**2)/(C_1_s*nu_d[i]+(0.75*C_2_s*R*g*d_s**3)**0.5)
        v_settle_S[i] = v_settle_S_c[i]*((TP_Variables.A_Mud_S+TP_Variables.A_Peat_S)/TP_Variables.A_S) + v_settle_S_s[i]*((TP_Variables.A_Sand_S + TP_Variables.A_Rock_S)/TP_Variables.A_S)
        
        
        J_des_M_N[i] = TP_MBFR.Des_flux(Γ_M_N[i],Mass_sed_M_N,TP_Variables.K_des_M)
        J_des_S_N[i] = TP_MBFR.Des_flux(Γ_S_N[i],Mass_sed_S_N,TP_Variables.K_des_S)
        J_des_R_N[i] = TP_MBFR.Des_flux(Γ_R_N[i],Mass_sed_R_N,TP_Variables.K_des_R)
        J_des_P_N[i] = TP_MBFR.Des_flux(Γ_P_N[i],Mass_sed_P_N,TP_Variables.K_des_P)
        J_des_M_S[i] = TP_MBFR.Des_flux(Γ_M_S[i],Mass_sed_M_S,TP_Variables.K_des_M)
        J_des_S_S[i] = TP_MBFR.Des_flux(Γ_S_S[i],Mass_sed_S_S,TP_Variables.K_des_S)
        J_des_R_S[i] = TP_MBFR.Des_flux(Γ_R_S[i],Mass_sed_R_S,TP_Variables.K_des_R)
        J_des_P_S[i] = TP_MBFR.Des_flux(Γ_P_S[i],Mass_sed_P_S,TP_Variables.K_des_P)
        
        J_ads_M_N[i] = TP_MBFR.Ads_flux(DIP_pore_M_N[i],Γ_M_N[i],Mass_sed_M_N,TP_Variables.K_ads_M,TP_Variables.Γ_inf)
        J_ads_S_N[i] = TP_MBFR.Ads_flux(DIP_pore_S_N[i],Γ_S_N[i],Mass_sed_S_N,TP_Variables.K_ads_S,TP_Variables.Γ_inf)
        J_ads_R_N[i] = TP_MBFR.Ads_flux(DIP_pore_R_N[i],Γ_R_N[i],Mass_sed_R_N,TP_Variables.K_ads_R,TP_Variables.Γ_inf)
        J_ads_P_N[i] = TP_MBFR.Ads_flux(DIP_pore_P_N[i],Γ_P_N[i],Mass_sed_P_N,TP_Variables.K_ads_P,TP_Variables.Γ_inf)
        J_ads_M_S[i] = TP_MBFR.Ads_flux(DIP_pore_M_S[i],Γ_M_S[i],Mass_sed_M_S,TP_Variables.K_ads_M,TP_Variables.Γ_inf)
        J_ads_S_S[i] = TP_MBFR.Ads_flux(DIP_pore_S_S[i],Γ_S_S[i],Mass_sed_S_S,TP_Variables.K_ads_S,TP_Variables.Γ_inf)
        J_ads_R_S[i] = TP_MBFR.Ads_flux(DIP_pore_R_S[i],Γ_R_S[i],Mass_sed_R_S,TP_Variables.K_ads_R,TP_Variables.Γ_inf)
        J_ads_P_S[i] = TP_MBFR.Ads_flux(DIP_pore_P_S[i],Γ_P_S[i],Mass_sed_P_S,TP_Variables.K_ads_P,TP_Variables.Γ_inf)
        
        J_sedburial_M_N[i] = TP_MBFR.Sed_burial_flux(P_sed_M_N[i],TP_Variables.Bulk_density_M,TP_Variables.A_Mud_N,TP_Variables.v_burial_M,TP_Variables.Per_H2O_M)
        J_sedburial_S_N[i] = TP_MBFR.Sed_burial_flux(P_sed_S_N[i],TP_Variables.Bulk_density_S,TP_Variables.A_Sand_N,TP_Variables.v_burial_S,TP_Variables.Per_H2O_S)
        J_sedburial_R_N[i] = TP_MBFR.Sed_burial_flux(P_sed_R_N[i],TP_Variables.Bulk_density_R,TP_Variables.A_Rock_N,TP_Variables.v_burial_R,TP_Variables.Per_H2O_R)
        J_sedburial_P_N[i] = TP_MBFR.Sed_burial_flux(P_sed_P_N[i],TP_Variables.Bulk_density_P,TP_Variables.A_Peat_N,TP_Variables.v_burial_P,TP_Variables.Per_H2O_P)
        J_sedburial_M_S[i] = TP_MBFR.Sed_burial_flux(P_sed_M_S[i],TP_Variables.Bulk_density_M,TP_Variables.A_Mud_S,TP_Variables.v_burial_M,TP_Variables.Per_H2O_M)
        J_sedburial_S_S[i] = TP_MBFR.Sed_burial_flux(P_sed_S_S[i],TP_Variables.Bulk_density_S,TP_Variables.A_Sand_S,TP_Variables.v_burial_S,TP_Variables.Per_H2O_S)
        J_sedburial_R_S[i] = TP_MBFR.Sed_burial_flux(P_sed_R_S[i],TP_Variables.Bulk_density_R,TP_Variables.A_Rock_S,TP_Variables.v_burial_R,TP_Variables.Per_H2O_R)
        J_sedburial_P_S[i] = TP_MBFR.Sed_burial_flux(P_sed_P_S[i],TP_Variables.Bulk_density_P,TP_Variables.A_Peat_S,TP_Variables.v_burial_P,TP_Variables.Per_H2O_P)
        
        Sed_Resusp_M_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_M_N[i]*TP_Variables.A_Mud_N if W_SS[i] > Crtcl_ShStr else 0 #mg
        Sed_Resusp_S_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_S_N[i]*TP_Variables.A_Sand_N if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_R_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_R_N[i]*TP_Variables.A_Rock_N if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_P_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_P_N[i]*TP_Variables.A_Peat_N if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_M_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_M_S[i]*TP_Variables.A_Mud_S if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_S_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_S_S[i]*TP_Variables.A_Sand_S if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_R_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_R_S[i]*TP_Variables.A_Rock_S if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_P_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_P_S[i]*TP_Variables.A_Peat_S if W_SS[i] > Crtcl_ShStr else 0
        
        
        J_Γburial_M_N[i] = TP_MBFR.Sor_P_burialflux(Γ_M_N[i],TP_Variables.Bulk_density_M,TP_Variables.A_Mud_N,TP_Variables.v_burial_M,TP_Variables.Per_H2O_M)
        J_Γburial_S_N[i] = TP_MBFR.Sor_P_burialflux(Γ_S_N[i],TP_Variables.Bulk_density_S,TP_Variables.A_Sand_N,TP_Variables.v_burial_S,TP_Variables.Per_H2O_S)
        J_Γburial_R_N[i] = TP_MBFR.Sor_P_burialflux(Γ_R_N[i],TP_Variables.Bulk_density_R,TP_Variables.A_Rock_N,TP_Variables.v_burial_R,TP_Variables.Per_H2O_R)
        J_Γburial_P_N[i] = TP_MBFR.Sor_P_burialflux(Γ_P_N[i],TP_Variables.Bulk_density_P,TP_Variables.A_Peat_N,TP_Variables.v_burial_P,TP_Variables.Per_H2O_P)
        J_Γburial_M_S[i] = TP_MBFR.Sor_P_burialflux(Γ_M_S[i],TP_Variables.Bulk_density_M,TP_Variables.A_Mud_S,TP_Variables.v_burial_M,TP_Variables.Per_H2O_M)
        J_Γburial_S_S[i] = TP_MBFR.Sor_P_burialflux(Γ_S_S[i],TP_Variables.Bulk_density_S,TP_Variables.A_Sand_S,TP_Variables.v_burial_S,TP_Variables.Per_H2O_S)
        J_Γburial_R_S[i] = TP_MBFR.Sor_P_burialflux(Γ_R_S[i],TP_Variables.Bulk_density_R,TP_Variables.A_Rock_S,TP_Variables.v_burial_R,TP_Variables.Per_H2O_R)
        J_Γburial_P_S[i] = TP_MBFR.Sor_P_burialflux(Γ_P_S[i],TP_Variables.Bulk_density_P,TP_Variables.A_Peat_S,TP_Variables.v_burial_P,TP_Variables.Per_H2O_P)
        
        Γ_M_N[i+1] = TP_MBFR.Sor_P_conc(J_ads_M_N[i],J_des_M_N[i],J_Γburial_M_N[i],Γ_M_N[i],Mass_sed_M_N) if TP_MBFR.Sor_P_conc(J_ads_M_N[i],J_des_M_N[i],J_Γburial_M_N[i],Γ_M_N[i],Mass_sed_M_N) > 0 else 0
        Γ_S_N[i+1] = TP_MBFR.Sor_P_conc(J_ads_S_N[i],J_des_S_N[i],J_Γburial_S_N[i],Γ_S_N[i],Mass_sed_S_N) if TP_MBFR.Sor_P_conc(J_ads_S_N[i],J_des_S_N[i],J_Γburial_S_N[i],Γ_S_N[i],Mass_sed_S_N) > 0 else 0 
        Γ_R_N[i+1] = TP_MBFR.Sor_P_conc(J_ads_R_N[i],J_des_R_N[i],J_Γburial_R_N[i],Γ_R_N[i],Mass_sed_R_N) if TP_MBFR.Sor_P_conc(J_ads_R_N[i],J_des_R_N[i],J_Γburial_R_N[i],Γ_R_N[i],Mass_sed_R_N) > 0 else 0 
        Γ_P_N[i+1] = TP_MBFR.Sor_P_conc(J_ads_P_N[i],J_des_P_N[i],J_Γburial_P_N[i],Γ_P_N[i],Mass_sed_P_N) if TP_MBFR.Sor_P_conc(J_ads_P_N[i],J_des_P_N[i],J_Γburial_P_N[i],Γ_P_N[i],Mass_sed_P_N) > 0 else 0 
        Γ_M_S[i+1] = TP_MBFR.Sor_P_conc(J_ads_M_S[i],J_des_M_S[i],J_Γburial_M_S[i],Γ_M_S[i],Mass_sed_M_S) if TP_MBFR.Sor_P_conc(J_ads_M_S[i],J_des_M_S[i],J_Γburial_M_S[i],Γ_M_S[i],Mass_sed_M_S) > 0 else 0 
        Γ_S_S[i+1] = TP_MBFR.Sor_P_conc(J_ads_S_S[i],J_des_S_S[i],J_Γburial_S_S[i],Γ_S_S[i],Mass_sed_S_S) if TP_MBFR.Sor_P_conc(J_ads_S_S[i],J_des_S_S[i],J_Γburial_S_S[i],Γ_S_S[i],Mass_sed_S_S) > 0 else 0 
        Γ_R_S[i+1] = TP_MBFR.Sor_P_conc(J_ads_R_S[i],J_des_R_S[i],J_Γburial_R_S[i],Γ_R_S[i],Mass_sed_R_S) if TP_MBFR.Sor_P_conc(J_ads_R_S[i],J_des_R_S[i],J_Γburial_R_S[i],Γ_R_S[i],Mass_sed_R_S) > 0 else 0 
        Γ_P_S[i+1] = TP_MBFR.Sor_P_conc(J_ads_P_S[i],J_des_P_S[i],J_Γburial_P_S[i],Γ_P_S[i],Mass_sed_P_S) if TP_MBFR.Sor_P_conc(J_ads_P_S[i],J_des_P_S[i],J_Γburial_P_S[i],Γ_P_S[i],Mass_sed_P_S) > 0 else 0 
        
        J_decomp_M_N[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_M, P_sed_M_N[i], Mass_sed_M_N)
        J_decomp_S_N[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_S, P_sed_S_N[i], Mass_sed_S_N)
        J_decomp_R_N[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_R, P_sed_R_N[i], Mass_sed_R_N)
        J_decomp_P_N[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_P, P_sed_P_N[i], Mass_sed_P_N)
        J_decomp_M_S[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_M, P_sed_M_S[i], Mass_sed_M_S)
        J_decomp_S_S[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_S, P_sed_S_S[i], Mass_sed_S_S)
        J_decomp_R_S[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_R, P_sed_R_S[i], Mass_sed_R_S)
        J_decomp_P_S[i] = TP_MBFR.J_decomp(TP_Variables.K_decomp_P, P_sed_P_S[i], Mass_sed_P_S)
        
        P_sed_M_N[i+1] = TP_MBFR.P_sed(Lake_O_A_M_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.K_decomp_M,v_settle_N[i]) - Sed_Resusp_M_N[i]/Mass_sed_M_N if TP_MBFR.P_sed(Lake_O_A_M_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.K_decomp_M,v_settle_N[i]) - Sed_Resusp_M_N[i]/Mass_sed_M_N > 0 else 0
        P_sed_S_N[i+1] = TP_MBFR.P_sed(Lake_O_A_S_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.K_decomp_S,v_settle_N[i]) - Sed_Resusp_S_N[i]/Mass_sed_S_N if TP_MBFR.P_sed(Lake_O_A_S_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.K_decomp_S,v_settle_N[i]) - Sed_Resusp_S_N[i]/Mass_sed_S_N > 0 else 0
        P_sed_R_N[i+1] = TP_MBFR.P_sed(Lake_O_A_R_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.K_decomp_R,v_settle_N[i]) - Sed_Resusp_R_N[i]/Mass_sed_R_N if TP_MBFR.P_sed(Lake_O_A_R_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.K_decomp_R,v_settle_N[i]) - Sed_Resusp_R_N[i]/Mass_sed_R_N > 0 else 0
        P_sed_P_N[i+1] = TP_MBFR.P_sed(Lake_O_A_P_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.K_decomp_P,v_settle_N[i]) - Sed_Resusp_P_N[i]/Mass_sed_P_N if TP_MBFR.P_sed(Lake_O_A_P_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.K_decomp_P,v_settle_N[i]) - Sed_Resusp_P_N[i]/Mass_sed_P_N > 0 else 0
        P_sed_M_S[i+1] = TP_MBFR.P_sed(Lake_O_A_M_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.K_decomp_M,v_settle_S[i]) - Sed_Resusp_M_S[i]/Mass_sed_M_S if TP_MBFR.P_sed(Lake_O_A_M_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.K_decomp_M,v_settle_S[i]) - Sed_Resusp_M_S[i]/Mass_sed_M_S > 0 else 0
        P_sed_S_S[i+1] = TP_MBFR.P_sed(Lake_O_A_S_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.K_decomp_S,v_settle_S[i]) - Sed_Resusp_S_S[i]/Mass_sed_S_S if TP_MBFR.P_sed(Lake_O_A_S_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.K_decomp_S,v_settle_S[i]) - Sed_Resusp_S_S[i]/Mass_sed_S_S > 0 else 0
        P_sed_R_S[i+1] = TP_MBFR.P_sed(Lake_O_A_R_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.K_decomp_R,v_settle_S[i]) - Sed_Resusp_R_S[i]/Mass_sed_R_S if TP_MBFR.P_sed(Lake_O_A_R_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.K_decomp_R,v_settle_S[i]) - Sed_Resusp_R_S[i]/Mass_sed_R_S > 0 else 0
        P_sed_P_S[i+1] = TP_MBFR.P_sed(Lake_O_A_P_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.K_decomp_P,v_settle_S[i]) - Sed_Resusp_P_S[i]/Mass_sed_P_S if TP_MBFR.P_sed(Lake_O_A_P_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.K_decomp_P,v_settle_S[i]) - Sed_Resusp_P_S[i]/Mass_sed_P_S > 0 else 0
        
        DIP_pore_M_N[i+1] = TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_N[i],DIP_Lake_N[i],J_des_M_N[i],J_ads_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.v_diff_M,TP_Variables.A_Mud_N,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) if TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_N[i],DIP_Lake_N[i],J_des_M_N[i],J_ads_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.v_diff_M,TP_Variables.A_Mud_N,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) > 0 else 0
        DIP_pore_S_N[i+1] = TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_N[i],DIP_Lake_N[i],J_des_S_N[i],J_ads_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.v_diff_S,TP_Variables.A_Sand_N,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) if TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_N[i],DIP_Lake_N[i],J_des_S_N[i],J_ads_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.v_diff_S,TP_Variables.A_Sand_N,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) > 0 else 0
        DIP_pore_R_N[i+1] = TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_N[i],DIP_Lake_N[i],J_des_R_N[i],J_ads_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.v_diff_R,TP_Variables.A_Rock_N,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) if TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_N[i],DIP_Lake_N[i],J_des_R_N[i],J_ads_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.v_diff_R,TP_Variables.A_Rock_N,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) > 0 else 0
        DIP_pore_P_N[i+1] = TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_N[i],DIP_Lake_N[i],J_des_P_N[i],J_ads_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.v_diff_P,TP_Variables.A_Peat_N,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) if TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_N[i],DIP_Lake_N[i],J_des_P_N[i],J_ads_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.v_diff_P,TP_Variables.A_Peat_N,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) > 0 else 0
        DIP_pore_M_S[i+1] = TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_S[i],DIP_Lake_S[i],J_des_M_S[i],J_ads_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.v_diff_M,TP_Variables.A_Mud_S,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) if TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_S[i],DIP_Lake_S[i],J_des_M_S[i],J_ads_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.v_diff_M,TP_Variables.A_Mud_S,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) > 0 else 0
        DIP_pore_S_S[i+1] = TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_S[i],DIP_Lake_S[i],J_des_S_S[i],J_ads_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.v_diff_S,TP_Variables.A_Sand_S,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) if TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_S[i],DIP_Lake_S[i],J_des_S_S[i],J_ads_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.v_diff_S,TP_Variables.A_Sand_S,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) > 0 else 0
        DIP_pore_R_S[i+1] = TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_S[i],DIP_Lake_S[i],J_des_R_S[i],J_ads_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.v_diff_R,TP_Variables.A_Rock_S,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) if TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_S[i],DIP_Lake_S[i],J_des_R_S[i],J_ads_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.v_diff_R,TP_Variables.A_Rock_S,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) > 0 else 0
        DIP_pore_P_S[i+1] = TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_S[i],DIP_Lake_S[i],J_des_P_S[i],J_ads_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.v_diff_P,TP_Variables.A_Peat_S,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) if TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_S[i],DIP_Lake_S[i],J_des_P_S[i],J_ads_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.v_diff_P,TP_Variables.A_Peat_S,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) > 0 else 0
            
        Settling_P_N[i] = TP_MBFR.Sett_P(TP_Lake_N[i], DIP_Lake_N[i], Lake_O_A_N[i], Lake_O_Storage_N[i], v_settle_N[i]) if TP_MBFR.Sett_P(TP_Lake_N[i], DIP_Lake_N[i], Lake_O_A_N[i], Lake_O_Storage_N[i], v_settle_N[i]) >0 else 0
        Settling_P_S[i] = TP_MBFR.Sett_P(TP_Lake_S[i], DIP_Lake_S[i], Lake_O_A_S[i], Lake_O_Storage_S[i], v_settle_S[i]) if TP_MBFR.Sett_P(TP_Lake_S[i], DIP_Lake_S[i], Lake_O_A_S[i], Lake_O_Storage_S[i], v_settle_S[i]) >0 else 0
        
        P_diff_M_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_M, DIP_pore_M_N[i], DIP_Lake_N[i], Θ_M, TP_Variables.A_Mud_N,Lake_O_Storage_N[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_M, DIP_pore_M_N[i], DIP_Lake_N[i], Θ_M, TP_Variables.A_Mud_N,Lake_O_Storage_N[i]) >0 else 0
        P_diff_S_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_S, DIP_pore_S_N[i], DIP_Lake_N[i], Θ_S, TP_Variables.A_Sand_N,Lake_O_Storage_N[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_S, DIP_pore_S_N[i], DIP_Lake_N[i], Θ_S, TP_Variables.A_Sand_N,Lake_O_Storage_N[i]) >0 else 0
        P_diff_R_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_R, DIP_pore_R_N[i], DIP_Lake_N[i], Θ_R, TP_Variables.A_Rock_N,Lake_O_Storage_N[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_R, DIP_pore_R_N[i], DIP_Lake_N[i], Θ_R, TP_Variables.A_Rock_N,Lake_O_Storage_N[i]) >0 else 0
        P_diff_P_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_P, DIP_pore_P_N[i], DIP_Lake_N[i], Θ_P, TP_Variables.A_Peat_N,Lake_O_Storage_N[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_P, DIP_pore_P_N[i], DIP_Lake_N[i], Θ_P, TP_Variables.A_Peat_N,Lake_O_Storage_N[i]) >0 else 0
        P_diff_M_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_M, DIP_pore_M_S[i], DIP_Lake_S[i], Θ_M, TP_Variables.A_Mud_S,Lake_O_Storage_S[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_M, DIP_pore_M_S[i], DIP_Lake_S[i], Θ_M, TP_Variables.A_Mud_S,Lake_O_Storage_S[i]) >0 else 0
        P_diff_S_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_S, DIP_pore_S_S[i], DIP_Lake_S[i], Θ_S, TP_Variables.A_Sand_S,Lake_O_Storage_S[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_S, DIP_pore_S_S[i], DIP_Lake_S[i], Θ_S, TP_Variables.A_Sand_S,Lake_O_Storage_S[i]) >0 else 0
        P_diff_R_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_R, DIP_pore_R_S[i], DIP_Lake_S[i], Θ_R, TP_Variables.A_Rock_S,Lake_O_Storage_S[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_R, DIP_pore_R_S[i], DIP_Lake_S[i], Θ_R, TP_Variables.A_Rock_S,Lake_O_Storage_S[i]) >0 else 0
        P_diff_P_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_P, DIP_pore_P_S[i], DIP_Lake_S[i], Θ_P, TP_Variables.A_Peat_S,Lake_O_Storage_S[i]) if TP_MBFR.Diff_P(TP_Variables.v_diff_P, DIP_pore_P_S[i], DIP_Lake_S[i], Θ_P, TP_Variables.A_Peat_S,Lake_O_Storage_S[i]) >0 else 0
                
        
        TP_Lake_N[i+1] = TP_MBFR.TP_Lake_N(L_ext_M[i],Atm_Dep_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_N[i],DIP_pore_S_N[i],DIP_pore_R_N[i],DIP_pore_P_N[i],DIP_Lake_N[i],Q_N2S[i],Lake_O_A_N[i],TP_Lake_N[i],Lake_O_Storage_N[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_N[i]) + ((Sed_Resusp_M_N[i]+Sed_Resusp_S_N[i]+Sed_Resusp_R_N[i]+Sed_Resusp_P_N[i])/Lake_O_Storage_N[i]) if TP_MBFR.TP_Lake_N(L_ext_M[i],Atm_Dep_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_N[i],DIP_pore_S_N[i],DIP_pore_R_N[i],DIP_pore_P_N[i],DIP_Lake_N[i],Q_N2S[i],Lake_O_A_N[i],TP_Lake_N[i],Lake_O_Storage_N[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_N[i])+ ((Sed_Resusp_M_N[i]+Sed_Resusp_S_N[i]+Sed_Resusp_R_N[i]+Sed_Resusp_P_N[i])/Lake_O_Storage_N[i]) > 0 else 0
        TP_Lake_S[i+1] = TP_MBFR.TP_Lake_S(Atm_Dep_S[i],Q_N2S[i],TP_Lake_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_S[i],DIP_pore_S_S[i],DIP_pore_R_S[i],DIP_pore_P_S[i],DIP_Lake_S[i],Q_O_M[i],Lake_O_A_S[i],TP_Lake_S[i],Lake_O_Storage_S[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_S[i]) + ((Sed_Resusp_M_S[i]+Sed_Resusp_S_S[i]+Sed_Resusp_R_S[i]+Sed_Resusp_P_S[i])/Lake_O_Storage_S[i]) if TP_MBFR.TP_Lake_S(Atm_Dep_S[i],Q_N2S[i],TP_Lake_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_S[i],DIP_pore_S_S[i],DIP_pore_R_S[i],DIP_pore_P_S[i],DIP_Lake_S[i],Q_O_M[i],Lake_O_A_S[i],TP_Lake_S[i],Lake_O_Storage_S[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_S[i])+ ((Sed_Resusp_M_S[i]+Sed_Resusp_S_S[i]+Sed_Resusp_R_S[i]+Sed_Resusp_P_S[i])/Lake_O_Storage_S[i]) > 0 else 0
        TP_Lake_Mean[i+1] = ((TP_Lake_N[i+1] + TP_Lake_S[i+1])/2)
        Suspended_Sed_N[i] = ((Sed_Resusp_M_N[i]+Sed_Resusp_S_N[i]+Sed_Resusp_R_N[i]+Sed_Resusp_P_N[i])/Lake_O_Storage_N[i])
        Suspended_Sed_S[i] = ((Sed_Resusp_M_S[i]+Sed_Resusp_S_S[i]+Sed_Resusp_R_S[i]+Sed_Resusp_P_S[i])/Lake_O_Storage_S[i])
        Ext_Load_Rate_N[i] = (L_ext_M[i]+Atm_Dep_N[i])/Lake_O_Storage_N[i]
        Ext_Load_Rate_S[i] = (Atm_Dep_S[i]+Q_N2S[i]*TP_Lake_N[i])/Lake_O_Storage_S[i]
        Load_Out_N[i] = (Q_N2S[i]*TP_Lake_N[i])/Lake_O_Storage_N[i]
        Load_Out_S[i] = (Q_O_M[i]*TP_Lake_S[i])/Lake_O_Storage_S[i]
        P_Load_Cal[i] = M_var.Outlet1USREG[i]*0.028316847*3600*24*TP_Lake_S[i] #mg
        P_Load_StL[i] = M_var.Outlet2USRG[i]*0.028316847*3600*24*TP_Lake_S[i] #mg
        P_Load_South[i] = M_var.TotRegSo[i]*1233.48*TP_Lake_S[i]
          
        
        
        
    Output_df = pd.DataFrame(date_rng_2, columns=['Date']) 

    Output_df['Stage_LO'] =  M_var.Lake_Stage[2:]
    Output_df['S308_Q'] = M_var.Outlet2USRG[2:]
    Output_df['S77_Q'] = M_var.Outlet1USREG[2:]
    Output_df['Storage'] = M_var.Storage[2:]
    Output_df['Cut_back'] = M_var.Cut_back[2:]
    Output_df['NOx'] = NO_MEAN[:-2]
    Output_df['NOx_Load_Cal'] = NO_Load_Cal[:-2]/1E9 #tons
    Output_df['NOx_Load_StL'] = NO_Load_StL[:-2]/1E9 #tons
    Output_df['NOx_Load_South'] = NO_Load_South[:-2]/1E9 #tons
    Output_df['Chla_Lake'] = Sim_Chla[:-2]
    Output_df['Chla_Load_Cal'] = Chla_Load_Cal[:-2]/1E6 #kgs
    Output_df['Chla_Load_StL'] = Chla_Load_StL[:-2]/1E6 #kgs
    Output_df['Chla_Load_South'] = Chla_Load_South[:-2]/1E6 #kgs
  
    Output_df['P_Lake'] = TP_Lake_Mean[:-2]
    Output_df['P_Load_Cal'] = P_Load_Cal[:-2]/1E9 #tons
    Output_df['P_Load_StL'] = P_Load_StL[:-2]/1E9 #tons
    Output_df['P_Load_South'] = P_Load_South[:-2]/1E9 #tons

    return(Output_df)

Exported_File = LOONE_HydNO()
Exported_File.drop(index=Exported_File.index[-1],axis=0,inplace=True)
Exported_File.drop(index=Exported_File.index[-1],axis=0,inplace=True)
Exported_File['Stage_LO'] = Exported_File['Stage_LO'].astype(float)
Exported_File['Storage']=Exported_File['Storage'].astype(float)
Exported_File['S308_Q'] = Exported_File['S308_Q'].astype(float)
Exported_File['S77_Q'] = Exported_File['S77_Q'].astype(float)
Exported_File['Cut_back']=Exported_File['Cut_back'].astype(float)

Exported_File['NOx']=pd.to_numeric(Exported_File['NOx'])
Exported_File['NOx_Load_Cal']=pd.to_numeric(Exported_File['NOx_Load_Cal'])
Exported_File['NOx_Load_StL']=pd.to_numeric(Exported_File['NOx_Load_StL'])
Exported_File['NOx_Load_South']=pd.to_numeric(Exported_File['NOx_Load_South'])

Exported_File['Chla_Lake']=pd.to_numeric(Exported_File['Chla_Lake'])
Exported_File['Chla_Load_Cal']=pd.to_numeric(Exported_File['Chla_Load_Cal'])
Exported_File['Chla_Load_StL']=pd.to_numeric(Exported_File['Chla_Load_StL'])
Exported_File['Chla_Load_South']=pd.to_numeric(Exported_File['Chla_Load_South'])

Exported_File['P_Lake']=pd.to_numeric(Exported_File['P_Lake'])
Exported_File['P_Load_Cal']=pd.to_numeric(Exported_File['P_Load_Cal'])
Exported_File['P_Load_StL']=pd.to_numeric(Exported_File['P_Load_StL'])
Exported_File['P_Load_South']=pd.to_numeric(Exported_File['P_Load_South'])

Exported_File = Exported_File.set_index('Date')
Exported_File.index = pd.to_datetime(Exported_File.index, unit = 'ns')
Exported_File_Mean = Exported_File.resample('M').mean()
Exported_File_Sum = Exported_File.resample('M').sum()

# Exported_File_Mean.to_csv('./Outputs/Exported_File_Chla_Opt_Mean_%s.csv'%Pre_defined_Variables.Schedule)
# Exported_File_Sum.to_csv('./Outputs/Exported_File_Chla_Opt_Sum_%s.csv'%Pre_defined_Variables.Schedule)

#Temp
Exported_File_Mean.to_csv('./Outputs/Exported_File_NOx_Chla_TP_Opt_Mean_%s_Scen1.csv'%Pre_defined_Variables.Schedule)
Exported_File_Sum.to_csv('./Outputs/Exported_File_NOx_Chla_TP_Opt_Sum_%s_Scen1.csv'%Pre_defined_Variables.Schedule)

