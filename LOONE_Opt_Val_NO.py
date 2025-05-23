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
    Results_data = pd.read_csv('./Outputs/Dec_Var.csv')
    Results = Results_data['Value']

    NO_1 = Results[0]
    NO_2 = Results[1]
    S77_DV = Results[2:14]
    S77_DV = S77_DV.reset_index(drop=True)
    S308_DV = Results[14:26]
    S308_DV = S308_DV.reset_index(drop=True)

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
    Q_in = pd.read_csv('./Model_Data_Filled_20082023/LO_Inflows_20082023.csv')
    NOx_In = pd.read_csv('./Model_Data_Filled_20082023/Daily_NOx_External_Loads_20082023.csv') #mg
    S65E_NO_data = pd.read_csv('./Model_Data_Filled_20082023/S65E_NO_Interpolated.csv') #mg/m3
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
    Temp_data = pd.read_csv('./Model_Data_Filled_20082023/Temp_Avg_20082023.csv')
    DO_data = pd.read_csv('./Model_Data_Filled_20082023/LO_Avg_DO_20082023.csv')
    RAD_data = pd.read_csv('./Model_Data_Filled_20082023/LO_RADT_20082023.csv')
    Storage = pd.read_csv('./Model_Data_Filled_20082023/Average_LO_Storage_20082023.csv')
    # Storage_deviation = pd.read_csv('./Model_Data_Filled_20082023/Storage_Dev_20082023.csv')
    Chla_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_Merged_Chla.csv') #microgram/L
    Chla_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_Merged_Chla.csv') #microgram/L
    Photoperiod = pd.read_csv('./Model_Data_Filled_20082023/PhotoPeriod_20082023.csv') 
    Photoperiod['Date'] = pd.to_datetime(Photoperiod['Date'])
    LO_DIP_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_OP.csv') #mg/m3
    LO_DIP_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_OP.csv') #mg/m3
    LO_DIN_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_DIN.csv') #mg/m3
    LO_DIN_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_DIN.csv') #mg/m3
    # RAD = RAD_data['RAD'].astype(float) * 0.484583*1000 #This one gives the best results for NO simulations
    RAD = RAD_data['RAD'].astype(float) * 4.6*1000 

    Temp = Temp_data['Mean_T'].astype(float)
    DO = DO_data['DO'].astype(float)
    External_NO = NOx_In['NOx_Load_mg'].astype(float) #mg
    S65E_NO = S65E_NO_data['NO'].astype(float) #mg/m3
    # N-S Procedure
    N_Per = 0.43
    S_Per = 0.57
    
    volume  = Storage['Storage_m3'].astype(float) #m3
    volume_N = volume*N_Per
    volume_S = volume*S_Per
    # Storage_dev = Storage_deviation['DS_dev'].astype(float) #acft
    Stage = Storage['Stage_ft'].astype(float) * 0.3048 #m
    Q_I = Q_in['Inflows_cmd'].astype(float) #m3

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
    Cal_Par = Cal_Par_Opt['Par']
    
       
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
    kni0 = Cal_Par[0]
    kni = Cal_Par[1]
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
    kden0 =  Cal_Par[2]
    kden =  Cal_Par[3]
    Theta_deni = 1.06
    
    #Light Attenuation due to water Kw and algae Kc
    Kw = Cal_Par[4] #1.7
    Kc = Cal_Par[5] #0.029
    #Light limitation and inhibition coefficients K1 and K2
    K1 = Cal_Par[6] #100
    K2 = Cal_Par[7] #700
    # Maximum Phytoplankton growth rate
    G_max = Cal_Par[8] #1.5
    
    #Half Saturation Coefficients
    K_NH = Cal_Par[9] #0.1
    K_Nitr = Cal_Par[10]
    K_DIN = K_NH + K_Nitr
    KP = Cal_Par[11]
    K_TN = Cal_Par[12] #0.1
    YNOChla = Cal_Par[13] #0.1
    
    #Temperatures
    T_opt = Cal_Par[14] #15 #C
    T_min = Cal_Par[15] #10 #C
    T_max = Cal_Par[16] #30 #C
    # Dissolved Oxygen    
    KDO = Cal_Par[17]
    #Sediment Release
    S_NO = Cal_Par[18] #mg/m2/d
    
    Theta_NO = 1.06
    NO_Temp_Adj = Theta_NO**(Temp-20)  
    Area = 1730 *1E6 #m2
    Area_N = Area * N_Per
    Area_S = Area * S_Per
    
    X = len(Q_in.index)
    
    #Nitrogen and Phosphorus Limiting
    f_P = DIP_Obs/(KP+DIP_Obs)
    f_N = DIN_Obs/(K_DIN+DIN_Obs)
    # f_N_N = np.zeros(X,dtype = object)
    # f_N_S = np.zeros(X,dtype = object)
    # f_N = np.zeros(X,dtype = object)
    
    f_P_N = DIP_N_Obs/(KP+DIP_N_Obs)
    f_P_S = DIP_S_Obs/(KP+DIP_S_Obs)
    
    f_N_N = DIN_N_Obs/(K_DIN+DIN_N_Obs)
    f_N_S = DIN_S_Obs/(K_DIN+DIN_S_Obs)
    
    Q_O = np.zeros(n_rows,dtype = object)

    Q_I_M = np.zeros(n_rows,dtype = object)
    Q_O_M = np.zeros(n_rows,dtype = object)
    Q_N2S = np.zeros(n_rows,dtype = object)
    External_NO_M = np.zeros(n_rows,dtype = object)

    Nit_R_N = np.zeros(n_rows,dtype = object)
    Nit_R_N_T = np.zeros(n_rows,dtype = object)
    Nit_R_S = np.zeros(n_rows,dtype = object)
    Nit_R_S_T = np.zeros(n_rows,dtype = object)
    Denit_R_N = np.zeros(n_rows,dtype = object)
    Denit_R_N_T = np.zeros(n_rows,dtype = object)
    Denit_R_S = np.zeros(n_rows,dtype = object)
    Denit_R_S_T = np.zeros(n_rows,dtype = object)
    
    fT = np.zeros(n_rows,dtype = object)
    fL = np.zeros(n_rows,dtype = object)
    NO_N = np.zeros(n_rows,dtype = object)
    NO_S = np.zeros(n_rows,dtype = object)
    NO_MEAN = np.zeros(n_rows,dtype = object)
    
    NO_Load_Cal = np.zeros(n_rows,dtype = object)
    NO_Load_StL = np.zeros(n_rows,dtype = object)
    NO_Load_South = np.zeros(n_rows,dtype = object)

    
    #Water depth could be stage - B.L. (variable) or a constant average depth (e.g., 2.5 m)

    z = np.zeros(n_rows,dtype = object)
    # B_L = 1 #m
    # z = Stage - B_L #m

    for i in range(X):
    #     # z[i] = 2.1296864
        z[i] = 2.5
    
    fox_min = Cal_Par[19]
    DO_Cr_n = Cal_Par[20]
    DO_Opt_n = Cal_Par[21]
    a = Cal_Par[22]
    DO_Cr_d = Cal_Par[23]
    DO_Opt_d = Cal_Par[24]
    
    
    Grazing = 0
    Theta_G = 1.06
    Grazing_T = Grazing*Theta_G**(Temp-20)
    NO_N[0] = NOx_N_Obs[0]
    NO_S[0] = NOx_S_Obs[0]
    
    NO_MEAN[0] = (NO_N[0]+NO_S[0])*0.5
    
    # Map the monthly values to the DataFrame
    Nitro_Model_Output = pd.DataFrame(Q_in['date'],columns=['date'])
    Nitro_Model_Output['date'] = pd.to_datetime(Nitro_Model_Output['date'])
    print("LOONE Nitrogen Module is Running!")
    #####################################################################################################################
    M_var.Zone_Code[0] = LO_FNs.Zone_Code(M_var.Lake_Stage[0],df_WSMs['A'].iloc[0],df_WSMs['B'].iloc[0],df_WSMs['C'].iloc[0],df_WSMs['D3'].iloc[0],df_WSMs['D2'].iloc[0],df_WSMs['D1'].iloc[0],df_WSMs['D0'].iloc[0],df_WSMs['WSM1'].iloc[0])
    M_var.LO_Zone[0] = LO_FNs.LO_Zone(M_var.Zone_Code[0])
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
            elif (NO_S[i] <= NO_1) and (date_rng_6[i+2].month in [1,2,3,4,11,12]):
                M_var.Outlet2USRG[i+2] = S308_DV[(date_rng_6[i+2].month)-1]
            elif (NO_S[i] <= NO_2) and (date_rng_6[i+2].month in [5,6,7,8,9,10]):
                M_var.Outlet2USRG[i+2] = S308_DV[(date_rng_6[i+2].month)-1]
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
            elif (NO_S[i] <= NO_1) and (date_rng_6[i+2].month in [1,2,3,4,11,12]):
                M_var.Outlet1USREG[i+2] = S77_DV[(date_rng_6[i+2].month)-1]
            elif (NO_S[i] <= NO_2) and (date_rng_6[i+2].month in [5,6,7,8,9,10]):
                M_var.Outlet1USREG[i+2] = S77_DV[(date_rng_6[i+2].month)-1]
            else:
                M_var.Outlet1USREG[i+2] = 0
        M_var.Outlet1DS[i+2] = LO_FNs.Outlet1DS(M_var.Outlet1USREG[i+2],M_var.Outlet1USEWS[i+2],M_var.ROwest[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S79'],Pre_defined_Variables.Outlet1DS_Switch)
        M_var.TotRegEW[i+2] = (M_var.Outlet1USREG[i+2] + M_var.Outlet2USRG[i+2])*1.9835
        M_var.Choose_WCA[i+2] = LO_FNs.Choose_WCA(Data.SFWMM_Daily_Outputs.at[i+2, 'RegWCA'],Pre_defined_Variables.Option_RegWCA,Pre_defined_Variables.Constant_RegWCA)
        M_var.RegWCA[i+2] = min(Pre_defined_Variables.MaxCap_RegWCA , Pre_defined_Variables.Multiplier_RegWCA*M_var.Choose_WCA[i+2])
        M_var.Choose_L8C51[i+2] = LO_FNs.Choose_L8C51(Data.SFWMM_Daily_Outputs.at[i+2, 'RegL8C51'],Pre_defined_Variables.Option_RegL8C51,Pre_defined_Variables.Constant_RegL8C51)
        M_var.RegL8C51[i+2] = min(Pre_defined_Variables.MaxCap_RegL8C51 , Pre_defined_Variables.Multiplier_RegL8C51*M_var.Choose_L8C51[i+2])
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
    
        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48 #m3/d
            Q_I_M[i] = Q_I[i]
            External_NO_M[i] = External_NO[i]         
    
        Q_N2S[i] = (Q_I_M[i]*1 + Q_O_M[i]*0)   
    
        t = np.linspace(1,2,num = 2)
        Nit_R_N[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0,kni,NH4_N_Obs[i],K_NH,DO[i],KDO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_N_T[i] = Nit_R_N[i]*Theta_ni**(Temp[i]-20) 
        Nit_R_S[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0,kni,NH4_S_Obs[i],K_NH,DO[i],KDO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_S_T[i] = Nit_R_S[i]*Theta_ni**(Temp[i]-20) 
    
        Denit_R_N[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0,kden,NO_N[i],K_Nitr,DO[i],KDO,DO_Cr_d,DO_Opt_d)
        Denit_R_N_T[i] = Denit_R_N[i]*Theta_deni**(Temp[i]-20) 
        Denit_R_S[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0,kden,NO_S[i],K_Nitr,DO[i],KDO,DO_Cr_d,DO_Opt_d)
        Denit_R_S_T[i] = Denit_R_S[i]*Theta_deni**(Temp[i]-20) 
        # f_N_N[i] = (NO_N[i]+ NH4_N_Obs[i])/(NO_N[i]+ NH4_N_Obs[i]+K_DIN)
        # f_N_S[i] = (NO_S[i]+ NH4_S_Obs[i])/(NO_S[i]+ NH4_S_Obs[i]+K_DIN)
        # f_N[i] = (NO_MEAN[i]+ NH4_Obs[i])/(NO_MEAN[i]+ NH4_Obs[i]+K_DIN)
        # f_N_N[i] = (NO_N[i]+ NH4_N[i])/(NO_N[i]+ NH4_N[i]+K_DIN)
        # f_N_S[i] = (NO_S[i]+ NH4_S[i])/(NO_S[i]+ NH4_S[i]+K_DIN)
        # f_N[i] = (NO_MEAN[i]+ NH4_Obs[i])/(NO_MEAN[i]+ NH4_Obs[i]+K_DIN)
    
    
        fT[i] = f_T_alt1(Temp[i],T_opt,T_min,T_max)
        # fL[i] = f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw,Kc,Sim_Chla[i],z[i],K1,K2)*1.5 if Nitro_Model_Output['date'].iloc[i].month in (6,7,8,9,10) else f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw,Kc,Sim_Chla[i],z[i],K1,K2)*0.75
        fL[i] = f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw,Kc,Chla[i],z[i],K1,K2)

        DfEq_Res_N = odeint(NOx_N_DiffEq,NO_N[i],t,args=(External_NO[i],Atm_Deposition_N,Q_N2S[i],Nit_R_N_T[i],Denit_R_N_T[i],NH4_N_Obs[i],volume_N[i],G_max,fT[i],fL[i],f_P_N[i],f_N_N[i],K_NH,K_TN,YNOChla,Chla_N[i],S_NO,Area_N,NO_Temp_Adj[i],))
        NO_N[i+1] = DfEq_Res_N[:,0][1]
        DfEq_Res_S = odeint(NOx_S_DiffEq,NO_S[i],t,args=(Atm_Deposition_S,Q_N2S[i],Q_O_M[i],NO_N[i],Nit_R_S_T[i],Denit_R_S_T[i],NH4_S_Obs[i],volume_S[i],G_max,fT[i],fL[i],f_P_S[i],f_N_S[i],K_NH,K_TN,YNOChla,Chla_S[i],S_NO,Area_S,NO_Temp_Adj[i],))
        NO_S[i+1] = DfEq_Res_S[:,0][1]
        NO_MEAN[i+1] = (NO_N[i+1] + NO_S[i+1])*0.5
          
        
        NO_Load_Cal[i] = M_var.Outlet1USREG[i]*0.028316847*3600*24*NO_S[i] #mg/d P
        NO_Load_StL[i] = M_var.Outlet2USRG[i]*0.028316847*3600*24*NO_S[i] #mg/d P
        NO_Load_South[i] = M_var.TotRegSo[i]*1233.48 *NO_S[i] #mg/d P

        Output_df = pd.DataFrame(date_rng_2, columns=['Date']) #1/1/2008-12/31/2018

        Output_df['Stage_LO'] =  M_var.Lake_Stage[2:]
        Output_df['S308_Q'] = M_var.Outlet2USRG[2:]
        Output_df['S77_Q'] = M_var.Outlet1USREG[2:]
        Output_df['Storage'] = M_var.Storage[2:]
        Output_df['Cut_back'] = M_var.Cut_back[2:]
        Output_df['NO_Lake'] = NO_MEAN[:-2]
        Output_df['NO_Load_Cal'] = NO_Load_Cal[:-2]/1E9 #tons
        Output_df['NO_Load_StL'] = NO_Load_StL[:-2]/1E9 #tons
        Output_df['NO_Load_South'] = NO_Load_South[:-2]/1E9 #tons
  
    return(Output_df)

Exported_File = LOONE_HydNO()
Exported_File.drop(index=Exported_File.index[-1],axis=0,inplace=True)
Exported_File.drop(index=Exported_File.index[-1],axis=0,inplace=True)
Exported_File['Stage_LO'] = Exported_File['Stage_LO'].astype(float)
Exported_File['Storage']=Exported_File['Storage'].astype(float)
Exported_File['S308_Q'] = Exported_File['S308_Q'].astype(float)
Exported_File['S77_Q'] = Exported_File['S77_Q'].astype(float)
Exported_File['Cut_back']=Exported_File['Cut_back'].astype(float)
Exported_File['NO_Lake']=pd.to_numeric(Exported_File['NO_Lake'])
Exported_File['NO_Load_Cal']=pd.to_numeric(Exported_File['NO_Load_Cal'])
Exported_File['NO_Load_StL']=pd.to_numeric(Exported_File['NO_Load_StL'])
Exported_File['NO_Load_South']=pd.to_numeric(Exported_File['NO_Load_South'])

Exported_File = Exported_File.set_index('Date')
Exported_File.index = pd.to_datetime(Exported_File.index, unit = 'ns')
Exported_File_Mean = Exported_File.resample('M').mean()
Exported_File_Sum = Exported_File.resample('M').sum()

Exported_File_Mean.to_csv('./Outputs/Exported_File_NO_Opt_Mean_%s.csv'%Pre_defined_Variables.Schedule)
Exported_File_Sum.to_csv('./Outputs/Exported_File_NO_Opt_Sum_%s.csv'%Pre_defined_Variables.Schedule)

