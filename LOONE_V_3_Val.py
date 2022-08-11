# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 18:44:37 2021

@author: osama
"""

#This Script incorporates the Comprehensive LOONE Model!
Working_Path = 'C:/Osama_PC/LOONE/Model/LOONE_Model'
import os
import pandas as pd
from datetime import datetime
import numpy as np
from scipy import interpolate
from calendar import monthrange  
os.chdir('%s'%Working_Path) 
from Pre_defined_Variables_July4_22 import Pre_defined_Variables 
from LO_FNs_July4_22 import LO_FNs
from Stg_Sto_Ar import Stg_Sto_Ar 
from LONINO_FNs_July4_22 import LONINO_FNs
from Dec_Tree_FNs import Dec_Tree_FNs
from WCA_Stages_Cls import WCA_Stages_Cls
from Additional_Fncs import Add_Fn
from THC_Class import THC_Class
from Data import Data
from TP_Variables_Regions import TP_Variables
import TP_Mass_Balance_Functions_Regions as TP_MBFR

def LOONE_HydNut(): 
    # Based on the defined Start and End year, month, and day on the Pre_defined_Variables File, Startdate and enddate are defined. 
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    startdate = datetime(year, month, day).date() 
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    begdateCS = datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.enddate_entry)
    enddate = datetime(year, month, day).date()
        
#############################################################################################################################      
    Results_data = pd.read_csv('./Outputs/Opt_Decision_Var.csv')
    Results = Results_data['Value']
    P_1 = Results[0]
    P_2 = Results[1] 
    S77_DV = Results[2:14]    
    S308_DV = Results[14:26]
    #First, I interpolated each Water Shortage Management (WSMs) and each Regulation Schedule Breakpoint Zone (D, C, B, and A).
    #Set time frame for model run such that it starts on the defined startdate but ends on 1/1/(endyear+1)
    date_rng_1 = pd.date_range(start = startdate, end = '1/1/%d'%(Pre_defined_Variables.endyear+1), freq= 'D')
    #Create a data frame with a date column
    df_WSMs = pd.DataFrame(date_rng_1, columns =['date'])
    WSM_length = len(df_WSMs.index)
    WSM_Count = np.zeros(WSM_length)
    Oper_Zones = list(Data.WSMs_RSBKs)
    Oper_Zones.remove('Date')
    Oper_Zones.remove('Day')
    for i in (Oper_Zones):
        globals()[i] = np.zeros(WSM_length) 
    for i in range(WSM_length):
        if df_WSMs['date'].iloc[i].month == 1 and df_WSMs['date'].iloc[i].day == 1:
            J = 1
        elif df_WSMs['date'].iloc[i].month == 2 and df_WSMs['date'].iloc[i].day == 29:
            J = J
        else: 
            J = J+1
        WSM_Count[i] = J
        for j in Oper_Zones:
            globals()[j][i] = interpolate.interp1d(Data.WSMs_RSBKs['Day'], Data.WSMs_RSBKs['%s'%j], kind = 'linear')(WSM_Count[i])
            
    df_WSMs['count'] = WSM_Count
    for j in Oper_Zones:
        df_WSMs['%s'%j] = globals()[j]
    if Pre_defined_Variables.Opt_NewTree == 1:
        df_WSMs['C-b'] = Data.WSMs_RSBKs['C-b_NewTree']
    else:
        df_WSMs['C-b'] = Data.WSMs_RSBKs['C-b_NoNewTree']
    df_WSMs.drop(['C-b_NewTree','C-b_NoNewTree'],axis=1, inplace=True)
# Add one row in top of the dataframe (i.e. Dec/31/previous year) where the values = values of the original first row in the dataframe!
    First_row = pd.DataFrame({'WSM4':df_WSMs['WSM4'].iloc[0], 'WSM3':df_WSMs['WSM3'].iloc[0], 'WSM2':df_WSMs['WSM2'].iloc[0], 'WSM1':df_WSMs['WSM1'].iloc[0], 'D0':df_WSMs['D0'].iloc[0], 'D1':df_WSMs['D1'].iloc[0], 'D2':df_WSMs['D2'].iloc[0], 'D3':df_WSMs['D3'].iloc[0], 'C':df_WSMs['C'].iloc[0], 'B':df_WSMs['B'].iloc[0], 'A':df_WSMs['A'].iloc[0], 'C-b':df_WSMs['C-b'].iloc[0]}, index = [0])
    df_WSMs = pd.concat([First_row, df_WSMs]).reset_index(drop= True)
#############################################################################
     #The Following Code interpolates daily LOSA demand from weekly data for 6 differnet datasets where the user defines the LOSA demand that will be used based on a Code (1:6).
    #Set time frame for model run
    date_rng_2 = pd.date_range(start=startdate, end = enddate, freq= 'D')
    #Create a data frame with a date column
    LOSA_dmd = pd.DataFrame(date_rng_2, columns =['date'])
    
    N = []
    Wk = []
    #Generate a count list
    for i in LOSA_dmd['date']:
        if i.month == 1 and i.day == 1:
            n = 0
        else:
            n = n + 1
        N.append(n)
    LOSA_dmd['count'] = N
    #Calculate the week number for all rows in the data frame
    for i in LOSA_dmd['count']:
        if i > 363:
            J = 52
        else:
            J = int(i/7)+1
        Wk.append(J)
    LOSA_dmd['Week_num'] = Wk
    dd = []
    #Calculate daily water demand
    for i in LOSA_dmd['Week_num']:
        D = ((Data.Weekly_dmd['C%s'%Pre_defined_Variables.Code].iloc[i-1])/7)*(Pre_defined_Variables.Multiplier/100)
        dd.append(D)
    LOSA_dmd['Daily_demand'] = dd
    ##############################################################################################
    # I determine daily values for the Tributary conditions and Seasonal and Multi-Seasonal LONINO classes 
    #using a weekly Trib. Condition data and Monthly LONINO data. 
    
    #Generate weekly time step date column where frequency is 'W-Fri' to start on 01/01/2008.
    #FIXME: Always check here for start date, end date, and frequency to match with the Trib. Condition weekly data obtained.
    year, month, day = map(int, Pre_defined_Variables.enddate_TC)
    enddate_TC = datetime(year, month, day).date() 
    date_rng_3 = pd.date_range(start=startdate, end = enddate_TC, freq= 'W-Fri')
    
    #Generate the Tributary Condition Dataframe.
    Trib_Cond_df = pd.DataFrame(date_rng_3, columns =['date'])
    TC_Count = len(Trib_Cond_df.index)
    #Identify Weekly Rainfall (RF) classes based on the RF values.
    RF_Cls = np.zeros(TC_Count)
    S65E_Cls = np.zeros(TC_Count)
    Palmer_Cls = np.zeros(TC_Count)
    NetInflow_Cls = np.zeros(TC_Count)
    Max_RF_S65E = np.zeros(TC_Count)
    Max_Palmer_NetInf = np.zeros(TC_Count)
    for i in range(TC_Count):
    #Identify Weekly RF classes based on the RF values.
        RF_Cls[i] = LONINO_FNs.RF_Cls(Data.Wkly_Trib_Cond['NetRF'].iloc[i])
    #Identify Weekly S65E classes based on the S65E values.
        S65E_Cls[i] = LONINO_FNs.S65E_Cls(Data.Wkly_Trib_Cond['S65E'].iloc[i])
    #Identify Weekly Palmer Index classes based on the Palmer Index values.
        Palmer_Cls[i] = LONINO_FNs.Palmer_Cls(Data.Wkly_Trib_Cond['Palmer'].iloc[i])
    #Identify Weekly NetInflow classes based on the Netinflow values.
        NetInflow_Cls[i] = LONINO_FNs.NetInflow_Cls(Data.Wkly_Trib_Cond['NetInf'].iloc[i])
    #Calculate maximum value between RF and S65E
        Max_RF_S65E[i] =  max(RF_Cls[i], S65E_Cls[i])
    #Calculate Maximum value between Palmer and NetInflow 
        Max_Palmer_NetInf[i] =  max(Palmer_Cls[i], NetInflow_Cls[i])
    #Determine Tributary Condition Index either the maximum of Palmer and Netinf or the maximum of RF and S65E based on the user defined TCI.
    if Pre_defined_Variables.TCI == 1:
        Trib_Cond_df['TCI'] = Max_Palmer_NetInf
    else:
        Trib_Cond_df['TCI'] = Max_RF_S65E
    #Generate a monthly time step date column
    date_rng_4 = pd.date_range(start=startdate, end = enddate, freq= 'MS')
    #Create a LONINO Dataframe
    LONINO_df = pd.DataFrame(date_rng_4, columns =['date'])
    LONINO_Count = len(LONINO_df.index)
    Seas = np.zeros(LONINO_Count)
    M_Seas = np.zeros(LONINO_Count)
    for i in range(LONINO_Count):
        Seas[i] = Data.LONINO_Seas_data['%s'%LONINO_df['date'].iloc[i].month].iloc[LONINO_df['date'].iloc[i].year-Pre_defined_Variables.startyear]
        M_Seas[i] = Data.LONINO_Mult_Seas_data['%s'%LONINO_df['date'].iloc[i].month].iloc[LONINO_df['date'].iloc[i].year-Pre_defined_Variables.startyear]
    LONINO_df['LONINO_Seas'] = Seas
    LONINO_df['LONINO_Mult_Seas'] = M_Seas
    LONINO_Seas_cls = np.zeros(Pre_defined_Variables.Month_N)
    LONINO_M_Seas_cls = np.zeros(Pre_defined_Variables.Month_N)
    for i in range(Pre_defined_Variables.Month_N):
        LONINO_Seas_cls[i] = LONINO_FNs.LONINO_Seas_cls(LONINO_df['LONINO_Seas'].iloc[i])
        LONINO_M_Seas_cls[i] = LONINO_FNs.LONINO_M_Seas_cls(LONINO_df['LONINO_Mult_Seas'].iloc[i])
    LONINO_df['LONINO_Seasonal_Cls'] = LONINO_Seas_cls
    LONINO_df['LONINO_Mult_Seasonal_Cls'] = LONINO_M_Seas_cls
    #Genarate a daily date range 5 
    date_rng_5 = pd.date_range(start = startdate, end = enddate, freq ='D')
    TC_LONINO_df = pd.DataFrame(date_rng_5, columns = ['Date'])
    row_nm = len(TC_LONINO_df.index)        
    Trib_Cond = np.zeros(row_nm)
    for i in range(row_nm):
        Trib_Cond[i] = Trib_Cond_df['TCI'].iloc[int(i/7)]
    TC_LONINO_df['Tributary_Condition'] = Trib_Cond
    data_S_MS = [LONINO_df['date'], LONINO_df['LONINO_Seasonal_Cls'], LONINO_df['LONINO_Mult_Seasonal_Cls']]
    headers_S_MS = ['date', 'LONINO_Seasonal_Class', 'LONINO_MSeasonal_Class']
    LONINO_Seas_MSeas_df = pd.concat(data_S_MS, axis=1, keys=headers_S_MS)
    LONINO_Seas_MSeas_df = LONINO_Seas_MSeas_df.set_index('date')
    LONINO_Seas_MSeas_df = LONINO_Seas_MSeas_df.reindex(date_rng_5, method = 'ffill')
    LONINO_Seas_MSeas_df = LONINO_Seas_MSeas_df.reset_index()
    TC_LONINO_df['LONINO_Seasonal_Classes']= LONINO_Seas_MSeas_df['LONINO_Seasonal_Class']
    TC_LONINO_df['LONINO_MultiSeasonal_Classes']= LONINO_Seas_MSeas_df['LONINO_MSeasonal_Class']
    #############################################################################################
    WCA_Stages_df = WCA_Stages_Cls(TC_LONINO_df)
    #A dataframe to determine eachday's season (Months 11,12,1,2 are Season 1, Months 3,4,5 are season 2, Months 6,7 are season 3, Months 8,9,10 are season 4 )
    Seasons = pd.DataFrame(date_rng_5, columns =['date'])
    Seas_Count = len(Seasons.index)
    Daily_Seasons = np.zeros(Seas_Count)
    Mon = np.zeros(Seas_Count)
    for i in range(Seas_Count):
        if Seasons['date'].iloc[i].month > 2 and Seasons['date'].iloc[i].month < 6:
            S = 2
        elif Seasons['date'].iloc[i].month > 5 and Seasons['date'].iloc[i].month < 8:
            S = 3
        elif Seasons['date'].iloc[i].month > 7 and Seasons['date'].iloc[i].month < 11:
            S = 4
        else:
            S = 1
        Daily_Seasons[i] = S
        Mon[i] = Seasons['date'].iloc[i].month
    Seasons['Season'] = Daily_Seasons 
    Seasons['Month'] = Mon
##################################################################################################################   
    #This following Script runs the main model daily simulations.
    #Set time frame for model run
    date_rng_6 = pd.date_range(start='12/30/%d'%(Pre_defined_Variables.startyear-1), end = enddate, freq= 'D')
    #Create a data frame with a date column
    LO_Model = pd.DataFrame(date_rng_6, columns =['date'])
    #Create Net_Inflow Column in (ac-ft)
    LO_Model['Net_Inflow'] = Data.NetInf_Input['Netflows_acft']
    #Determine number of time steps
    n_rows = len(LO_Model.index)
    #
    LO_Model['LOSA_dmd_SFWMM'] = Data.SFWMM_W_dmd['LOSA_dmd'] * (Pre_defined_Variables.Mult_LOSA/100)
    #
    LO_Model['C44RO'] = Data.C44_Runoff['C44RO']
    ##################################
    #This following section calculates the parameter states (trib conds, stage tests, seasonal & multi-seasonal LONINO) and sets a 4-digit code for the branch.  The branch code is used with the Routing sheet to determine release rates.
    DecTree_df = pd.DataFrame(date_rng_5, columns = ['Date'])
    Counter = len(DecTree_df)
    DecTree_df['Zone_B_MetFcast'] = TC_LONINO_df['LONINO_Seasonal_Classes']
#######################################################################################################################################################################
    #Create a dataframe that includes Monthly Mean Basin Runoff & BaseFlow-Runoff & Runoff-Baseflow (cfs)
    date_rng_11 = pd.date_range(start=startdate, end = enddate, freq= 'MS')
    date_rng_11d = pd.date_range(start=startdate, end = enddate, freq= 'D')
    date_rng_11d.name = 'Date'
    Basin_RO = pd.DataFrame(date_rng_11, columns =['date'])
    
    #Define Caloosahatchee and St. Lucie Baseflows.
    Cal_Est_baseflow = Data.S77_RegRelRates['Zone_D0'].iloc[0]   
    StlEst_baseflow = Data.S80_RegRelRates['Zone_D0'].iloc[0]
    #Calculta number of months in the timeseries data.
    num_B_R = len(Basin_RO.index)
    BS_C43RO = np.zeros(num_B_R)
    BS_C44RO = np.zeros(num_B_R)
    C44RO_SLTRIB = np.zeros(num_B_R)
    C44RO_BS = np.zeros(num_B_R)
    Num_days = np.zeros(num_B_R)
    for i in range(num_B_R) :
        Num_days[i] = monthrange(Basin_RO['date'].iloc[i].year, Basin_RO['date'].iloc[i].month)[1] #no. of days in each time step month.
        BS_C43RO[i] = max(0, (Cal_Est_baseflow - Data.C43RO['C43RO'].iloc[i]))
        BS_C44RO[i] = max(0, (StlEst_baseflow - Data.C44RO['C44RO'].iloc[i]))
        C44RO_SLTRIB[i] = BS_C44RO[i] + Data.SLTRIB['SLTRIB_cfs'].iloc[i]
        C44RO_BS[i] = max(0, Data.C44RO['C44RO'].iloc[i] - StlEst_baseflow)*Num_days[i]
    Basin_RO['Ndays'] = Num_days
    Basin_RO['C43RO'] = Data.C43RO['C43RO']
    Basin_RO['BS-C43RO'] = BS_C43RO
    Basin_RO['C44RO'] = Data.C44RO['C44RO']
    Basin_RO['BS-C44RO'] = BS_C44RO
    Basin_RO['SLTRIB'] = Data.SLTRIB['SLTRIB_cfs']
    Basin_RO['C44RO_SLTRIB'] = C44RO_SLTRIB
    Basin_RO['C44RO-BS'] = C44RO_BS
    #TODO#I stopped here for MANUAL
    LO_Model['C43RO'] = Data.C43RO_Daily['C43RO']
    #date_year = pd.date_range(start='1/1/2008', end = '12/1/2020', freq= 'YS')
    #Alternating Lake Okeechobee Damaging Discharges
    #ALODQ = pd.DataFrame(date_year, columns =['date'])
    #num_ALODQ = len(ALODQ.index)
    #for i in range(1,num_ALODQ):
    S80avgL1 = Data.Pulses['S-80_L1'].mean()
    S80avgL2 = Data.Pulses['S-80_L2'].mean()
    S80avgL3 = Data.Pulses['S-80_L3'].mean()
    S77avgL1 = Data.Pulses['S-77_L1'].mean() #LORS
    S77avgL2 = Data.Pulses['S-77_L2'].mean() #LORS
    #S77avgL1 = Data.Pulses['S-77_L1_WSE_Run25'].mean() #RUN25 or WSE
    #S77avgL2 = Data.Pulses['S-77_L2_WSE_Run25'].mean() #RUN25 or WSE
    S77avgL3 = Data.Pulses['S-77_L3'].mean()
    #Code the VLOOKUP(date, MonthlyBaseflowTable,5) for S80BS
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
    Count_AP = len(AdapProt_df.index) 
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
    #Generate a dataframe with same timeframe as AdapProt_df for each percentage (i.e. 10,20,25,30,40,45,50,60)
    #I assigned the value of each date to this date for each year of the time series (i.e. 2008-2020)
    Targ_Stg_df = pd.DataFrame(date_rng_5, columns = ['dates'])
    V10per = np.zeros(len(Targ_Stg_df))
    V20per = np.zeros(len(Targ_Stg_df))
    V25per = np.zeros(len(Targ_Stg_df))
    V30per = np.zeros(len(Targ_Stg_df))
    V40per = np.zeros(len(Targ_Stg_df))
    V45per = np.zeros(len(Targ_Stg_df))
    V50per = np.zeros(len(Targ_Stg_df))
    V60per = np.zeros(len(Targ_Stg_df))

    for i in range(len(Targ_Stg_df)):
        V10per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,10,Targ_Stg)   
        V20per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,20,Targ_Stg)
        V25per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,25,Targ_Stg)
        V30per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,30,Targ_Stg)   
        V40per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,40,Targ_Stg) 
        V45per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,45,Targ_Stg)
        V50per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,50,Targ_Stg)   
        V60per[i] = Add_Fn.Replicate(Targ_Stg_df['dates'].iloc[i].year, Targ_Stg_df['dates'].iloc[i].timetuple().tm_yday,60,Targ_Stg)   

    V10per_c = [x for x in V10per if ~np.isnan(x)]
    V20per_c = [x for x in V20per if ~np.isnan(x)]
    V25per_c = [x for x in V25per if ~np.isnan(x)]
    V30per_c = [x for x in V30per if ~np.isnan(x)]
    V40per_c = [x for x in V40per if ~np.isnan(x)]
    V45per_c = [x for x in V45per if ~np.isnan(x)]
    V50per_c = [x for x in V50per if ~np.isnan(x)]
    V60per_c = [x for x in V60per if ~np.isnan(x)]
    Targ_Stg_df['10%'] = V10per_c
    Targ_Stg_df['20%'] = V20per_c
    Targ_Stg_df['25%'] = V25per_c
    Targ_Stg_df['30%'] = V30per_c
    Targ_Stg_df['40%'] = V40per_c
    Targ_Stg_df['45%'] = V45per_c
    Targ_Stg_df['50%'] = V50per_c
    Targ_Stg_df['60%'] = V60per_c

    
    #Now calculate the Lake stage corresponding to a defined % chance 6/1 stage falls <11
    #Calculated using reverse routing of SFWMM-simulated Dsto (assumed no baseflow or ews) starting at 11' on 6/1.
    #Note: Values for June2-Sep30 are set to -999 since wet season probabilities are too uncertain.
    CalEst_baseflow = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    #Code the VLOOKUP(date, MonthlyBaseflowTable,3) for S79BS
    VLOOKUP2 = Basin_RO_Daily['BS-C43RO']
    VLOOKUP2_c = [x for x in VLOOKUP2 if ~np.isnan(x)]
####################################################################################################################
    #Read External TP loadings into Lake Okeechobee (mg) 
    Load_ext = pd.read_csv('./Data/LO_External_Loadings_3MLag_%s.csv'%Pre_defined_Variables.Schedule)
    #Lake Okeechobee inflows (m3/day)
    Q_in = pd.read_csv('./Data/LO_Inflows_BK_%s.csv'%Pre_defined_Variables.Schedule)

####################################################################################################################
    Stage_LO = np.zeros(n_rows,dtype = object)
    Stage_LO[0] = Pre_defined_Variables.begstageCS
    #Assume that Lake Okeechobee stage at second time step = initial stage!
    Stage_LO[1] = Pre_defined_Variables.begstageCS
    WSM_Zone = np.zeros(n_rows,dtype = object)
    Max_Supply = np.zeros(n_rows,dtype = object)
    Daily_Supply_1 = np.zeros(n_rows,dtype = object)
    LOSA_Supply = np.zeros(n_rows,dtype = object)
    Cut_back = np.zeros(n_rows,dtype = object)
    Dem_N_Sup = np.zeros(n_rows,dtype = object)
    NI_Supply = np.zeros(n_rows,dtype = object)
    Zone_Code = np.zeros(n_rows,dtype = object)
    LO_Zone = np.zeros(n_rows,dtype = object)
    Zone_D_Trib = np.zeros(Counter)
    Zone_D_stage = np.zeros(Counter)
    Zone_D_Seas = np.zeros(Counter)
    Zone_D_MSeas = np.zeros(Counter)
    Zone_D_Branch_Code = np.zeros(Counter)
    Zone_D_Rel_Code = np.zeros(Counter)
    Zone_C_Trib = np.zeros(Counter)
    Zone_C_Seas = np.zeros(Counter)
    Zone_C_MSeas = np.zeros(Counter)
    Zone_C_MetFcast = np.zeros(Counter)
    Zone_C_Branch_Code = np.zeros(Counter)
    Zone_C_Rel_Code = np.zeros(Counter)
    Zone_B_Trib = np.zeros(Counter)
    Zone_B_Stage = np.zeros(Counter)
    Zone_B_Seas = np.zeros(Counter)
    Zone_B_Branch_Code = np.zeros(Counter)
    Zone_B_Rel_Code = np.zeros(Counter)
    DecTree_Relslevel = np.zeros(n_rows,dtype = object)
    DecTree_Relslevel[0] = np.nan
    DecTree_Relslevel[1] = np.nan
    DayFlags = np.zeros(n_rows,dtype = object)
    if startdate.month == LO_Model['date'].iloc[2].month and startdate.day == LO_Model['date'].iloc[2].day:
        X1 = 'SimDay1'
    elif begdateCS.year == LO_Model['date'].iloc[2].year and begdateCS.month == LO_Model['date'].iloc[2].month and begdateCS.day == LO_Model['date'].iloc[2].day:
        X1 = 'CS start date'
    else:
        X1 = LO_Model['date'].iloc[2]
    DayFlags[2] = X1
    PlsDay = np.zeros(n_rows,dtype = object)
    Release_Level = np.zeros(n_rows,dtype = object)
    dh_7days = np.zeros(n_rows,dtype = object)
    ZoneCodeminus1Code = np.zeros(n_rows,dtype = object) 
    ZoneCodeCode = np.zeros(n_rows,dtype = object) 
    Fraction_of_Zone_height = np.zeros(n_rows,dtype = object)
    ReLevelCode_1 = np.zeros(n_rows,dtype = object)
    ReLevelCode_2 = np.zeros(n_rows,dtype = object)
    ReLevelCode_3_S80 = np.zeros(n_rows,dtype = object)
    S80_Mult = np.zeros(n_rows,dtype = object)
    S80_Mult_2 = np.zeros(n_rows,dtype = object)
    S80RS = np.zeros(n_rows,dtype = object)
    S308RG1 = np.zeros(n_rows,dtype = object)
    Sum_S308RG1 = np.zeros(n_rows,dtype = object)
    S80BS = np.zeros(n_rows,dtype = object)
    S308BK = np.zeros(n_rows,dtype = object)
    ROeast = np.zeros(n_rows,dtype = object)
    S308BS = np.zeros(n_rows,dtype = object)
    Sum_S308BK = np.zeros(n_rows,dtype = object)
    S308RG_Code = np.zeros(n_rows,dtype = object)
    S308RG = np.zeros(n_rows,dtype = object)
    S80 = np.zeros(n_rows,dtype = object)            
    ReLevelCode_3_S77 = np.zeros(n_rows,dtype = object)
    S77_Mult = np.zeros(n_rows,dtype = object)
    THC_Class_normal_or_above = np.zeros(Count_AP)
    Lake_O_Stage_AP = np.zeros(Count_AP)
    Lake_O_Schedule_Zone = np.zeros(Count_AP)
    LStgCorres = np.zeros(Count_AP)
    LowChance_Check = np.zeros(Count_AP)
    S77RS_AP = np.zeros(Count_AP)
    S77BS_AP = np.zeros(Count_AP)
    S77RS_Pre_AP_S77_Baseflow = np.zeros(Count_AP)
    Forecast_D_Sal = np.zeros(Count_AP)
    n30d_mavg = np.zeros(Count_AP)
    n30davgForecast = np.zeros(Count_AP)
    LORS08_bf_rel = np.zeros(Count_AP)
    LDS_LC6_1 = np.zeros(Count_AP)
    S_O = np.zeros(Count_AP)
    All_4 = np.zeros(Count_AP)
    Sabf = np.zeros(Count_AP)
    Swbf = np.zeros(Count_AP)
    Swbu = np.zeros(Count_AP)
    All_4andStage = np.zeros(Count_AP)
    All_4andStagein = np.zeros(Count_AP)
    P_AP_BF_Stg = np.zeros(Count_AP)
    Logic_test_1 = np.zeros(Count_AP)
    Post_Ap_Baseflow = np.zeros(Count_AP)
    S77RSplusPreAPS77bsf = np.zeros(Count_AP)
    AndEstNeedsLakeWater = np.zeros(Count_AP)
    Choose_PAPEWS_1 = np.zeros(Count_AP)
    Choose_PAPEWS_2 = np.zeros(Count_AP)
    AndLowChance61stagelessth11 = np.zeros(Count_AP)
    ATHCnora = np.zeros(Count_AP)
    Post_AP_EWS = np.zeros(Count_AP)
    Post_AP_Baseflow_EWS_cfs = np.zeros(Count_AP)
    S77BSAP = np.zeros(n_rows,dtype = object)
    S77_Mult_2 = np.zeros(n_rows,dtype = object)
    S77RS = np.zeros(n_rows,dtype = object)
    Sum_S77RS = np.zeros(n_rows,dtype = object)
    S77BK = np.zeros(n_rows,dtype = object)
    ROwest = np.zeros(n_rows,dtype = object)
    S79BS = np.zeros(n_rows,dtype = object)
    S77BS = np.zeros(n_rows,dtype = object)
    S77EWS = np.zeros(n_rows,dtype = object)
    S77REG = np.zeros(n_rows,dtype = object)
    S79 = np.zeros(n_rows,dtype = object) 
    TotRegEW = np.zeros(n_rows,dtype = object) 
    Choose_WCA = np.zeros(n_rows,dtype = object)
    RegWCA = np.zeros(n_rows,dtype = object)
    Choose_L8C51 = np.zeros(n_rows,dtype = object)
    RegL8C51 = np.zeros(n_rows,dtype = object)
    TotRegSo = np.zeros(n_rows,dtype = object) 
    Stage2ar = np.zeros(n_rows,dtype = object)
    Stage2marsh = np.zeros(n_rows,dtype = object) 
    RF = np.zeros(n_rows,dtype = object)
    ET = np.zeros(n_rows,dtype = object) 
    Choose_WSA_1 = np.zeros(n_rows,dtype = object) 
    Choose_WSA_2 = np.zeros(n_rows,dtype = object) 
    WSA_MIA = np.zeros(n_rows,dtype = object) 
    WSA_NNR = np.zeros(n_rows,dtype = object) 
    DSto = np.zeros(n_rows,dtype = object) 
    StartStorage = Stg_Sto_Ar.stg2sto(Pre_defined_Variables.startstage,0)
    Storage = np.zeros(n_rows,dtype = object)
    Storage[0] = StartStorage
    Storage[1] = StartStorage
    
    ##Here, I will insert the Storage Deviaiton Values as Input!
    Storage_dev = Data.Stroage_dev_df['DS_dev'] 
    #Create a Choose Function for AP Post Baseflow
    if Pre_defined_Variables.Opt_AdapProt == 1:
        C = 450
    elif Pre_defined_Variables.Opt_AdapProt == 2:
        C = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    Choose_1 = C
##############################################################################################################
    L_ext = Load_ext['TP_Loads_In_mg'] #mg
    Atm_Dep_N = TP_Variables.N_Per * Load_ext['Atm_Loading_mg']
    Atm_Dep_S = TP_Variables.S_Per * Load_ext['Atm_Loading_mg']
    # Q_Out = pd.read_csv('./Data/Outflows_consd_20082018.csv')
    # C_rain = 10.417 #TP Rainfall Concentration (µg P L-1 = mg P /m3) 
    # L_drdep = 0.0385 # mg P / m2 / day 
    # Atm_Dep_N = TP_Variables.N_Per * (C_rain*RF_Vol*1233.48 + L_drdep*LO_Area*4046.85642)
    # Atm_Dep_S = TP_Variables.S_Per * (C_rain*RF_Vol*1233.48 + L_drdep*LO_Area*4046.85642)
    # Atm_Dep_N = TP_Variables.N_Per*(18/365)*LO_Area*4046.85642 #Based on data presented by Curtis Pollman, the Lake Okeechobee Technical Advisory Committee (2000) recommended that 18 mgP/m2-yr is an appropriate atmospheric loading of phosphorus over the open lake. 
    # Atm_Dep_S = TP_Variables.S_Per*(18/365)*LO_Area*4046.85642
    #Read Shear Stress driven by Wind Speed
    Wind_ShearStr = pd.read_csv('./Data/WindShearStress_%s.csv'%Pre_defined_Variables.Schedule)
    W_SS = Wind_ShearStr['ShearStress'] #Dyne/cm2
    nu_ts = pd.read_csv('./Data/nu_%s.csv'%Pre_defined_Variables.Schedule)
    LO_BL = 0.5 # m (Bed Elevation of LO)
    # LO_WD = pd.to_numeric(Stage_Storage['Stage_m'])-LO_BL
    g = 9.8 #m/s2 gravitational acceleration
    Cal_Res = pd.read_csv('C:/Osama_PC/LOONE/Model/LOONE_Model/Data/nondominated_Sol_var.csv')
    Par = Cal_Res['Par']
    d_c = Par[20] # m (particle diameter 10 microm /1E6 to convert to m) clay
    d_s = Par[21] # m sand
    nu_d = nu_ts['nu']
    # LO_Temp = 1.0034/1E6 # m2/s (kinematic viscosity of water at T = 20 C) 
    # water_density = 1 # g/cm3
    # a = 20.0
    # n = 0.9
    # b = 2.5
    # m = 1.2
    R = 1.65 #submerged specific gravity (1.65 for quartz in water)
    C_1_c = Par[16]
    C_2_c = Par[17]
    C_1_s = Par[18]
    C_2_s = Par[19]
    #Parameters associated with sediment resuspension
    E_0 = 1E-4
    E_1 = 2
    E_2 = 3
    Crtcl_ShStr = Par[22] #0.32 #Dyne/cm2
    Td = Par[23] #days
    n_rows = len(Load_ext.index)
    L_ext_M = np.zeros(n_rows,dtype = object)
    Q_N2S = np.zeros(n_rows,dtype = object)
    # Stage_LO = Stage_Storage['Stage_ft']
    # Storage = Stage_Storage['Storage_acft']
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
    
    # TP_N_to_S = np.zeros(n_rows,dtype = object)
    # TP_Out = np.zeros(n_rows,dtype = object)
    # L_Ext_mgperm3 = np.zeros(n_rows,dtype = object)
    Q_I = Q_in['Flow_cmd']
    Q_I_M = np.zeros(n_rows,dtype = object)
    Q_O = np.zeros(n_rows,dtype = object)
    # Indust_O = pd.read_csv('./Data/INDUST_Outflow_20082018.csv')
    # Q_O = Q_Out['Total_Outflows_acft'] * 1233.48 + Indust_O['INDUST_cmd']
    Q_O_M = np.zeros(n_rows,dtype = object)
    
    P_Load_Cal = np.zeros(n_rows,dtype = object)
    P_Load_StL = np.zeros(n_rows,dtype = object)
    P_Load_South = np.zeros(n_rows,dtype = object)
    
    #Ferguson, R. I., and Church, M. (2004).
    # v_settle_N_c = (R*g*d_c**2)/(C_1_c*nu+(0.75*C_2_c*R*g*d_c**3)**0.5)
    # v_settle_N_s = (R*g*d_s**2)/(C_1_s*nu+(0.75*C_2_s*R*g*d_s**3)**0.5)
    # v_settle_N = v_settle_N_c*((TP_Variables.A_Mud_N+TP_Variables.A_Peat_N)/TP_Variables.A_N) + v_settle_N_s*((TP_Variables.A_Sand_N + TP_Variables.A_Rock_N)/TP_Variables.A_N)

    # v_settle_S_c = (R*g*d_c**2)/(C_1_c*nu+(0.75*C_2_c*R*g*d_c**3)**0.5)
    # v_settle_S_s = (R*g*d_s**2)/(C_1_s*nu+(0.75*C_2_s*R*g*d_s**3)**0.5)
    # v_settle_S = v_settle_S_c*((TP_Variables.A_Mud_S+TP_Variables.A_Peat_S)/TP_Variables.A_S) + v_settle_S_s*((TP_Variables.A_Sand_S + TP_Variables.A_Rock_S)/TP_Variables.A_S)

    v_settle_N_c = np.zeros(n_rows,dtype = object)
    v_settle_N_s = np.zeros(n_rows,dtype = object)
    v_settle_N = np.zeros(n_rows,dtype = object)
    v_settle_S_c = np.zeros(n_rows,dtype = object)
    v_settle_S_s = np.zeros(n_rows,dtype = object)
    v_settle_S = np.zeros(n_rows,dtype = object)

    # v_settle_N = np.zeros(n_rows,dtype = object)
    # v_settle_S = np.zeros(n_rows,dtype = object)   
    #####################################################################################################
    ##Initial Values##
    #S.A. is calculated based on the Lake's previous time step Stage, but for the S.A. at i=0 I used same time step Stage!
    Stage2ar[0] = Stg_Sto_Ar.stg2ar(Stage_LO[0],0)
    Stage2ar[1] = Stg_Sto_Ar.stg2ar(Stage_LO[1],0)
    Storage[0] = StartStorage #ac-ft
    Storage[1] = Stg_Sto_Ar.stg2sto(Stage_LO[1],0) #ac-ft
    Q_O[0] = 0
    Q_O[1] = 46485 #cmd
    #TP_MassBalanceModel Initial Values.
    TP_Lake_N[0] = 225 #mg/m3
    TP_Lake_S[0] = 275 #mg/m3
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
    
######################################################################################################################################    
    Zone_Code[0] = LO_FNs.Zone_Code(Stage_LO[0],df_WSMs['A'].iloc[0],df_WSMs['B'].iloc[0],df_WSMs['C'].iloc[0],df_WSMs['D3'].iloc[0],df_WSMs['D2'].iloc[0],df_WSMs['D1'].iloc[0],df_WSMs['D0'].iloc[0],df_WSMs['WSM1'].iloc[0])
    LO_Zone[0] = LO_FNs.LO_Zone(Zone_Code[0])
    for i in range(n_rows-2):
        WSM_Zone[i+2] = LO_FNs.WSM_Zone(Stage_LO[i+1],df_WSMs.at[i+1, 'WSM4'],df_WSMs.at[i+1, 'WSM3'],df_WSMs.at[i+1, 'WSM2'],df_WSMs.at[i+1, 'WSM1'])
    #Calculate Daily Maximum Water Supply
    # Note that in LOSA_dmd we used (i) because this file starts from 1/1/2008 so i at this point =0.
    #Cutbacks are determined based on the WSM Zone. 
        Max_Supply[i+2] = LO_FNs.Max_Supply(WSM_Zone[i+2],LOSA_dmd.at[i, 'Daily_demand'],Pre_defined_Variables.Z1_cutback,Pre_defined_Variables.Z2_cutback,Pre_defined_Variables.Z3_cutback,Pre_defined_Variables.Z4_cutback)
    #Actual Daily Water supply
        LOSA_Supply[i+2] = LO_FNs.LOSA_Supply(WSM_Zone[i+2],LO_Model.at[i+2, 'LOSA_dmd_SFWMM'],Max_Supply[i+2],Pre_defined_Variables.Opt_LOSAws)
    # NetInflow - LOSA Supply
        NI_Supply[i+2] = LO_Model.at[i+2, 'Net_Inflow'] - LOSA_Supply[i+2]
    #TODO Note: for the pass statement, We will read the Daily Water supply from the SFWMM as an input.
    #Calculate the cutback where Cutback = Demand - Supply
        ctbk = LO_Model.at[i+2, 'LOSA_dmd_SFWMM'] - LOSA_Supply[i+2]        
        Cut_back[i+2] = ctbk
    #Calculate percentage of the demand that is not supplied for each day
        if LO_Model.at[i+2, 'LOSA_dmd_SFWMM'] == 0:
            DNS = 0
        else:
            DNS = (Cut_back[i+2] / LO_Model.at[i+2, 'LOSA_dmd_SFWMM'])*100
        Dem_N_Sup[i+2] = DNS
    # Calculate the Zone Code
    #Note that to calculate the Zone Code in Dec 31 2020 we needed the WSM and breakpoint zones in 1/1/2021!
    #Note Also that i = 0 in Stage indicates Dec 30 1964 while i = 0 in df_WSMs indicates Dec 31 1964!
        Zone_Code[i+1] = LO_FNs.Zone_Code(Stage_LO[i+1],df_WSMs.at[i+1, 'A'],df_WSMs.at[i+1, 'B'],df_WSMs.at[i+1, 'C'],df_WSMs.at[i+1, 'D3'],df_WSMs.at[i+1, 'D2'],df_WSMs.at[i+1, 'D1'],df_WSMs.at[i+1, 'D0'],df_WSMs.at[i+1, 'WSM1'])
    #Generate the Zone Column based on the corresponding Zone Code.
        LO_Zone[i+1] = LO_FNs.LO_Zone(Zone_Code[i+1])
        #
        Zone_D_Trib[i] = Dec_Tree_FNs.Zone_D_Trib(TC_LONINO_df.at[i, 'Tributary_Condition'],Pre_defined_Variables.Opt_NewTree)
        #
        Zone_D_stage[i] = Dec_Tree_FNs.Zone_D_stage(Stage_LO[i+1],df_WSMs.at[i, 'C-b'])
        #
        Zone_D_Seas[i] = Dec_Tree_FNs.Zone_D_Seas(TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],Zone_D_Trib[i],Pre_defined_Variables.Opt_NewTree)
        #
        Zone_D_MSeas[i] = Dec_Tree_FNs.Zone_D_MSeas(TC_LONINO_df.at[i, 'LONINO_MultiSeasonal_Classes'])
        #
        ZDBC = Zone_D_Trib[i]*1000 + Zone_D_stage[i]*100 + Zone_D_Seas[i]*10 + Zone_D_MSeas[i]*1
        Zone_D_Branch_Code[i] = ZDBC
        #
        Zone_D_Rel_Code[i] = Dec_Tree_FNs.Zone_D_Rel_Code(Zone_D_Branch_Code[i],Pre_defined_Variables.Opt_DecTree)  
        #
        Zone_C_Trib[i] = Dec_Tree_FNs.Zone_C_Trib(TC_LONINO_df.at[i, 'Tributary_Condition'],Pre_defined_Variables.Opt_NewTree)
        #
        Zone_C_Seas[i] = Dec_Tree_FNs.Zone_C_Seas(TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],Pre_defined_Variables.Opt_NewTree)  
        #
        Zone_C_MSeas[i] = Dec_Tree_FNs.Zone_C_MSeas(TC_LONINO_df.at[i, 'LONINO_MultiSeasonal_Classes'])
        #
        Zone_C_MetFcast[i] = Dec_Tree_FNs.Zone_C_MetFcast(Zone_C_Seas[i],TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],Pre_defined_Variables.Zone_C_MetFcast_Indicator)
        #
        ZCBC = Zone_C_Trib[i]*1000 + Zone_C_MetFcast[i]*100 + Zone_C_Seas[i]*10 + Zone_C_MSeas[i]*1
        Zone_C_Branch_Code[i] = ZCBC
        #
        Zone_C_Rel_Code[i] = Dec_Tree_FNs.Zone_C_Rel_Code(Zone_C_Branch_Code[i],Pre_defined_Variables.Opt_DecTree)   
        #
        Zone_B_Trib[i] = Dec_Tree_FNs.Zone_B_Trib(TC_LONINO_df.at[i, 'Tributary_Condition'],Pre_defined_Variables.Opt_NewTree)
        #
        Zone_B_Stage[i] = Dec_Tree_FNs.Zone_B_Stage(Stage_LO[i+1],Seasons.at[i, 'Season'])
        #
        Zone_B_Seas[i] = Dec_Tree_FNs.Zone_B_Seas(TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'])
        #
        ZBBC = Zone_B_Trib[i]*1000 + Zone_B_Stage[i]*100 + DecTree_df.at[i, 'Zone_B_MetFcast']*10 + Zone_B_Seas[i]*1
        Zone_B_Branch_Code[i] = ZBBC
        #
        Zone_B_Rel_Code[i] = Dec_Tree_FNs.Zone_B_Rel_Code(Zone_B_Branch_Code[i],Pre_defined_Variables.Opt_DecTree)
        #Generate the DecTree Release level
        DecTree_Relslevel[i+2] = LO_FNs.DecTree_Relslevel(Zone_Code[i+1],Zone_D_Rel_Code[i],Zone_C_Rel_Code[i],Zone_B_Rel_Code[i])
        #Calculate DayFlags
        if i >= 3:
            if startdate.month == LO_Model.at[i, 'date'].month and startdate.day == LO_Model.at[i, 'date'].day and (Pre_defined_Variables.CSflag == 0 or startdate.year == LO_Model.at[i, 'date'].year):
                X2 = 'SimDay1'
            else:
                X2 = LO_Model.at[i, 'date'].date()
            DayFlags[i] = X2
        #
        PlsDay[i+2] = LO_FNs.PlsDay(DayFlags[i+2],DecTree_Relslevel[i+2],Pre_defined_Variables.PlsDay_Switch)
        #Calculate the Release Level Column
        Release_Level[i+2] = LO_FNs.Release_Level(Stage_LO[i+1],TC_LONINO_df.at[i, 'Tributary_Condition'],PlsDay[i+2],Zone_Code[i+1],DecTree_Relslevel[i+2],Pre_defined_Variables.MaxQstgTrigger)
        #Calculating the 7day dh Column!
        if i >= 6:
            dh = Stage_LO[i+1] - Stage_LO[i-6]
            dh_7days[i+1] = dh
    #Generate the fraction of zone height column
    #First I will generate A Code column based on the (ZoneCode-1)
        ZoneCodeminus1Code[i+1] = LO_FNs.ZoneCodeminus1Code(Zone_Code[i+1],df_WSMs.at[i+1, 'WSM1'],df_WSMs.at[i+1, 'D0'],df_WSMs.at[i+1, 'D1'],df_WSMs.at[i+1, 'D2'],df_WSMs.at[i+1, 'D3'],df_WSMs.at[i+1, 'C'],df_WSMs.at[i+1, 'B'],df_WSMs.at[i+1, 'A'])
    #Then, I will generate A Code column based on the (ZoneCode)
        ZoneCodeCode[i+1] = LO_FNs.ZoneCodeCode(Zone_Code[i+1],df_WSMs.at[i+1, 'WSM1'],df_WSMs.at[i+1, 'D0'],df_WSMs.at[i+1, 'D1'],df_WSMs.at[i+1, 'D2'],df_WSMs.at[i+1, 'D3'],df_WSMs.at[i+1, 'C'],df_WSMs.at[i+1, 'B'],df_WSMs.at[i+1, 'A'])
    #Now, Calculate the Fraction Zone Height Column!
        Fraction_of_Zone_height[i+1] = LO_FNs.Fraction_of_Zone_height(Zone_Code[i+1],Stage_LO[i+1],ZoneCodeminus1Code[i+1],ZoneCodeCode[i+1])
    #Define three Codes based on the RelsLevel that will be used to calculate S80_Mult!
        ReLevelCode_1[i+2] = LO_FNs.ReLevelCode_1(Release_Level[i+2],Pre_defined_Variables.dstar_D1,Pre_defined_Variables.dstar_D2,Pre_defined_Variables.dstar_D3,Pre_defined_Variables.dstar_C,Pre_defined_Variables.dstar_B)
        #
        ReLevelCode_2[i+2] = LO_FNs.ReLevelCode_2(Release_Level[i+2],Pre_defined_Variables.astar_D1,Pre_defined_Variables.astar_D2,Pre_defined_Variables.astar_D3,Pre_defined_Variables.astar_C,Pre_defined_Variables.astar_B)
        #
        ReLevelCode_3_S80[i+2] = LO_FNs.ReLevelCode_3_S80(Release_Level[i+2],Pre_defined_Variables.bstar_S80_D1,Pre_defined_Variables.bstar_S80_D2,Pre_defined_Variables.bstar_S80_D3,Pre_defined_Variables.bstar_S80_C,Pre_defined_Variables.bstar_S80_B)
    #Multiplier derived from specified user option & parameters that depends on season, zone, and factors (dstar,astar & bstar).
        S80_Mult[i+2] = LO_FNs.S80_Mult(Seasons.at[i, 'Season'],Seasons.at[i, 'Month'],dh_7days[i+1],ReLevelCode_1[i+2],Fraction_of_Zone_height[i+1],ReLevelCode_2[i+2],ReLevelCode_3_S80[i+2],Pre_defined_Variables.Opt_QregMult)
    #Calculate the S80Mult2
    #Adjusts S80mult to allow 10-day pulses that start in May (or October) and end in June (or November) to use the same multiplier for the entire 10-day pulse.
    #Same assumption is made by the SFWMM.
        S80_Mult_2[i+2] = LO_FNs.S80_Mult_2(LO_Model.at[i+2, 'date'].month,LO_Model.at[i+2, 'date'].day,PlsDay[i+2],S80_Mult[i+2-PlsDay[i+2]],S80_Mult[i+2],Pre_defined_Variables.Opt_QregMult)
    #Calculate Reg Schedule Release without BaseFlow (S80RS)
        
        S80RS[i+2] = LO_FNs.S80RS(Release_Level[i+2],Data.S80_RegRelRates.at[0, 'Zone_D1'],S80avgL1,Data.Pulses.at[PlsDay[i+2]-1 if PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-80_L1'],S80_Mult_2[i+2],Data.CE_SLE_turns.at[LO_Model.at[i+2, 'date'].year-Pre_defined_Variables.startyear, 'SLEturn'],Data.S80_RegRelRates.at[0, 'Zone_D2'],S80avgL2,Data.Pulses.at[PlsDay[i+2]-1 if PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-80_L2'],Data.S80_RegRelRates.at[0, 'Zone_D3'],S80avgL3,Data.Pulses.at[PlsDay[i+2]-1 if PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-80_L3'],Data.S80_RegRelRates.at[0, 'Zone_C'],Data.S80_RegRelRates.at[0, 'Zone_B'],Data.S80_RegRelRates.at[0, 'Zone_A'])
    #Calculate S308RG1
        S308RG1[i+2] = max(0,S80RS[i+2]-LO_Model.at[i+2, 'C44RO'])
    #Calculale accumulated S308RG1
        Sum_S308RG1[i+2] = LO_FNs.Sum_S308RG1(LO_Model.at[i+2, 'date'].day,S308RG1[i+2])
    #S80 Baseflow Target
    #accounts for basin runoff and accumulated S308reg.  
    #If either is enough to meet baseflow requirement for the month, then S80BS target = 0. 
    #NOTE: Starting with v6.10, the MultiSeasonal LONINO constraint on Baseflow was removed.  
    #The constraint tested if the MSLONINO was dry, then no Baseflow.  
    #The constraint existed because the SFWMM had to use the enviro ws code to approximate the baseflow provision of LORS08.
    #In 2013 the SFWMM was modified to use the regulatory discharge code for baseflow.  
    #So the constraint was removed from the LOOPS Model (v6.10 and after) to be consistent with the new version of the SFWMM (circa Jan 2013).
    #Note that relaxing this constraint did not affect the baseflow for S80 since it was set to zero for LORS08 consistency testing with the SFWMM.
    #See similar note for S79.
    #Now, calculate the S80BS
        S80BS[i+2] = LO_FNs.S80BS(Release_Level[i+2],Sum_S308RG1[i+2],VLOOKUP1_c[i],StlEst_baseflow,Pre_defined_Variables.Option_S80Baseflow)
    #Calculate S308BK
    #C44 Basin runoff can backflow if the FIRST TWO conditions are true (3rd condition was needed for older versions of SFWMM which handled baseflow using the environmental ws code):
    #(1) lake stage falls below backflow line,
    #(2) no Lake reg or baseflow is released,
    #(3) accum backflow for the month is < max backflow
    #   where max backflow = monthly runoff vol - S80 baseflow target
        S308BK[i+2] = LO_FNs.S308BK(Stage_LO[i+1],df_WSMs.at[i+1, 'D1'],S308RG[i+1],LO_Model.at[i+2, 'C44RO'],Data.SFWMM_Daily_Outputs.at[i+2, 'S308BK'],Pre_defined_Variables.Opt_S308,Pre_defined_Variables.S308BK_Const,Pre_defined_Variables.S308_BK_Thr)
    #Calculate ROeast column.
        ROeast[i+2] = LO_Model.at[i+2, 'C44RO'] - S308BK[i+2]
    # Calculate S308BS (Baseflow Component needed from LOK)
        S308BS[i+2] = LO_FNs.S308BS(S80BS[i+2],S308RG1[i+2],ROeast[i+2],Pre_defined_Variables.Option_S80Baseflow)
    #Calculate Monthly sum of S308BK
        Sum_S308BK[i+2] = LO_FNs.Sum_S308BK(LO_Model.at[i+2, 'date'].day,S308BK[i+2])
    #Calculate S308RG
    #Define a Choose function to be used to calculate S308RG.
        S308RG_Code[i+2] = LO_FNs.S308RG_Code(S308RG1[i+2],S308BS[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S308RG'],Data.SFWMM_Daily_Outputs.at[i+2, 'STEST'],Pre_defined_Variables.Option_RegS77S308)
    #
        # S308RG[i+2] = LO_FNs.S308RG(S308RG_Code[i+2],Data.SFWMM_Daily_Outputs.at[i+2,'S308RG'],Data.SFWMM_Daily_Outputs.at[i+2,'STEST'],Pre_defined_Variables.Opt_S308,Pre_defined_Variables.S308RG_Const)
        if Stage_LO[i+1] >= 18:
            S308RG[i+2] = 7200
        elif Stage_LO[i+1] <= 8:
            S308RG[i+2] = 0
        elif (TP_Lake_S[i] <= P_1) and (date_rng_6[i+2].month in [1,2,3,4,11,12]):
            S308RG[i+2] = S308_DV.iloc[(date_rng_6[i+2].month)-1]
        elif (TP_Lake_S[i] <= P_2) and (date_rng_6[i+2].month in [5,6,7,8,9,10]):
            S308RG[i+2] = S308_DV.iloc[(date_rng_6[i+2].month)-1]
        else:
            S308RG[i+2] = 0

    #Calculate S80 Column
    #Note that S308EWS should be included in this Equation if used!
        S80[i+2] = LO_FNs.S80(ROeast[i+2],S308RG[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S80'],Pre_defined_Variables.S80_Const)
        #FIXME Here I will just add 1 at 1/1/2008 to mock LOOPS because of S308EWS
        #input('Please, Check S80 value at 1/1/2008 due to S308EWS ignorance')
        #Calculate S77mult Column.
        #Multiplier derived from specified user option & parameters (see notes in Assumptions sheet). 
        #Multiplier depends on season, zone, and factors (dstar,astar & bstar).
        #Note:  this logic was added to v5.96.  
        #Previous versions, particularly v5.5 and earlier, applied the multiplier with formula conditions within the S80RS and S77RS formulas.
        #The new logic determines the multiplier in separate (preceding) column formulas.
        #There are some slight differences between the methods:
        #1. Old method applied multiplier to zones D1,D2,D3 and C.  
        #New method is user-input for release levels corresponding to D1,D2,D3,C & B.
        #2. If stage fell below D1 during 10-day pulse, the old method did not apply multiplier to flows during times stage was in D0.  
        #New method applies multiplier to release type and reduces L1 pulses even for latter days of 10-day pulse when stage falls within D0.
        #We will use the same three Codes based on the RelsLevel that we used with S80Mult to calculate S77_Mult!
        #Except we will adjust the third code to fit S77!
        ReLevelCode_3_S77[i+2] = LO_FNs.ReLevelCode_3_S77(Release_Level[i+2],Pre_defined_Variables.bstar_S77_D1,Pre_defined_Variables.bstar_S77_D2,Pre_defined_Variables.bstar_S77_D3,Pre_defined_Variables.bstar_S77_C,Pre_defined_Variables.bstar_S77_B)
        #
        S77_Mult[i+2] = LO_FNs.S77_Mult(Seasons.at[i, 'Season'],Seasons.at[i, 'Month'],dh_7days[i+1],ReLevelCode_1[i+2],Fraction_of_Zone_height[i+1],ReLevelCode_2[i+2],ReLevelCode_3_S77[i+2],Pre_defined_Variables.Opt_QregMult)
        #Calculate S77mult2
        #Adjusts S77mult to allow 10-day pulses that start in May (or October) and end in June (or November) 
        #to use the same multiplier for the entire 10-day pulse.  Same assumption is made by the SFWMM.
        S77_Mult_2[i+2] = LO_FNs.S77_Mult_2(LO_Model.at[i+2, 'date'].month,LO_Model.at[i+2, 'date'].day,PlsDay[i+2],S77_Mult[i+2-PlsDay[i+2]],S77_Mult[i+2],Pre_defined_Variables.Opt_QregMult)
    #Calculate S77RS (Reg Schedule Release without BaseFlow).
        S77RS[i+2] = LO_FNs.S77RS(Release_Level[i+2],Data.S77_RegRelRates.at[0, 'Zone_D1'],S77avgL1,Data.Pulses.at[PlsDay[i+2]-1 if PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-77_L1'],S77_Mult_2[i+2],LO_Model.at[i+2, 'C43RO'],Data.CE_SLE_turns.at[LO_Model.at[i+2, 'date'].year-Pre_defined_Variables.startyear, 'CEturn'],Data.S77_RegRelRates.at[0, 'Zone_D2'],S77avgL2,Data.Pulses.at[PlsDay[i+2]-1 if PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-77_L2'],Zone_Code[i+1],Data.S77_RegRelRates.at[0, 'Zone_D3'],S77avgL3,Data.Pulses.at[PlsDay[i+2]-1 if PlsDay[i+2]-1>=0 else len(Data.Pulses)-1, 'S-77_L3'],Data.S77_RegRelRates.at[0, 'Zone_C'],Data.S77_RegRelRates.at[0, 'Zone_B'],Data.S77_RegRelRates.at[0, 'Zone_A'],Pre_defined_Variables.Opt_S79RG)
    #Calculale accumulated S77RS
        Sum_S77RS[i+2] = LO_FNs.Sum_S77RS(LO_Model.at[i+2, 'date'].day,S77RS[i+2])
        #S-77 backflow logic:
        #Backflow can occur if:
        #1. LOK stage < threshold (eg, 11.1 ft)
        #2. Basin runoff > 0
        #3. No Lake O regulatory discharge at S-77
        #4. No Baseflow or EWS for Caloos Estuary at S-77 (previous day)
        #Note, the 4th condition was added beginning with v6.27.  The previous day is used for baseflow and EWS to avoid a circular reference since the Adaptive Protocol calculations occur after backflow is computed.  Future logic can be smarter, but will require revising order of operations.
        #Pre-v6.27 formula was:
        #=IF($AT$16=1,IF(AND(BO20<$AT$17,AR21=0),MAX(0,0.3*AO21-5),0),SFWMMdata4LOOPS!M11)
        S77BK[i+2] = LO_FNs.S77BK(Stage_LO[i+1],S77RS[i+2],S77BSAP[i+1],S77EWS[i+1],LO_Model.at[i+2, 'C43RO'],Data.SFWMM_Daily_Outputs.at[i+2, 'S77BK'],Pre_defined_Variables.S77BK_Switch,Pre_defined_Variables.S77BK_Threshold)
        #ROwest
        ROwest[i+2] = LO_Model.at[i+2, 'C43RO'] - S77BK[i+2]
        #FIXME
        #S79 Baseflow Target
        #accounts for basin runoff and accumulated S77reg.  If either is enough to meet baseflow requirement for the month, then S79BS target = 0. 
        #NOTE: Starting with v6.10, the MultiSeasonal LONINO constraint on Baseflow was removed.  
        #The constraint tested if the MSLONINO was dry, then no Baseflow.  The constraint existed because the SFWMM had to use the enviro ws code to approximate the baseflow provision of LORS08.  
        #In 2013 the SFWMM was modified to use the regulatory discharge code for baseflow.  So the constraint was removed from the LOOPS Model (v6.10 and after) to be consistent with the new version of the SFWMM (circa Jan 2013).
        #Note that relaxing this constraint reduced the number of days of zero baseflow by almost 600 days (about 5% less zero baseflow days).
        #Now, calculate the S79BS
        S79BS[i+2] = LO_FNs.S79BS(Release_Level[i+2],Sum_S77RS[i+2],VLOOKUP2_c[i],CalEst_baseflow,Pre_defined_Variables.Option_S77Baseflow)
    #S77BS
        S77BS[i+2] = LO_FNs.S77BS(S79BS[i+2],S77RS[i+2],ROwest[i+2],Pre_defined_Variables.Option_S77Baseflow)
    #############################################################################################
        #Define THC Class Normal or above
        if i < (n_rows-2):
            Post_Ap_Baseflow[i] = THC_Class(i,THC_Class_normal_or_above,Lake_O_Stage_AP,Lake_O_Schedule_Zone,LStgCorres,LowChance_Check,S77RS_AP,S77BS_AP,
              S77RS_Pre_AP_S77_Baseflow,Forecast_D_Sal,n30d_mavg,n30davgForecast,LORS08_bf_rel,LDS_LC6_1,S_O,All_4,
              Sabf,Swbf,Swbu,All_4andStage,All_4andStagein,P_AP_BF_Stg,Logic_test_1,Post_Ap_Baseflow,S77RSplusPreAPS77bsf,
              AndEstNeedsLakeWater,AndLowChance61stagelessth11,ATHCnora,Choose_PAPEWS_1,Choose_PAPEWS_2,Post_AP_EWS,
              Post_AP_Baseflow_EWS_cfs,AdapProt_df,Stage_LO,Zone_Code,df_WSMs,Targ_Stg_df,S77RS,S77BS,Data.Estuary_needs_water,
              Choose_1,WSM_Zone)['Post_Ap_Baseflow']
            Post_AP_EWS[i] = THC_Class(i,THC_Class_normal_or_above,Lake_O_Stage_AP,Lake_O_Schedule_Zone,LStgCorres,LowChance_Check,S77RS_AP,S77BS_AP,
              S77RS_Pre_AP_S77_Baseflow,Forecast_D_Sal,n30d_mavg,n30davgForecast,LORS08_bf_rel,LDS_LC6_1,S_O,All_4,
              Sabf,Swbf,Swbu,All_4andStage,All_4andStagein,P_AP_BF_Stg,Logic_test_1,Post_Ap_Baseflow,S77RSplusPreAPS77bsf,
              AndEstNeedsLakeWater,AndLowChance61stagelessth11,ATHCnora,Choose_PAPEWS_1,Choose_PAPEWS_2,Post_AP_EWS,
              Post_AP_Baseflow_EWS_cfs,AdapProt_df,Stage_LO,Zone_Code,df_WSMs,Targ_Stg_df,S77RS,S77BS,Data.Estuary_needs_water,
              Choose_1,WSM_Zone)['Post_AP_EWS']
        #######################################################################################################################################
        #Baseflow Component needed from LOK as modified using Adaptive Protocol logic.
        #Refer to sheet AdapProt.
        #NOTE: THE EWS COMPONENT WAS  SPLIT OUT FOR ACCOUNTING PURPOSES SINCE EWS IS NOT A REG RELEASE.  MAY2013.
        S77BSAP[i+2] = LO_FNs.S77BSAP(S77BS[i+2],Post_Ap_Baseflow[i],Pre_defined_Variables.Opt_AdapProt)
        #THE EWS COMPONENT WAS  SPLIT OUT FOR ACCOUNTING PURPOSES SINCE EWS IS NOT A REG RELEASE.  MAY2013.
        S77EWS[i+2] = LO_FNs.S77EWS(Post_AP_EWS[i],Data.SFWMM_Daily_Outputs.at[i+2, 'CAEST'],Pre_defined_Variables.S77EWS_Switch,Pre_defined_Variables.Opt_AdapProt)
        #S77REG
        #v6.10 (2013) modification excludes env ws release (previously lumped with S77BSAP).  
        #That new accounting variable, S77EWS, is handled separate from regulatory discharges.
        # S77REG[i+2] = LO_FNs.S77REG(S77RS[i+2],S77BSAP[i+2],Data.SFWMM_Daily_Outputs.at[i+2,'S77RG'],Pre_defined_Variables.S77REG_Switch,Pre_defined_Variables.Option_RegS77S308)
        if Stage_LO[i+1] >= 18:
            S77REG[i+2] = 7800
        elif Stage_LO[i+1] <= 8:
            S77REG[i+2] = 0
        elif (TP_Lake_S[i] <= P_1) and (date_rng_6[i+2].month in [1,2,3,4,11,12]):
            S77REG[i+2] = S77_DV.iloc[(date_rng_6[i+2].month)-1]
        elif (TP_Lake_S[i] <= P_2) and (date_rng_6[i+2].month in [5,6,7,8,9,10]):
            S77REG[i+2] = S77_DV.iloc[(date_rng_6[i+2].month)-1]
        else:
            S77REG[i+2] = 0

        S79[i+2] = LO_FNs.S79(S77REG[i+2],S77EWS[i+2],ROwest[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'S79'],Pre_defined_Variables.S79_Switch)
        #TotRegEW
        TotRegEW[i+2] = (S77REG[i+2] + S308RG[i+2])*1.9835
        #RegWCA
        #Define Choose function for RegWCA
        Choose_WCA[i+2] = LO_FNs.Choose_WCA(Data.SFWMM_Daily_Outputs.at[i+2, 'RegWCA'],Pre_defined_Variables.Option_RegWCA,Pre_defined_Variables.Constant_RegWCA)
    #Now define the RegWCA
        RegWCA[i+2] = min(Pre_defined_Variables.MaxCap_RegWCA , Pre_defined_Variables.Multiplier_RegWCA*Choose_WCA[i+2])
    #Define Choose function for RegL8C51
        Choose_L8C51[i+2] = LO_FNs.Choose_L8C51(Data.SFWMM_Daily_Outputs.at[i+2, 'RegL8C51'],Pre_defined_Variables.Option_RegL8C51,Pre_defined_Variables.Constant_RegL8C51)
    #Now define the RegWCA
        RegL8C51[i+2] = min(Pre_defined_Variables.MaxCap_RegL8C51 , Pre_defined_Variables.Multiplier_RegL8C51*Choose_L8C51[i+2])
        #TotRegSo
        #See the Setup Sheet.  User can choose from several options for testing southward releases.  
        #The common option is to use the same values that were simulated by the SFWMM that were used to derive the Net Inflow term.
        TotRegSo[i+2] = (RegWCA[i+2] + RegL8C51[i+2]) * 1.9835
    # Calculate stg2ar(ac)
        Stage2ar[i+2] = Stg_Sto_Ar.stg2ar(Stage_LO[i+1],0)
    # Calculate stg2mar(ac)
        Stage2marsh[i+2] = Stg_Sto_Ar.stg2mar(Stage_LO[i+1],0)
    #RF(af)
        RF[i+2] = Data.RF_Vol.at[i+2, 'RF_acft']
    #ET(af)
        ET[i+2] = LO_FNs.ET(Data.SFWMM_Daily_Outputs.at[i+2, 'et_dry'],Stage2ar[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'et_litoral'],Stage2marsh[i+2],Data.SFWMM_Daily_Outputs.at[i+2, 'et_open'],Data.ET_Vol.at[i+2, 'ETVol_acft'],Pre_defined_Variables.ET_Switch)
        #WS Augmentation
        #Water Supply Augmentation Logic:
        #IF stage is in Zone WSA2, THEN available flow to Lake is limited by structure capacity for WSA2.
        #IF stage is in Zone WSA1, THEN available flow to Lake is limited by structure capacity for WSA1.
        #IF stage is not in Zone WSA1 or WSA2, then no WSA flow to the Lake.
        #If flood control pumping (Qfc) is occurring on this same day, then WS Augmentation (Qws) is limited.
        #Qws=max(0,min(CAP,runoff)-Qfc)
        #And if WCA-3A stage is below prescribed elevation (eg, -2 feet below the traditional 9.5-10.5' regulation schedule), 
        #then no WS Augmentation (Qws) since this EAA runoff may benefit the Everglades.
        #AND, if the Seasonal LONINO is higher (wetter) than a user-specified value, 
        #then no WS Augmentation (Qws) to reduce chance of increasing near-future high Lake stages and regulatory discharges.
        #MIA (cfs)
        # EAA to Lake O Structure Inflow Capacities (cfs)
        #Choose_1 and 2 for WS Augmentation calculations    
        Choose_WSA_1[i+2] = LO_FNs.Choose_WSA_1(df_WSMs.at[i+2, 'WSM1'],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSAtrig2,Pre_defined_Variables.WSAoff2)   
        Choose_WSA_2[i+2] = LO_FNs.Choose_WSA_2(df_WSMs.at[i+2, 'WSM1'],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSAtrig1,Pre_defined_Variables.WSAoff1)    
        WSA_MIA[i+2] = LO_FNs.WSA_MIA(WCA_Stages_df.at[i, 'Are WCA stages too low?'],TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],Stage_LO[i+1],Choose_WSA_1[i+2],Data.EAA_MIA_RUNOFF.at[i, 'MIA'],Data.EAA_MIA_RUNOFF.at[i, 'S3PMP'],Choose_WSA_2[i+2],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSA_THC,Pre_defined_Variables.MIAcap2,Pre_defined_Variables.MIAcap1)
        #NNR WSA    
        WSA_NNR[i+2] = LO_FNs.WSA_NNR(WCA_Stages_df.at[i, 'Are WCA stages too low?'],TC_LONINO_df.at[i, 'LONINO_Seasonal_Classes'],Stage_LO[i+1],Choose_WSA_1[i+2],Data.EAA_MIA_RUNOFF.at[i, 'NNR'],Data.EAA_MIA_RUNOFF.at[i, 'S2PMP'],Choose_WSA_2[i+2],Pre_defined_Variables.Opt_WSA,Pre_defined_Variables.WSA_THC,Pre_defined_Variables.NNRcap2,Pre_defined_Variables.NNRcap1)
        #Simulated DSto (af)
        DSto[i+2] = NI_Supply[i+2] + RF[i+2] - ET[i+2] + 1.9835*(S308BK[i+2]\
                     + S77BK[i+2] + WSA_MIA[i+2] + WSA_NNR[i+2]\
                     - S77EWS[i+2]) - TotRegEW[i+2] - TotRegSo[i+2] + Storage_dev[i+2]
        #Storage (af)
        #Storage at first date
        Storage[i+2] = LO_FNs.Storage(DayFlags[i+2],Storage[i],StartStorage,Storage[i+1],DSto[i+2])
        #Lake Okeechobee Stage!
        Stage_LO[i+2] = LO_FNs.Stage_LO(Stg_Sto_Ar.stg2sto(Storage[i+2],1),Data.SFWMM_Daily_Outputs.at[i+2, 'EOD Stg(ft,NGVD)'],Pre_defined_Variables.Option_Stage)
#################################################################################################################################################################      
        Q_O[i+2] = (S77EWS[i+2] *0.028316847 + ((TotRegEW[i+2] + TotRegSo[i+2])/70.0456)) * 3600 * 24
  
        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48 #m3/d
            Q_O_M[i] = Q_O[i]
            L_ext_M[i] = L_ext[i] + Q_I_M[i] * TP_Lake_N[i]            
        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48 #m3/d
            Q_I_M[i] = Q_I[i]
            L_ext_M[i] = L_ext[i]
        Q_N2S[i] = (Q_I_M[i] + Q_O_M[i])/2
        Stage2ar[i+2] = Stg_Sto_Ar.stg2ar(Stage_LO[i+2],0)
        LO_WD[i] = Stage_LO[i]*0.3048 - LO_BL
        Lake_O_Storage_N[i] = Storage[i] * TP_Variables.N_Per * 4046.85642 * 0.305 #m3
        Lake_O_Storage_S[i] = Storage[i] * TP_Variables.S_Per * 4046.85642 * 0.305 #m3
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
       
        Sed_Resusp_M_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_M_N[i] if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_S_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_S_N[i] if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_R_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_R_N[i] if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_P_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_P_N[i] if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_M_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_M_S[i] if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_S_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_S_S[i] if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_R_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_R_S[i] if W_SS[i] > Crtcl_ShStr else 0
        Sed_Resusp_P_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_P_S[i] if W_SS[i] > Crtcl_ShStr else 0
       
        P_sed_M_N[i+1] = TP_MBFR.P_sed(Lake_O_A_M_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.K_decomp_M,v_settle_N[i]) - Sed_Resusp_M_N[i]*Lake_O_Storage_N[i]/Mass_sed_M_N if TP_MBFR.P_sed(Lake_O_A_M_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.K_decomp_M,v_settle_N[i]) - Sed_Resusp_M_N[i]*Lake_O_Storage_N[i]/Mass_sed_M_N > 0 else 0
        P_sed_S_N[i+1] = TP_MBFR.P_sed(Lake_O_A_S_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.K_decomp_S,v_settle_N[i]) - Sed_Resusp_S_N[i]*Lake_O_Storage_N[i]/Mass_sed_S_N if TP_MBFR.P_sed(Lake_O_A_S_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.K_decomp_S,v_settle_N[i]) - Sed_Resusp_S_N[i]*Lake_O_Storage_N[i]/Mass_sed_S_N > 0 else 0
        P_sed_R_N[i+1] = TP_MBFR.P_sed(Lake_O_A_R_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.K_decomp_R,v_settle_N[i]) - Sed_Resusp_R_N[i]*Lake_O_Storage_N[i]/Mass_sed_R_N if TP_MBFR.P_sed(Lake_O_A_R_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.K_decomp_R,v_settle_N[i]) - Sed_Resusp_R_N[i]*Lake_O_Storage_N[i]/Mass_sed_R_N > 0 else 0
        P_sed_P_N[i+1] = TP_MBFR.P_sed(Lake_O_A_P_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.K_decomp_P,v_settle_N[i]) - Sed_Resusp_P_N[i]*Lake_O_Storage_N[i]/Mass_sed_P_N if TP_MBFR.P_sed(Lake_O_A_P_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.K_decomp_P,v_settle_N[i]) - Sed_Resusp_P_N[i]*Lake_O_Storage_N[i]/Mass_sed_P_N > 0 else 0
        P_sed_M_S[i+1] = TP_MBFR.P_sed(Lake_O_A_M_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.K_decomp_M,v_settle_S[i]) - Sed_Resusp_M_S[i]*Lake_O_Storage_S[i]/Mass_sed_M_S if TP_MBFR.P_sed(Lake_O_A_M_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.K_decomp_M,v_settle_S[i]) - Sed_Resusp_M_S[i]*Lake_O_Storage_S[i]/Mass_sed_M_S > 0 else 0
        P_sed_S_S[i+1] = TP_MBFR.P_sed(Lake_O_A_S_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.K_decomp_S,v_settle_S[i]) - Sed_Resusp_S_S[i]*Lake_O_Storage_S[i]/Mass_sed_S_S if TP_MBFR.P_sed(Lake_O_A_S_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.K_decomp_S,v_settle_S[i]) - Sed_Resusp_S_S[i]*Lake_O_Storage_S[i]/Mass_sed_S_S > 0 else 0
        P_sed_R_S[i+1] = TP_MBFR.P_sed(Lake_O_A_R_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.K_decomp_R,v_settle_S[i]) - Sed_Resusp_R_S[i]*Lake_O_Storage_S[i]/Mass_sed_R_S if TP_MBFR.P_sed(Lake_O_A_R_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.K_decomp_R,v_settle_S[i]) - Sed_Resusp_R_S[i]*Lake_O_Storage_S[i]/Mass_sed_R_S > 0 else 0
        P_sed_P_S[i+1] = TP_MBFR.P_sed(Lake_O_A_P_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.K_decomp_P,v_settle_S[i]) - Sed_Resusp_P_S[i]*Lake_O_Storage_S[i]/Mass_sed_P_S if TP_MBFR.P_sed(Lake_O_A_P_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.K_decomp_P,v_settle_S[i]) - Sed_Resusp_P_S[i]*Lake_O_Storage_S[i]/Mass_sed_P_S > 0 else 0
        
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
        
        DIP_pore_M_N[i+1] = TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_N[i],DIP_Lake_N[i],J_des_M_N[i],J_ads_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.v_diff_M,TP_Variables.A_Mud_N,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) if TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_N[i],DIP_Lake_N[i],J_des_M_N[i],J_ads_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.v_diff_M,TP_Variables.A_Mud_N,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) > 0 else 0
        DIP_pore_S_N[i+1] = TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_N[i],DIP_Lake_N[i],J_des_S_N[i],J_ads_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.v_diff_S,TP_Variables.A_Sand_N,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) if TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_N[i],DIP_Lake_N[i],J_des_S_N[i],J_ads_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.v_diff_S,TP_Variables.A_Sand_N,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) > 0 else 0
        DIP_pore_R_N[i+1] = TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_N[i],DIP_Lake_N[i],J_des_R_N[i],J_ads_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.v_diff_R,TP_Variables.A_Rock_N,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) if TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_N[i],DIP_Lake_N[i],J_des_R_N[i],J_ads_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.v_diff_R,TP_Variables.A_Rock_N,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) > 0 else 0
        DIP_pore_P_N[i+1] = TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_N[i],DIP_Lake_N[i],J_des_P_N[i],J_ads_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.v_diff_P,TP_Variables.A_Peat_N,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) if TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_N[i],DIP_Lake_N[i],J_des_P_N[i],J_ads_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.v_diff_P,TP_Variables.A_Peat_N,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) > 0 else 0
        DIP_pore_M_S[i+1] = TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_S[i],DIP_Lake_S[i],J_des_M_S[i],J_ads_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.v_diff_M,TP_Variables.A_Mud_S,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) if TP_MBFR.DIP_pore(Θ_M,DIP_pore_M_S[i],DIP_Lake_S[i],J_des_M_S[i],J_ads_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.v_diff_M,TP_Variables.A_Mud_S,TP_Variables.K_decomp_M,TP_Variables.v_burial_M) > 0 else 0
        DIP_pore_S_S[i+1] = TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_S[i],DIP_Lake_S[i],J_des_S_S[i],J_ads_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.v_diff_S,TP_Variables.A_Sand_S,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) if TP_MBFR.DIP_pore(Θ_S,DIP_pore_S_S[i],DIP_Lake_S[i],J_des_S_S[i],J_ads_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.v_diff_S,TP_Variables.A_Sand_S,TP_Variables.K_decomp_S,TP_Variables.v_burial_S) > 0 else 0
        DIP_pore_R_S[i+1] = TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_S[i],DIP_Lake_S[i],J_des_R_S[i],J_ads_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.v_diff_R,TP_Variables.A_Rock_S,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) if TP_MBFR.DIP_pore(Θ_R,DIP_pore_R_S[i],DIP_Lake_S[i],J_des_R_S[i],J_ads_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.v_diff_R,TP_Variables.A_Rock_S,TP_Variables.K_decomp_R,TP_Variables.v_burial_R) > 0 else 0
        DIP_pore_P_S[i+1] = TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_S[i],DIP_Lake_S[i],J_des_P_S[i],J_ads_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.v_diff_P,TP_Variables.A_Peat_S,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) if TP_MBFR.DIP_pore(Θ_P,DIP_pore_P_S[i],DIP_Lake_S[i],J_des_P_S[i],J_ads_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.v_diff_P,TP_Variables.A_Peat_S,TP_Variables.K_decomp_P,TP_Variables.v_burial_P) > 0 else 0
        
        Settling_P_N[i] = TP_MBFR.Sett_P(TP_Lake_N[i], DIP_Lake_N[i], Lake_O_A_N[i], Lake_O_Storage_N[i], v_settle_N[i])
        Settling_P_S[i] = TP_MBFR.Sett_P(TP_Lake_S[i], DIP_Lake_S[i], Lake_O_A_S[i], Lake_O_Storage_S[i], v_settle_S[i])
        
        P_diff_M_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_M, DIP_pore_M_N[i], DIP_Lake_N[i], Θ_M, TP_Variables.A_Mud_N,Lake_O_Storage_N[i])
        P_diff_S_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_S, DIP_pore_S_N[i], DIP_Lake_N[i], Θ_S, TP_Variables.A_Sand_N,Lake_O_Storage_N[i])
        P_diff_R_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_R, DIP_pore_R_N[i], DIP_Lake_N[i], Θ_R, TP_Variables.A_Rock_N,Lake_O_Storage_N[i])
        P_diff_P_N[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_P, DIP_pore_P_N[i], DIP_Lake_N[i], Θ_P, TP_Variables.A_Peat_N,Lake_O_Storage_N[i])
        P_diff_M_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_M, DIP_pore_M_S[i], DIP_Lake_S[i], Θ_M, TP_Variables.A_Mud_S,Lake_O_Storage_S[i])
        P_diff_S_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_S, DIP_pore_S_S[i], DIP_Lake_S[i], Θ_S, TP_Variables.A_Sand_S,Lake_O_Storage_S[i])
        P_diff_R_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_R, DIP_pore_R_S[i], DIP_Lake_S[i], Θ_R, TP_Variables.A_Rock_S,Lake_O_Storage_S[i])
        P_diff_P_S[i] = TP_MBFR.Diff_P(TP_Variables.v_diff_P, DIP_pore_P_S[i], DIP_Lake_S[i], Θ_P, TP_Variables.A_Peat_S,Lake_O_Storage_S[i])
                
        # TP_N_to_S[i] = TP_MBFR.P_N_to_S(Q_N2S[i], TP_Lake_N[i], Lake_O_Storage_N[i])
        # TP_Out[i] = TP_MBFR.P_Out(Q_O_M[i], TP_Lake_S[i], Lake_O_Storage_S[i])
        
        TP_Lake_N[i+1] = TP_MBFR.TP_Lake_N(L_ext_M[i],Atm_Dep_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_N[i],DIP_pore_S_N[i],DIP_pore_R_N[i],DIP_pore_P_N[i],DIP_Lake_N[i],Q_N2S[i],Lake_O_A_N[i],TP_Lake_N[i],Lake_O_Storage_N[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_N[i]) + (Sed_Resusp_M_N[i]+Sed_Resusp_S_N[i]+Sed_Resusp_R_N[i]+Sed_Resusp_P_N[i]) if TP_MBFR.TP_Lake_N(L_ext_M[i],Atm_Dep_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_N[i],DIP_pore_S_N[i],DIP_pore_R_N[i],DIP_pore_P_N[i],DIP_Lake_N[i],Q_N2S[i],Lake_O_A_N[i],TP_Lake_N[i],Lake_O_Storage_N[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_N[i])+ (Sed_Resusp_M_N[i]+Sed_Resusp_S_N[i]+Sed_Resusp_R_N[i]+Sed_Resusp_P_N[i]) > 0 else 0
        TP_Lake_S[i+1] = TP_MBFR.TP_Lake_S(Atm_Dep_S[i],Q_N2S[i],TP_Lake_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_S[i],DIP_pore_S_S[i],DIP_pore_R_S[i],DIP_pore_P_S[i],DIP_Lake_S[i],Q_O_M[i],Lake_O_A_S[i],TP_Lake_S[i],Lake_O_Storage_S[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_S[i]) + (Sed_Resusp_M_S[i]+Sed_Resusp_S_S[i]+Sed_Resusp_R_S[i]+Sed_Resusp_P_S[i]) if TP_MBFR.TP_Lake_S(Atm_Dep_S[i],Q_N2S[i],TP_Lake_N[i],Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_S[i],DIP_pore_S_S[i],DIP_pore_R_S[i],DIP_pore_P_S[i],DIP_Lake_S[i],Q_O_M[i],Lake_O_A_S[i],TP_Lake_S[i],Lake_O_Storage_S[i],TP_Variables.v_diff_M,TP_Variables.v_diff_S,TP_Variables.v_diff_R,TP_Variables.v_diff_P,v_settle_S[i])+ (Sed_Resusp_M_S[i]+Sed_Resusp_S_S[i]+Sed_Resusp_R_S[i]+Sed_Resusp_P_S[i]) > 0 else 0
        TP_Lake_Mean[i+1] = ((TP_Lake_N[i+1] + TP_Lake_S[i+1])/2)
        
        P_Load_Cal[i] = S77REG[i]*0.028316847*3600*24*TP_Lake_S[i] #mg/d P
        P_Load_StL[i] = S308RG[i]*0.028316847*3600*24*TP_Lake_S[i] #mg/d P
        P_Load_South[i] = TotRegSo[i]*1233.48*TP_Lake_S[i] #mg/d P
    
    Output_df = pd.DataFrame(date_rng_2, columns=['Date']) #1/1/2008-12/31/2018

    Output_df['Stage_LO'] = Stage_LO[2:]
    Output_df['S308_Q'] = S308RG[2:]
    Output_df['S77_Q'] = S77REG[2:]
    Output_df['Storage'] = Storage[2:]
    Output_df['Cut_back'] = Cut_back[2:]
    Output_df['P_Lake'] = TP_Lake_Mean
    # Output_df['P_Lake_N'] = TP_Lake_N
    # Output_df['P_Lake_S'] = TP_Lake_S
    # Output_df['DIP_pore_M_N'] = DIP_pore_M_N
    # Output_df['Q_N2S'] = Q_N2S
    # Output_df['Lake_O_A_N'] = Lake_O_A_N
    # Output_df['Lake_O_Storage_N'] = Lake_O_Storage_N
    # Output_df['Sed_Resusp_M_N'] = Sed_Resusp_M_N
    Output_df['P_Load_Cal'] = P_Load_Cal/1E9 #tons
    Output_df['P_Load_StL'] = P_Load_StL/1E9 #tons
    Output_df['P_Load_South'] = P_Load_South/1E9 #tons

    return(Output_df)
Exported_File = LOONE_HydNut()
Exported_File.drop(index=Exported_File.index[-1],axis=0,inplace=True)
Exported_File.drop(index=Exported_File.index[-1],axis=0,inplace=True)
Exported_File['Stage_LO'] = Exported_File['Stage_LO'].astype(float)
Exported_File['Storage']=Exported_File['Storage'].astype(float)
Exported_File['S308_Q'] = Exported_File['S308_Q'].astype(float)
Exported_File['S77_Q'] = Exported_File['S77_Q'].astype(float)
Exported_File['Cut_back']=Exported_File['Cut_back'].astype(float)
Exported_File['P_Lake']=pd.to_numeric(Exported_File['P_Lake'])
# Exported_File['P_Lake_N']=pd.to_numeric(Exported_File['P_Lake_N'])
# Exported_File['P_Lake_S']=pd.to_numeric(Exported_File['P_Lake_S'])
# Exported_File['DIP_pore_M_N']=pd.to_numeric(Exported_File['DIP_pore_M_N'])
# Exported_File['Q_N2S']=pd.to_numeric(Exported_File['Q_N2S'])
# Exported_File['Lake_O_A_N']=pd.to_numeric(Exported_File['Lake_O_A_N'])
# Exported_File['Lake_O_Storage_N']=pd.to_numeric(Exported_File['Lake_O_Storage_N'])
# Exported_File['Sed_Resusp_M_N']=pd.to_numeric(Exported_File['Sed_Resusp_M_N'])
Exported_File['P_Load_Cal']=pd.to_numeric(Exported_File['P_Load_Cal'])
Exported_File['P_Load_StL']=pd.to_numeric(Exported_File['P_Load_StL'])
Exported_File['P_Load_South']=pd.to_numeric(Exported_File['P_Load_South'])

Exported_File = Exported_File.set_index('Date')
Exported_File.index = pd.to_datetime(Exported_File.index, unit = 'ns')
Exported_File_Mean = Exported_File.resample('M').mean()
Exported_File_Sum = Exported_File.resample('M').sum()

# Exported_File.to_csv('./Outputs/Daily.csv')
Exported_File_Mean.to_csv('./Outputs/Exported_File_Opt_0809_Mean_%s.csv'%Pre_defined_Variables.Schedule)
Exported_File_Sum.to_csv('./Outputs/Exported_File_Opt_0809_Sum_%s.csv'%Pre_defined_Variables.Schedule)

