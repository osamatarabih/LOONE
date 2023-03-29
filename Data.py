# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 01:19:24 2022

@author: osama
"""
# READ All Required Data Inputs
import pandas as pd
import os
from Model_Config import Model_Config 
Working_Path = Model_Config.Working_Path
os.chdir('%s'%Working_Path) 
from Pre_defined_Variables import Pre_defined_Variables 

class Data:
    #Read SFWMM Daily Output File
    SFWMM_Daily_Outputs = pd.read_csv('./Data/SFWMM_Daily_Outputs_%s.csv'%Pre_defined_Variables.Schedule)
    #read the Water Shortrage Management and Regulation Schedule Break points csv file
    #FIXME [Just Make Sure you are using the WSM and Operation Zone Lake Eleveation Data Associated with the identified Operation Rules] 
    WSMs_RSBKs = pd.read_csv('./Data/WSMs_RSBPs_%s.csv'%Pre_defined_Variables.Schedule)
    #Read the weekly demand input file
    Weekly_dmd = pd.read_csv('./Data/LOSA_wkly_dmd_%s.csv'%Pre_defined_Variables.Schedule)
    #read Weekly tributary condition (NetRF, S65E runoff, Palmer index, and Netinflow) data for the study period
    Wkly_Trib_Cond = pd.read_csv('./Data/Trib_cond_wkly_data_%s.csv'%Pre_defined_Variables.Schedule)
    #Read Seasonal LONINO data (the latest values for use with LORS simulations) defined as Option 4 in the LOOPS Model (TC-LONINO Tab).
    LONINO_Seas_data = pd.read_csv('./Data/Seasonal_LONINO_%s.csv'%Pre_defined_Variables.Schedule)
    LONINO_Mult_Seas_data = pd.read_csv('./Data/Multi_Seasonal_LONINO_%s.csv'%Pre_defined_Variables.Schedule)
    #read Netinflow data(ac-ft/day)
    NetInf_Input = pd.read_csv('./Data/NetFlows_acft_%s.csv'%Pre_defined_Variables.Schedule)
    #read actual daily water demand (output of SFWMM)
    SFWMM_W_dmd = pd.read_csv('./Data/Water_dmd_%s.csv'%Pre_defined_Variables.Schedule)
    #read Rainfall Volume data
    RF_Vol = pd.read_csv('./Data/RF_Volume_%s.csv'%Pre_defined_Variables.Schedule)
    #read ET Vol data
    ET_Vol = pd.read_csv('./Data/ETVol_%s.csv'%Pre_defined_Variables.Schedule)
    #Read the C44 Runoff data which is output of SFWMM simulation.
    C44_Runoff = pd.read_csv('./Data/C44RO_%s.csv'%Pre_defined_Variables.Schedule)
    
    #Read the mean monthly sum Basin runoffs file (i.e., sum of flow volume for each month of the entire period)
    Sum_Basin_RO = pd.read_csv('./Data/Basin_RO_inputs_%s.csv'%Pre_defined_Variables.Schedule)
    #Read Mean Monthly basin runoffs (cfs)
    C43RO = pd.read_csv('./Data/C43RO_Monthly_%s.csv'%Pre_defined_Variables.Schedule)
    C44RO = pd.read_csv('./Data/C44RO_Monthly_%s.csv'%Pre_defined_Variables.Schedule)
    SLTRIB = pd.read_csv('./Data/SLTRIB_Monthly_%s.csv'%Pre_defined_Variables.Schedule)
    #Read S80 and S77 Regulatory release rates for LORS2008 (Inputs from ActiveSchedule Tab in LOOPS Model!)
    S77_RegRelRates = pd.read_csv('./Data/S77_RegRelRates_%s.csv'%Pre_defined_Variables.Schedule)
    S80_RegRelRates = pd.read_csv('./Data/S80_RegRelRates_%s.csv'%Pre_defined_Variables.Schedule)
    #Read the daily C43Runoff
    C43RO_Daily = pd.read_csv('./Data/C43RO_%s.csv'%Pre_defined_Variables.Schedule)
    #FIXME: I will use the CE and SLE turns out of the LOOPS model as inputs here (to figure out how to calculate it later in this model!)
    CE_SLE_turns = pd.read_csv('./Data/CE_SLE_turns_inputs_%s.csv'%Pre_defined_Variables.Schedule)
    #Read the pulses input data (Tab Pulses in the LOOPS Spreadsheet Model!)
    Pulses =  pd.read_csv('./Data/Pulses_Inputs_%s.csv'%Pre_defined_Variables.Schedule)
    #Note that I entered values of (-9999) in days of the year that do not include data (i.e. June 2nd to September 30th)
    Targ_Stg_June_1st = pd.read_csv('./Data/Chance of June 1st Lake stage falling below 11.0ft_%s.csv'%Pre_defined_Variables.Schedule, parse_dates = ['Date'])
    Targ_Stg_May_1st = pd.read_csv('./Data/Chance of May 1st Lake stage falling below 11.0ft_%s.csv'%Pre_defined_Variables.Schedule, parse_dates = ['Date'])
    #FIXME
    #The Estuary needs water needs to be updated if needed!
    #Read the "Estuary needs water" for now (to be calculated later!)
    Estuary_needs_water = pd.read_csv('./Data/Estuary_needs_water_Input_%s.csv'%Pre_defined_Variables.Schedule)
    #Read EAA_MIA_RUNOFF in cfs
    EAA_MIA_RUNOFF = pd.read_csv('./Data/EAA_MIA_RUNOFF_Inputs_%s.csv'%Pre_defined_Variables.Schedule)
    #Read calculated Storage deviation daily values 
    Stroage_dev_df = pd.read_csv('./Data/Storage_Dev_%s.csv'%Pre_defined_Variables.Schedule)