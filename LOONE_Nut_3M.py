# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 18:44:37 2021

@author: osama
"""

#This Script incorporates the Comprehensive LOONE Model!
import os
import pandas as pd
from datetime import datetime
import numpy as np
from Model_Config import Model_Config
Working_Path = Model_Config.Working_Path
os.chdir('%s'%Working_Path) 
from Pre_defined_Variables import Pre_defined_Variables 
from Stg_Sto_Ar import Stg_Sto_Ar 
# from Data import Data
from TP_Variables_Regions import TP_Variables
import TP_Mass_Balance_Functions_Regions as TP_MBFR



def LOONE_Nut_NS(LOONE_Q_Outputs): 
    print("LOONE Nut Module is Running!")
    # Based on the defined Start and End year, month, and day on the Pre_defined_Variables File, Startdate and enddate are defined. 
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    startdate = datetime(year, month, day).date() 
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    begdateCS = datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.enddate_entry)
    enddate = datetime(year, month, day).date()
        
    date_rng_0 = pd.date_range(start = startdate, end = enddate, freq ='D')
    Load_ext = pd.read_csv('./Data/%s/ts_data/LO_External_Loadings_%s_3M.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    Q_in = pd.read_csv('./Data/%s/ts_data/LO_Inflows_BK_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    
   # Observed S77 S308 South
   # Outflows_Obs = pd.read_csv('./Data/%s/ts_data/Flow_df_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
   # S77_Q = Outflows_Obs['S77_Out']
   # S308_Q = Outflows_Obs['S308_Out']
   # TotRegSo = Outflows_Obs[['S351_Out','S354_Out','S352_Out','L8_Out']].sum(axis=1)/1233.48    #m3/day to acft
   # TotRegSo = Outflows_Obs['South_acft'] #acft/day

    #Obs Stage (ft) Storage (acft) 
    # Stage_Storage = pd.read_csv('./Data/%s/ts_data/Average_LO_Sto_Stg_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    # Stage_LO = Stage_Storage['Stage_ft'].astype(float)
    # Storage = Stage_Storage['Storage_acft'].astype(float)
    
    # Simulated Q
    S77_Q = LOONE_Q_Outputs['S77_Q'] * 0.0283168 * 86400 #cfs to cubic meters per day
    S308_Q = LOONE_Q_Outputs['S308_Q'] * 0.0283168 * 86400 #cfs to cubic meters per day
    TotRegSo = LOONE_Q_Outputs['TotRegSo'] #acft/day
    TotRegEW = LOONE_Q_Outputs['TotRegEW'] #acft/day
    
    # Simulated Stage (ft) Storage (acft)
    Stage_LO =  LOONE_Q_Outputs['Stage'] #ft
    Storage = LOONE_Q_Outputs['Storage'] #acft


    S65_P_Conc = pd.read_csv('./Data/%s/ts_data/S65E_TP_Conc_3M.csv'%Pre_defined_Variables.Schedule)
    S65_P = S65_P_Conc['TP_mg/m3'].astype(float)
    n_rows = len(date_rng_0)
    Storage_dev_data = pd.read_csv('./Data/%s/ts_data/Storage_Dev_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    Storage_dev = Storage_dev_data['DS_dev'] 
    L_ext = Load_ext['TP_Loads_In_mg'] #mg
    # LO_SA_Data = pd.read_csv('./Data/LORS20082023/ts_data_3M/LO_SA_20082018.csv')
    # LO_Area = LO_SA_Data['SA_acre'].astype(float)
    # RF_data = pd.read_csv('./Data/LORS20082023/ts_data_3M/RF_Volume_LORS2008.csv')
    # RF_Vol = RF_data['RF_acft'].astype(float)
    Atm_Dep_N = TP_Variables.N_Per * Load_ext['Atm_Loading_mg']
    Atm_Dep_S = TP_Variables.S_Per * Load_ext['Atm_Loading_mg']
    # C_rain = 10.417 #TP Rainfall Concentration (µg P L-1 = mg P /m3) 
    # L_drdep = 0.0385 # mg P / m2 / day 
    # Atm_Dep_N = TP_Variables.N_Per * (C_rain*RF_Vol*1233.48 + L_drdep*LO_Area*4046.85642)
    # Atm_Dep_S = TP_Variables.S_Per * (C_rain*RF_Vol*1233.48 + L_drdep*LO_Area*4046.85642)
    # Atm_Dep_N = TP_Variables.N_Per*(18/365)*LO_Area*4046.85642 #Based on data presented by Curtis Pollman, the Lake Okeechobee Technical Advisory Committee (2000) recommended that 18 mgP/m2-yr is an appropriate atmospheric loading of phosphorus over the open lake. 
    # Atm_Dep_S = TP_Variables.S_Per*(18/365)*LO_Area*4046.85642
    #Read Shear Stress driven by Wind Speed
    Wind_ShearStr = pd.read_csv('./Data/%s/ts_data/WindShearStress_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    W_SS = Wind_ShearStr['ShearStress'] #Dyne/cm2
    Current_Stress =  pd.read_csv('./Data/%s/ts_data/Current_ShearStress.csv'%Pre_defined_Variables.Schedule)
    Current_SS = Current_Stress['Current_Stress']
    Bottom_Stress = W_SS + Current_SS
    nu_ts = pd.read_csv('./Data/%s/ts_data/nu_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    LO_BL = 0.5 # m (Bed Elevation of LO)
    #Temp 
    Temp_data = pd.read_csv('./Data/%s/ts_data/Temp_Avg_%s.csv'%(Pre_defined_Variables.Schedule,Pre_defined_Variables.Schedule))
    Temp = Temp_data['Mean_T'].astype(float)

    # LO_WD = pd.to_numeric(Stage_Storage['Stage_m'])-LO_BL
    g = 9.8 #m/s2 gravitational acceleration
    Cal_Res = pd.read_csv('./Data/%s/nondominated_Sol_var.csv'%Pre_defined_Variables.Schedule)
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
    # Q_O = Q_Out['Total_Outflows_acft'].astype(float) * 4046.85642 * 0.305 #m3/d
    # Q_O = Q_Out['Outflows_cmd'].astype(float)
    # Q_O = (TotRegEW + TotRegSo) * 4046.85642 * 0.305 #m3/d
    
    ## Observed S77 S308 South
    Q_O = S77_Q + S308_Q + TotRegSo * 1233.48  #cmd
    Q_O_M = np.zeros(n_rows,dtype = object)
    Q_N2S = np.zeros(n_rows,dtype = object)

    P_Load_Cal = np.zeros(n_rows,dtype = object)
    P_Load_StL = np.zeros(n_rows,dtype = object)
    P_Load_South = np.zeros(n_rows,dtype = object)
    
    v_settle_N_c = np.zeros(n_rows,dtype = object)
    v_settle_N_s = np.zeros(n_rows,dtype = object)
    v_settle_N = np.zeros(n_rows,dtype = object)
    v_settle_S_c = np.zeros(n_rows,dtype = object)
    v_settle_S_s = np.zeros(n_rows,dtype = object)
    v_settle_S = np.zeros(n_rows,dtype = object)
    # Ext_Load_Rate_N = np.zeros(n_rows,dtype = object)
    # Ext_Load_Rate_S = np.zeros(n_rows,dtype = object)
    # Load_Out_N = np.zeros(n_rows,dtype = object)
    # Load_Out_S = np.zeros(n_rows,dtype = object)
    ##Initial Values##
    #S.A. is calculated based on the Lake's previous time step Stage, but for the S.A. at i=0 I used same time step Stage!
    StartStorage = Stg_Sto_Ar.stg2sto(Pre_defined_Variables.startstage,0)
    Stage2ar[0] = Stg_Sto_Ar.stg2ar(Stage_LO[0],0)
    Stage2ar[1] = Stg_Sto_Ar.stg2ar(Stage_LO[1],0)
    Storage[0] = StartStorage #ac-ft
    Storage[1] = Stg_Sto_Ar.stg2sto(Stage_LO[1],0) #ac-ft
    #TP_MassBalanceModel Initial Values.
    TP_Lake_N[0] = 300 #mg/m3
    TP_Lake_S[0] = 340 #mg/m3
    TP_Lake_Mean[0] = (TP_Lake_N[0] + TP_Lake_S[0])/2
    Γ_M_N[0] = 25 #mg/kg 
    Γ_S_N[0] = 25 #mg/kg 
    Γ_R_N[0] = 25 #mg/kg 
    Γ_P_N[0] = 25 #mg/kg 
    Γ_M_S[0] = 25 #mg/kg
    Γ_S_S[0] = 25 #mg/kg 
    Γ_R_S[0] = 25 #mg/kg 
    Γ_P_S[0] = 25 #mg/kg  
    
    # DIP_pore_M_N[0] = 3000#760 #mg/m3 
    # DIP_pore_S_N[0] = 1500#205 #mg/m3 
    # DIP_pore_R_N[0] = 1500#205 #mg/m3 
    # DIP_pore_P_N[0] = 1000#160 #mg/m3 
    # DIP_pore_M_S[0] = 3000#760 #mg/m3 
    # DIP_pore_S_S[0] = 1500#205 #mg/m3 
    # DIP_pore_R_S[0] = 1500#205 #mg/m3
    # DIP_pore_P_S[0] = 1000#160 #mg/m3 
    # P_sed_M_N[0] = 3000 #mg/kg 
    # P_sed_S_N[0] = 1000 #mg/kg  
    # P_sed_R_N[0] = 1000 #mg/kg 
    # P_sed_P_N[0] = 1600#mg/kg 
    # P_sed_M_S[0] = 3000 #mg/kg 
    # P_sed_S_S[0] = 1000 #mg/kg 
    # P_sed_R_S[0] = 1000 #mg/kg 
    # P_sed_P_S[0] = 1600 #mg/kg 
   
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
    
    # DIP_pore_M_N[0] = 15#760 #mg/m3 
    # DIP_pore_S_N[0] = 17#205 #mg/m3 
    # DIP_pore_R_N[0] = 16.5#205 #mg/m3 
    # DIP_pore_P_N[0] = 111#160 #mg/m3 
    # DIP_pore_M_S[0] = 7.5#760 #mg/m3 
    # DIP_pore_S_S[0] = 20#205 #mg/m3 
    # DIP_pore_R_S[0] = 40#205 #mg/m3
    # DIP_pore_P_S[0] = 14#160 #mg/m3 
    # P_sed_M_N[0] = 910 #mg/kg 
    # P_sed_S_N[0] = 423 #mg/kg  
    # P_sed_R_N[0] = 362 #mg/kg 
    # P_sed_P_N[0] = 253 #mg/kg 
    # P_sed_M_S[0] = 1026 #mg/kg 
    # P_sed_S_S[0] = 120 #mg/kg 
    # P_sed_R_S[0] = 144 #mg/kg 
    # P_sed_P_S[0] = 132 #mg/kg 
    
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
    P_Lake_df = pd.DataFrame(date_rng_0, columns=['Date']) #1/1/2008-12/31/2018

    for i in range(n_rows):  
        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48 #m3/d
            Q_O_M[i] = Q_O[i]
            L_ext_M[i] = L_ext[i] + Q_I_M[i] * S65_P[i]            
        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48 #m3/d
            Q_I_M[i] = Q_I[i]
            L_ext_M[i] = L_ext[i]
        # Q_N2S[i] = (Q_I_M[i]*TP_Variables.N_Per + Q_O_M[i]*TP_Variables.S_Per)   
        Q_N2S[i] = (Q_I_M[i]*1 + Q_O_M[i]*0)   

    for i in range(n_rows-2):
        
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
       
        # Sed_Resusp_M_N[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_M_N[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        # Sed_Resusp_S_N[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_S_N[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        # Sed_Resusp_R_N[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_R_N[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        # Sed_Resusp_P_N[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_P_N[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        # Sed_Resusp_M_S[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_M_S[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        # Sed_Resusp_S_S[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_S_S[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        # Sed_Resusp_R_S[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_R_S[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        # Sed_Resusp_P_S[i] = ((E_0/Td**E_1)*((Bottom_Stress[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10/LO_WD[i]*P_sed_P_S[i] if Bottom_Stress[i] > Crtcl_ShStr else 0
        
        # Sed_Resusp_M_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_M_N[i]*TP_Variables.A_Mud_N/Lake_O_Storage_N[i] if W_SS[i] > Crtcl_ShStr else 0
        # Sed_Resusp_S_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_S_N[i]*TP_Variables.A_Sand_N/Lake_O_Storage_N[i] if W_SS[i] > Crtcl_ShStr else 0
        # Sed_Resusp_R_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_R_N[i]*TP_Variables.A_Rock_N/Lake_O_Storage_N[i] if W_SS[i] > Crtcl_ShStr else 0
        # Sed_Resusp_P_N[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_P_N[i]*TP_Variables.A_Peat_N/Lake_O_Storage_N[i] if W_SS[i] > Crtcl_ShStr else 0
        # Sed_Resusp_M_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_M_S[i]*TP_Variables.A_Mud_S/Lake_O_Storage_S[i] if W_SS[i] > Crtcl_ShStr else 0
        # Sed_Resusp_S_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_S_S[i]*TP_Variables.A_Sand_S/Lake_O_Storage_S[i] if W_SS[i] > Crtcl_ShStr else 0
        # Sed_Resusp_R_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_R_S[i]*TP_Variables.A_Rock_S/Lake_O_Storage_S[i] if W_SS[i] > Crtcl_ShStr else 0
        # Sed_Resusp_P_S[i] = ((E_0/Td**E_1)*((W_SS[i]-Crtcl_ShStr)/Crtcl_ShStr)**E_2)*10*P_sed_P_S[i]*TP_Variables.A_Peat_S/Lake_O_Storage_S[i] if W_SS[i] > Crtcl_ShStr else 0

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
        
        # if P_Lake_df['Date'].iloc[i] == datetime(2030, 1, 1):
        #     P_sed_M_N[0] = 910 #mg/kg 
        #     P_sed_S_N[0] = 423 #mg/kg  
        #     P_sed_R_N[0] = 362 #mg/kg 
        #     P_sed_P_N[0] = 253 #mg/kg 
        #     P_sed_M_S[0] = 1026 #mg/kg 
        #     P_sed_S_S[0] = 120 #mg/kg 
        #     P_sed_R_S[0] = 144 #mg/kg 
        #     P_sed_P_S[0] = 132 #mg/kg 
        #     DIP_pore_M_N[0] = 15#760 #mg/m3 
        #     DIP_pore_S_N[0] = 17#205 #mg/m3 
        #     DIP_pore_R_N[0] = 16.5#205 #mg/m3 
        #     DIP_pore_P_N[0] = 111#160 #mg/m3 
        #     DIP_pore_M_S[0] = 7.5#760 #mg/m3 
        #     DIP_pore_S_S[0] = 20#205 #mg/m3 
        #     DIP_pore_R_S[0] = 40#205 #mg/m3
        #     DIP_pore_P_S[0] = 14#160 #mg/m3 
            
        # else:
            # P_sed_M_N[i+1] = TP_MBFR.P_sed(Lake_O_A_M_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.K_decomp_M,v_settle_N[i]) - Sed_Resusp_M_N[i]*Lake_O_Storage_N[i]/Mass_sed_M_N if TP_MBFR.P_sed(Lake_O_A_M_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_M_N[i],P_sed_M_N[i],Mass_sed_M_N,TP_Variables.K_decomp_M,v_settle_N[i]) - Sed_Resusp_M_N[i]*Lake_O_Storage_N[i]/Mass_sed_M_N > 0 else 0
            # P_sed_S_N[i+1] = TP_MBFR.P_sed(Lake_O_A_S_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.K_decomp_S,v_settle_N[i]) - Sed_Resusp_S_N[i]*Lake_O_Storage_N[i]/Mass_sed_S_N if TP_MBFR.P_sed(Lake_O_A_S_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_S_N[i],P_sed_S_N[i],Mass_sed_S_N,TP_Variables.K_decomp_S,v_settle_N[i]) - Sed_Resusp_S_N[i]*Lake_O_Storage_N[i]/Mass_sed_S_N > 0 else 0
            # P_sed_R_N[i+1] = TP_MBFR.P_sed(Lake_O_A_R_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.K_decomp_R,v_settle_N[i]) - Sed_Resusp_R_N[i]*Lake_O_Storage_N[i]/Mass_sed_R_N if TP_MBFR.P_sed(Lake_O_A_R_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_R_N[i],P_sed_R_N[i],Mass_sed_R_N,TP_Variables.K_decomp_R,v_settle_N[i]) - Sed_Resusp_R_N[i]*Lake_O_Storage_N[i]/Mass_sed_R_N > 0 else 0
            # P_sed_P_N[i+1] = TP_MBFR.P_sed(Lake_O_A_P_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.K_decomp_P,v_settle_N[i]) - Sed_Resusp_P_N[i]*Lake_O_Storage_N[i]/Mass_sed_P_N if TP_MBFR.P_sed(Lake_O_A_P_N[i],TP_Lake_N[i],DIP_Lake_N[i],J_sedburial_P_N[i],P_sed_P_N[i],Mass_sed_P_N,TP_Variables.K_decomp_P,v_settle_N[i]) - Sed_Resusp_P_N[i]*Lake_O_Storage_N[i]/Mass_sed_P_N > 0 else 0
            # P_sed_M_S[i+1] = TP_MBFR.P_sed(Lake_O_A_M_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.K_decomp_M,v_settle_S[i]) - Sed_Resusp_M_S[i]*Lake_O_Storage_S[i]/Mass_sed_M_S if TP_MBFR.P_sed(Lake_O_A_M_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_M_S[i],P_sed_M_S[i],Mass_sed_M_S,TP_Variables.K_decomp_M,v_settle_S[i]) - Sed_Resusp_M_S[i]*Lake_O_Storage_S[i]/Mass_sed_M_S > 0 else 0
            # P_sed_S_S[i+1] = TP_MBFR.P_sed(Lake_O_A_S_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.K_decomp_S,v_settle_S[i]) - Sed_Resusp_S_S[i]*Lake_O_Storage_S[i]/Mass_sed_S_S if TP_MBFR.P_sed(Lake_O_A_S_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_S_S[i],P_sed_S_S[i],Mass_sed_S_S,TP_Variables.K_decomp_S,v_settle_S[i]) - Sed_Resusp_S_S[i]*Lake_O_Storage_S[i]/Mass_sed_S_S > 0 else 0
            # P_sed_R_S[i+1] = TP_MBFR.P_sed(Lake_O_A_R_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.K_decomp_R,v_settle_S[i]) - Sed_Resusp_R_S[i]*Lake_O_Storage_S[i]/Mass_sed_R_S if TP_MBFR.P_sed(Lake_O_A_R_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_R_S[i],P_sed_R_S[i],Mass_sed_R_S,TP_Variables.K_decomp_R,v_settle_S[i]) - Sed_Resusp_R_S[i]*Lake_O_Storage_S[i]/Mass_sed_R_S > 0 else 0
            # P_sed_P_S[i+1] = TP_MBFR.P_sed(Lake_O_A_P_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.K_decomp_P,v_settle_S[i]) - Sed_Resusp_P_S[i]*Lake_O_Storage_S[i]/Mass_sed_P_S if TP_MBFR.P_sed(Lake_O_A_P_S[i],TP_Lake_S[i],DIP_Lake_S[i],J_sedburial_P_S[i],P_sed_P_S[i],Mass_sed_P_S,TP_Variables.K_decomp_P,v_settle_S[i]) - Sed_Resusp_P_S[i]*Lake_O_Storage_S[i]/Mass_sed_P_S > 0 else 0
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
        # Suspended_Sed_N[i] = ((Sed_Resusp_M_N[i]+Sed_Resusp_S_N[i]+Sed_Resusp_R_N[i]+Sed_Resusp_P_N[i])/Lake_O_Storage_N[i])
        # Suspended_Sed_S[i] = ((Sed_Resusp_M_S[i]+Sed_Resusp_S_S[i]+Sed_Resusp_R_S[i]+Sed_Resusp_P_S[i])/Lake_O_Storage_S[i])
        # Ext_Load_Rate_N[i] = (L_ext_M[i]+Atm_Dep_N[i])/Lake_O_Storage_N[i]
        # Ext_Load_Rate_S[i] = (Atm_Dep_S[i]+Q_N2S[i]*TP_Lake_N[i])/Lake_O_Storage_S[i]
        # Load_Out_N[i] = (Q_N2S[i]*TP_Lake_N[i])/Lake_O_Storage_N[i]
        # Load_Out_S[i] = (Q_O_M[i]*TP_Lake_S[i])/Lake_O_Storage_S[i]
        # P_Load_Cal[i] = S77_Q[i]*0.028316847*3600*24*TP_Lake_S[i] #mg
        # P_Load_StL[i] = S308_Q[i]*0.028316847*3600*24*TP_Lake_S[i] #mg
        # P_Load_South[i] = TotRegSo[i]*1233.48*TP_Lake_S[i]
        #Obs S77 S308 South
        P_Load_Cal[i] = S77_Q[i]*TP_Lake_S[i] #mg
        P_Load_StL[i] = S308_Q[i]*TP_Lake_S[i] #mg
        P_Load_South[i] = TotRegSo[i]*1233.48*TP_Lake_S[i]

    
    P_Loads_df = pd.DataFrame(date_rng_0, columns=['Date']) #1/1/2008-12/31/2018
    P_Lake_df = pd.DataFrame(date_rng_0, columns=['Date']) #1/1/2008-12/31/2018

    P_Loads_df['P_Load_Cal'] = pd.to_numeric(P_Load_Cal)/1E9 #tons
    P_Loads_df['P_Load_StL'] = pd.to_numeric(P_Load_StL)/1E9 #tons
    P_Loads_df['P_Load_South'] = pd.to_numeric(P_Load_South)/1E9 #tons
    P_Lake_df['P_Lake'] = pd.to_numeric(TP_Lake_Mean)
    P_Lake_df['TP_Lake_S'] = pd.to_numeric(TP_Lake_S)
    P_Lake_df['TP_Lake_N'] = pd.to_numeric(TP_Lake_N)
    P_Lake_df['Water Temp'] = pd.to_numeric(Temp) # C

    # P_Lake_df['Suspended_Sed_N'] = pd.to_numeric(Suspended_Sed_N)
    # P_Lake_df['Suspended_Sed_S'] = pd.to_numeric(Suspended_Sed_S)
    # P_Lake_df['Settling_P_N'] = pd.to_numeric(Settling_P_N)
    # P_Lake_df['Settling_P_S'] = pd.to_numeric(Settling_P_S)

    P_Loads_df = P_Loads_df.set_index('Date')
    P_Loads_df.index = pd.to_datetime(P_Loads_df.index, unit = 'ns')
    P_Loads_M = P_Loads_df.resample('M').sum()
    P_Loads_M = P_Loads_M.reset_index()
    P_Lake_df = P_Lake_df.set_index('Date')
    P_Lake_M = P_Lake_df.resample('M').mean()
    P_Lake_df = P_Lake_df.reset_index()

    # return(P_Loads_M)
    
    
    # Smr_Mnth_StL = []
    # Smr_Mnth_Cal = []
    # for i in range(len(P_Loads_M.index)):
    #     if P_Loads_M['Date'].iloc[i].month in [5,6,7,8,9,10]:
    #         Smr_Mnth_StL.append(P_Loads_M['P_Load_StL'].iloc[i])
    #         Smr_Mnth_Cal.append(P_Loads_M['P_Load_Cal'].iloc[i])
    # Smr_Mnth_StL_arr = np.asarray(Smr_Mnth_StL)
    # Smr_Mnth_Cal_arr = np.asarray(Smr_Mnth_Cal)
    # if Model_Config.Sim_type == 0 or Model_Config.Sim_type == 1 :
    #     return[P_Loads_M,P_Lake_M,Smr_Mnth_StL_arr,Smr_Mnth_Cal_arr]
    # else:
    #     return[Smr_Mnth_StL_arr,Smr_Mnth_Cal_arr,P_Lake_df]
    
    Algae_Opt_Mnth_StL = []
    Algae_Opt_Mnth_Cal = []
    for i in range(len(P_Loads_M.index)):
        if P_Lake_M['Water Temp'].iloc[i] >= 25: 
            Algae_Opt_Mnth_StL.append(P_Loads_M['P_Load_StL'].iloc[i])
            Algae_Opt_Mnth_Cal.append(P_Loads_M['P_Load_Cal'].iloc[i])
    Algae_Opt_Mnth_StL_arr = np.asarray(Algae_Opt_Mnth_StL)
    Algae_Opt_Mnth_Cal_arr = np.asarray(Algae_Opt_Mnth_Cal)
    if Model_Config.Sim_type == 0 or Model_Config.Sim_type == 1 :
        return[P_Loads_M,P_Lake_M,Algae_Opt_Mnth_StL_arr,Algae_Opt_Mnth_Cal_arr]
    else:
        return[Algae_Opt_Mnth_StL_arr,Algae_Opt_Mnth_Cal_arr,P_Lake_df]