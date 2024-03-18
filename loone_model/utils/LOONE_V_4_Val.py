# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 18:44:37 2021

@author: osama
"""

# This Script incorporates the Comprehensive LOONE Model!
Working_Path = "C:/Osama_PC/LOONE/Model/LOONE_Model"
import os
import pandas as pd
from datetime import datetime
import numpy as np
from calendar import monthrange

os.chdir("%s" % Working_Path)
from data.pre_defined_variables import Pre_defined_Variables
from data.model_variables import M_var
from loone_config.Model_Config import Model_Config
import utils.lo_functions
import utils.stg_sto_ar
import utils.dec_tree_functions
from utils.wca_stages_class import WCA_Stages_Cls
import utils.additional_functions
import utils.thc_class
from data.data import Data
import utils.df_wsms
import utils.trib_hc
from data.tp_variables_regions import TP_Variables
import utils.tp_mass_balance_functions_regions as TP_MBFR


def LOONE_HydNut():
    # Based on the defined Start and End year, month, and day on the Pre_defined_Variables File, Startdate and enddate are defined.
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    startdate = datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    begdateCS = datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.enddate_entry)
    enddate = datetime(year, month, day).date()

    #############################################################################################################################
    Results_data = pd.read_csv("./Outputs/Opt_Decision_Var.csv")
    Results = Results_data["Value"]
    P_1 = Results[0]
    P_2 = Results[1]
    S77_DV = Results[2:14]
    S308_DV = Results[14:26]
    # First, I interpolated each Water Shortage Management (WSMs) and each Regulation Schedule Breakpoint Zone (D, C, B, and A).
    # Set time frame for model run such that it starts on the defined startdate but ends on 1/1/(endyear+1)
    date_rng_1 = pd.date_range(
        start=startdate, end="1/1/%d" % (Pre_defined_Variables.endyear + 1), freq="D"
    )
    # Create a data frame with a date column
    if Model_Config.Sim_type == 0 or Model_Config.Sim_type == 1:
        utils.df_wsms.WSMs()
        df_WSMs = pd.read_csv("./Data/df_WSMs.csv")
    else:
        df_WSMs = pd.read_csv("./Data/df_WSMs.csv")
    #############################################################################
    # The Following Code interpolates daily LOSA demand from weekly data for 6 differnet datasets where the user defines the LOSA demand that will be used based on a Code (1:6).
    # Set time frame for model run
    date_rng_2 = pd.date_range(start=startdate, end=enddate, freq="D")
    # Create a data frame with a date column
    Water_dmd = pd.DataFrame(date_rng_2, columns=["date"])

    N = []
    Wk = []
    # Generate a count list
    for i in Water_dmd["date"]:
        if i.month == 1 and i.day == 1:
            n = 0
        else:
            n = n + 1
        N.append(n)
    Water_dmd["count"] = N
    # Calculate the week number for all rows in the data frame
    for i in Water_dmd["count"]:
        if i > 363:
            J = 52
        else:
            J = int(i / 7) + 1
        Wk.append(J)
    Water_dmd["Week_num"] = Wk
    dd = []  # daily demand
    # Calculate daily water demand
    for i in Water_dmd["Week_num"]:
        D = ((Data.Weekly_dmd["C%s" % Pre_defined_Variables.Code].iloc[i - 1]) / 7) * (
            Pre_defined_Variables.Multiplier / 100
        )
        dd.append(D)
    Water_dmd["Daily_demand"] = dd
    ##############################################################################################
    # Determine Tributary Hydrologic Conditions
    TC_LONINO_df = utils.trib_hc.Trib_HC()
    # Determine WCA Stages
    WCA_Stages_df = WCA_Stages_Cls(TC_LONINO_df)
    # A dataframe to determine eachday's season (Months 11,12,1,2 are Season 1, Months 3,4,5 are season 2, Months 6,7 are season 3, Months 8,9,10 are season 4 )
    date_rng_5 = pd.date_range(start=startdate, end=enddate, freq="D")
    Seasons = pd.DataFrame(date_rng_5, columns=["date"])
    Seas_Count = len(Seasons.index)
    for i in range(Seas_Count):
        if Seasons["date"].iloc[i].month > 2 and Seasons["date"].iloc[i].month < 6:
            S = 2
        elif Seasons["date"].iloc[i].month > 5 and Seasons["date"].iloc[i].month < 8:
            S = 3
        elif Seasons["date"].iloc[i].month > 7 and Seasons["date"].iloc[i].month < 11:
            S = 4
        else:
            S = 1
        M_var.Daily_Seasons[i] = S
        M_var.Mon[i] = Seasons["date"].iloc[i].month
    Seasons["Season"] = M_var.Daily_Seasons
    Seasons["Month"] = M_var.Mon
    ##################################################################################################################
    # This following Script runs the main model daily simulations.
    date_rng_6 = pd.date_range(
        start="12/30/%d" % (Pre_defined_Variables.startyear - 1), end=enddate, freq="D"
    )
    LO_Model = pd.DataFrame(date_rng_6, columns=["date"])
    LO_Model["Net_Inflow"] = Data.NetInf_Input["Netflows_acft"]
    n_rows = len(LO_Model.index)
    LO_Model["LOSA_dmd_SFWMM"] = Data.SFWMM_W_dmd["LOSA_dmd"] * (
        Pre_defined_Variables.Mult_LOSA / 100
    )
    LO_Model["C44RO"] = Data.C44_Runoff["C44RO"]
    ##################################
    DecTree_df = pd.DataFrame(date_rng_5, columns=["Date"])
    DecTree_df["Zone_B_MetFcast"] = TC_LONINO_df["LONINO_Seasonal_Classes"]
    # Create a dataframe that includes Monthly Mean Basin Runoff & BaseFlow-Runoff & Runoff-Baseflow (cfs)
    date_rng_11 = pd.date_range(start=startdate, end=enddate, freq="MS")
    date_rng_11d = pd.date_range(start=startdate, end=enddate, freq="D")
    date_rng_11d.name = "Date"
    Basin_RO = pd.DataFrame(date_rng_11, columns=["date"])
    # Baseflows
    Outlet1_baseflow = Data.S77_RegRelRates["Zone_D0"].iloc[0]
    Outlet2_baseflow = Data.S80_RegRelRates["Zone_D0"].iloc[0]
    # Calculta number of months in the timeseries data.
    num_B_R = len(Basin_RO.index)
    BS_C43RO = np.zeros(num_B_R)
    BS_C44RO = np.zeros(num_B_R)
    C44RO_SLTRIB = np.zeros(num_B_R)
    C44RO_BS = np.zeros(num_B_R)
    Num_days = np.zeros(num_B_R)
    for i in range(num_B_R):
        Num_days[i] = monthrange(
            Basin_RO["date"].iloc[i].year, Basin_RO["date"].iloc[i].month
        )[
            1
        ]  # no. of days in each time step month.
        BS_C43RO[i] = max(0, (Outlet1_baseflow - Data.C43RO["C43RO"].iloc[i]))
        BS_C44RO[i] = max(0, (Outlet2_baseflow - Data.C44RO["C44RO"].iloc[i]))
        C44RO_SLTRIB[i] = BS_C44RO[i] + Data.SLTRIB["SLTRIB_cfs"].iloc[i]
        C44RO_BS[i] = (
            max(0, Data.C44RO["C44RO"].iloc[i] - Outlet2_baseflow) * Num_days[i]
        )
    Basin_RO["Ndays"] = Num_days
    Basin_RO["C43RO"] = Data.C43RO["C43RO"]
    Basin_RO["BS-C43RO"] = BS_C43RO
    Basin_RO["C44RO"] = Data.C44RO["C44RO"]
    Basin_RO["BS-C44RO"] = BS_C44RO
    Basin_RO["SLTRIB"] = Data.SLTRIB["SLTRIB_cfs"]
    Basin_RO["C44RO_SLTRIB"] = C44RO_SLTRIB
    Basin_RO["C44RO-BS"] = C44RO_BS
    LO_Model["C43RO"] = Data.C43RO_Daily["C43RO"]
    S80avgL1 = Data.Pulses["S-80_L1_%s" % Pre_defined_Variables.Schedule].mean()
    S80avgL2 = Data.Pulses["S-80_L2_%s" % Pre_defined_Variables.Schedule].mean()
    S80avgL3 = Data.Pulses["S-80_L3_%s" % Pre_defined_Variables.Schedule].mean()
    S77avgL1 = Data.Pulses["S-77_L1_%s" % Pre_defined_Variables.Schedule].mean()  # LORS
    S77avgL2 = Data.Pulses["S-77_L2_%s" % Pre_defined_Variables.Schedule].mean()  # LORS
    S77avgL3 = Data.Pulses["S-77_L3_%s" % Pre_defined_Variables.Schedule].mean()
    Basin_RO = Basin_RO.set_index(["date"])
    Basin_RO.index = pd.to_datetime(Basin_RO.index)
    Basin_RO_Daily = Basin_RO.reindex(date_rng_11d, method="ffill")
    Basin_RO = Basin_RO.reset_index()
    VLOOKUP1 = Basin_RO_Daily["BS-C44RO"]
    VLOOKUP1_c = [x for x in VLOOKUP1 if ~np.isnan(x)]
    ##################################################################################################################
    # This following script contains the logic and calculations for the proposed Lake Okeechobee Adaptive Protocol.
    AdapProt_df = pd.DataFrame(date_rng_5, columns=["date"])
    # Calculate Late Dry Season (Apr-May) logic.
    Late_Dry_Season = []
    for i in AdapProt_df["date"]:
        if i.month > 3 and i.month < 6:
            L = True
        else:
            L = False
        Late_Dry_Season.append(L)
    AdapProt_df["Late_Dry_Season"] = Late_Dry_Season
    AdapProt_df["Tributary Hydrologic Condition"] = TC_LONINO_df["Tributary_Condition"]
    # Define "Low Chance" 6/1 stg<11'
    if Pre_defined_Variables.Opt_Date_Targ_Stg == 1:
        Targ_Stg = Data.Targ_Stg_June_1st
    else:
        Targ_Stg = Data.Targ_Stg_May_1st

    Targ_Stg_df = pd.DataFrame(date_rng_5, columns=["dates"])
    for i in range(len(Targ_Stg_df)):
        M_var.V10per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            10,
            Targ_Stg,
        )
        M_var.V20per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            20,
            Targ_Stg,
        )
        M_var.V25per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            25,
            Targ_Stg,
        )
        M_var.V30per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            30,
            Targ_Stg,
        )
        M_var.V40per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            40,
            Targ_Stg,
        )
        M_var.V45per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            45,
            Targ_Stg,
        )
        M_var.V50per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            50,
            Targ_Stg,
        )
        M_var.V60per[i] = utils.additional_functions.replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            60,
            Targ_Stg,
        )

    V10per_c = [x for x in M_var.V10per if ~np.isnan(x)]
    V20per_c = [x for x in M_var.V20per if ~np.isnan(x)]
    V25per_c = [x for x in M_var.V25per if ~np.isnan(x)]
    V30per_c = [x for x in M_var.V30per if ~np.isnan(x)]
    V40per_c = [x for x in M_var.V40per if ~np.isnan(x)]
    V45per_c = [x for x in M_var.V45per if ~np.isnan(x)]
    V50per_c = [x for x in M_var.V50per if ~np.isnan(x)]
    V60per_c = [x for x in M_var.V60per if ~np.isnan(x)]
    Targ_Stg_df["10%"] = V10per_c
    Targ_Stg_df["20%"] = V20per_c
    Targ_Stg_df["25%"] = V25per_c
    Targ_Stg_df["30%"] = V30per_c
    Targ_Stg_df["40%"] = V40per_c
    Targ_Stg_df["45%"] = V45per_c
    Targ_Stg_df["50%"] = V50per_c
    Targ_Stg_df["60%"] = V60per_c

    # Outlet1_baseflow = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    Outlet1_baseflow = 450  # cfs
    VLOOKUP2 = Basin_RO_Daily["BS-C43RO"]
    VLOOKUP2_c = [x for x in VLOOKUP2 if ~np.isnan(x)]
    ####################################################################################################################
    M_var.Lake_Stage[0] = Pre_defined_Variables.begstageCS
    M_var.Lake_Stage[1] = Pre_defined_Variables.begstageCS
    M_var.DecTree_Relslevel[0] = np.nan
    M_var.DecTree_Relslevel[1] = np.nan
    if (
        startdate.month == LO_Model["date"].iloc[2].month
        and startdate.day == LO_Model["date"].iloc[2].day
    ):
        X1 = "SimDay1"
    elif (
        begdateCS.year == LO_Model["date"].iloc[2].year
        and begdateCS.month == LO_Model["date"].iloc[2].month
        and begdateCS.day == LO_Model["date"].iloc[2].day
    ):
        X1 = "CS start date"
    else:
        X1 = LO_Model["date"].iloc[2]
    M_var.DayFlags[2] = X1
    StartStorage = utils.stg_sto_ar.stg2sto(Pre_defined_Variables.startstage, 0)
    M_var.Storage[0] = StartStorage
    M_var.Storage[1] = StartStorage
    # Flood = np.zeros(n_rows, dtype = object)
    ##Here, I will insert the Storage Deviaiton Values as Input!
    Storage_dev = Data.Stroage_dev_df["DS_dev"]
    # Create a Choose Function for AP Post Baseflow
    # if Pre_defined_Variables.Opt_AdapProt == 0:
    #     C = 450
    # elif Pre_defined_Variables.Opt_AdapProt == 1:
    #     C = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    # Choose_1 = C
    Choose_1 = 450  # cfs
    ####################################################################################################################
    Load_ext = pd.read_csv(
        "./Data/LO_External_Loadings_3MLag_%s.csv" % Pre_defined_Variables.Schedule
    )
    Q_in = pd.read_csv("./Data/LO_Inflows_BK_%s.csv" % Pre_defined_Variables.Schedule)

    ##############################################################################################################
    L_ext = Load_ext["TP_Loads_In_mg"]  # mg
    Atm_Dep_N = TP_Variables.N_Per * Load_ext["Atm_Loading_mg"]
    Atm_Dep_S = TP_Variables.S_Per * Load_ext["Atm_Loading_mg"]
    # Q_Out = pd.read_csv('./Data/Outflows_consd_20082018.csv')
    # C_rain = 10.417 #TP Rainfall Concentration (µg P L-1 = mg P /m3)
    # L_drdep = 0.0385 # mg P / m2 / day
    # Atm_Dep_N = TP_Variables.N_Per * (C_rain*RF_Vol*1233.48 + L_drdep*LO_Area*4046.85642)
    # Atm_Dep_S = TP_Variables.S_Per * (C_rain*RF_Vol*1233.48 + L_drdep*LO_Area*4046.85642)
    # Atm_Dep_N = TP_Variables.N_Per*(18/365)*LO_Area*4046.85642 #Based on data presented by Curtis Pollman, the Lake Okeechobee Technical Advisory Committee (2000) recommended that 18 mgP/m2-yr is an appropriate atmospheric loading of phosphorus over the open lake.
    # Atm_Dep_S = TP_Variables.S_Per*(18/365)*LO_Area*4046.85642
    # Read Shear Stress driven by Wind Speed
    Wind_ShearStr = pd.read_csv(
        "./Data/WindShearStress_%s.csv" % Pre_defined_Variables.Schedule
    )
    W_SS = Wind_ShearStr["ShearStress"]  # Dyne/cm2
    nu_ts = pd.read_csv("./Data/nu_%s.csv" % Pre_defined_Variables.Schedule)
    LO_BL = 0.5  # m (Bed Elevation of LO)
    # LO_WD = pd.to_numeric(Stage_Storage['Stage_m'])-LO_BL
    g = 9.8  # m/s2 gravitational acceleration
    Cal_Res = pd.read_csv(
        "C:/Osama_PC/LOONE/Model/LOONE_Model/Data/nondominated_Sol_var.csv"
    )
    Par = Cal_Res["Par"]
    d_c = Par[20]  # m (particle diameter 10 microm /1E6 to convert to m) clay
    d_s = Par[21]  # m sand
    nu_d = nu_ts["nu"]
    # LO_Temp = 1.0034/1E6 # m2/s (kinematic viscosity of water at T = 20 C)
    # water_density = 1 # g/cm3
    # a = 20.0
    # n = 0.9
    # b = 2.5
    # m = 1.2
    R = 1.65  # submerged specific gravity (1.65 for quartz in water)
    C_1_c = Par[16]
    C_2_c = Par[17]
    C_1_s = Par[18]
    C_2_s = Par[19]
    # Parameters associated with sediment resuspension
    E_0 = 1e-4
    E_1 = 2
    E_2 = 3
    Crtcl_ShStr = Par[22]  # 0.32 #Dyne/cm2
    Td = Par[23]  # days
    n_rows = len(Load_ext.index)
    L_ext_M = np.zeros(n_rows, dtype=object)
    Q_N2S = np.zeros(n_rows, dtype=object)
    # Stage_LO = Stage_Storage['Stage_ft']
    # Storage = Stage_Storage['Storage_acft']
    LO_WD = np.zeros(n_rows, dtype=object)
    Lake_O_Storage_N = np.zeros(n_rows, dtype=object)
    Lake_O_Storage_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_M_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_S_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_R_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_P_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_M_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_S_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_R_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_P_S = np.zeros(n_rows, dtype=object)
    DIP_Lake_N = np.zeros(n_rows, dtype=object)
    DIP_Lake_S = np.zeros(n_rows, dtype=object)
    TP_Lake_Mean = np.zeros(n_rows, dtype=object)
    J_des_M_N = np.zeros(n_rows, dtype=object)
    J_des_S_N = np.zeros(n_rows, dtype=object)
    J_des_R_N = np.zeros(n_rows, dtype=object)
    J_des_P_N = np.zeros(n_rows, dtype=object)
    J_des_M_S = np.zeros(n_rows, dtype=object)
    J_des_S_S = np.zeros(n_rows, dtype=object)
    J_des_R_S = np.zeros(n_rows, dtype=object)
    J_des_P_S = np.zeros(n_rows, dtype=object)
    J_ads_M_N = np.zeros(n_rows, dtype=object)
    J_ads_S_N = np.zeros(n_rows, dtype=object)
    J_ads_R_N = np.zeros(n_rows, dtype=object)
    J_ads_P_N = np.zeros(n_rows, dtype=object)
    J_ads_M_S = np.zeros(n_rows, dtype=object)
    J_ads_S_S = np.zeros(n_rows, dtype=object)
    J_ads_R_S = np.zeros(n_rows, dtype=object)
    J_ads_P_S = np.zeros(n_rows, dtype=object)
    P_sed_M_N = np.zeros(n_rows, dtype=object)
    P_sed_S_N = np.zeros(n_rows, dtype=object)
    P_sed_R_N = np.zeros(n_rows, dtype=object)
    P_sed_P_N = np.zeros(n_rows, dtype=object)
    P_sed_M_S = np.zeros(n_rows, dtype=object)
    P_sed_S_S = np.zeros(n_rows, dtype=object)
    P_sed_R_S = np.zeros(n_rows, dtype=object)
    P_sed_P_S = np.zeros(n_rows, dtype=object)
    J_sedburial_M_N = np.zeros(n_rows, dtype=object)
    J_sedburial_S_N = np.zeros(n_rows, dtype=object)
    J_sedburial_R_N = np.zeros(n_rows, dtype=object)
    J_sedburial_P_N = np.zeros(n_rows, dtype=object)
    J_sedburial_M_S = np.zeros(n_rows, dtype=object)
    J_sedburial_S_S = np.zeros(n_rows, dtype=object)
    J_sedburial_R_S = np.zeros(n_rows, dtype=object)
    J_sedburial_P_S = np.zeros(n_rows, dtype=object)
    J_Γburial_M_N = np.zeros(n_rows, dtype=object)
    J_Γburial_S_N = np.zeros(n_rows, dtype=object)
    J_Γburial_R_N = np.zeros(n_rows, dtype=object)
    J_Γburial_P_N = np.zeros(n_rows, dtype=object)
    J_Γburial_M_S = np.zeros(n_rows, dtype=object)
    J_Γburial_S_S = np.zeros(n_rows, dtype=object)
    J_Γburial_R_S = np.zeros(n_rows, dtype=object)
    J_Γburial_P_S = np.zeros(n_rows, dtype=object)
    Γ_M_N = np.zeros(n_rows, dtype=object)
    Γ_S_N = np.zeros(n_rows, dtype=object)
    Γ_R_N = np.zeros(n_rows, dtype=object)
    Γ_P_N = np.zeros(n_rows, dtype=object)
    Γ_M_S = np.zeros(n_rows, dtype=object)
    Γ_S_S = np.zeros(n_rows, dtype=object)
    Γ_R_S = np.zeros(n_rows, dtype=object)
    Γ_P_S = np.zeros(n_rows, dtype=object)
    DIP_pore_M_N = np.zeros(n_rows, dtype=object)
    DIP_pore_S_N = np.zeros(n_rows, dtype=object)
    DIP_pore_R_N = np.zeros(n_rows, dtype=object)
    DIP_pore_P_N = np.zeros(n_rows, dtype=object)
    DIP_pore_M_S = np.zeros(n_rows, dtype=object)
    DIP_pore_S_S = np.zeros(n_rows, dtype=object)
    DIP_pore_R_S = np.zeros(n_rows, dtype=object)
    DIP_pore_P_S = np.zeros(n_rows, dtype=object)
    TP_Lake_N = np.zeros(n_rows, dtype=object)
    TP_Lake_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_M_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_S_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_R_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_P_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_M_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_S_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_R_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_P_S = np.zeros(n_rows, dtype=object)

    J_decomp_M_N = np.zeros(n_rows, dtype=object)
    J_decomp_S_N = np.zeros(n_rows, dtype=object)
    J_decomp_R_N = np.zeros(n_rows, dtype=object)
    J_decomp_P_N = np.zeros(n_rows, dtype=object)
    J_decomp_M_S = np.zeros(n_rows, dtype=object)
    J_decomp_S_S = np.zeros(n_rows, dtype=object)
    J_decomp_R_S = np.zeros(n_rows, dtype=object)
    J_decomp_P_S = np.zeros(n_rows, dtype=object)

    Settling_P_N = np.zeros(n_rows, dtype=object)
    Settling_P_S = np.zeros(n_rows, dtype=object)

    P_diff_M_N = np.zeros(n_rows, dtype=object)
    P_diff_S_N = np.zeros(n_rows, dtype=object)
    P_diff_R_N = np.zeros(n_rows, dtype=object)
    P_diff_P_N = np.zeros(n_rows, dtype=object)
    P_diff_M_S = np.zeros(n_rows, dtype=object)
    P_diff_S_S = np.zeros(n_rows, dtype=object)
    P_diff_R_S = np.zeros(n_rows, dtype=object)
    P_diff_P_S = np.zeros(n_rows, dtype=object)

    # TP_N_to_S = np.zeros(n_rows,dtype = object)
    # TP_Out = np.zeros(n_rows,dtype = object)
    # L_Ext_mgperm3 = np.zeros(n_rows,dtype = object)
    Q_I = Q_in["Flow_cmd"]
    Q_I_M = np.zeros(n_rows, dtype=object)
    Q_O = np.zeros(n_rows, dtype=object)
    # Indust_O = pd.read_csv('./Data/INDUST_Outflow_20082018.csv')
    # Q_O = Q_Out['Total_Outflows_acft'] * 1233.48 + Indust_O['INDUST_cmd']
    Q_O_M = np.zeros(n_rows, dtype=object)

    P_Load_Cal = np.zeros(n_rows, dtype=object)
    P_Load_StL = np.zeros(n_rows, dtype=object)
    P_Load_South = np.zeros(n_rows, dtype=object)

    # Ferguson, R. I., and Church, M. (2004).
    # v_settle_N_c = (R*g*d_c**2)/(C_1_c*nu+(0.75*C_2_c*R*g*d_c**3)**0.5)
    # v_settle_N_s = (R*g*d_s**2)/(C_1_s*nu+(0.75*C_2_s*R*g*d_s**3)**0.5)
    # v_settle_N = v_settle_N_c*((TP_Variables.A_Mud_N+TP_Variables.A_Peat_N)/TP_Variables.A_N) + v_settle_N_s*((TP_Variables.A_Sand_N + TP_Variables.A_Rock_N)/TP_Variables.A_N)

    # v_settle_S_c = (R*g*d_c**2)/(C_1_c*nu+(0.75*C_2_c*R*g*d_c**3)**0.5)
    # v_settle_S_s = (R*g*d_s**2)/(C_1_s*nu+(0.75*C_2_s*R*g*d_s**3)**0.5)
    # v_settle_S = v_settle_S_c*((TP_Variables.A_Mud_S+TP_Variables.A_Peat_S)/TP_Variables.A_S) + v_settle_S_s*((TP_Variables.A_Sand_S + TP_Variables.A_Rock_S)/TP_Variables.A_S)

    v_settle_N_c = np.zeros(n_rows, dtype=object)
    v_settle_N_s = np.zeros(n_rows, dtype=object)
    v_settle_N = np.zeros(n_rows, dtype=object)
    v_settle_S_c = np.zeros(n_rows, dtype=object)
    v_settle_S_s = np.zeros(n_rows, dtype=object)
    v_settle_S = np.zeros(n_rows, dtype=object)

    # v_settle_N = np.zeros(n_rows,dtype = object)
    # v_settle_S = np.zeros(n_rows,dtype = object)
    #####################################################################################################
    ##Initial Values##
    # S.A. is calculated based on the Lake's previous time step Stage, but for the S.A. at i=0 I used same time step Stage!

    Q_O[0] = 0
    Q_O[1] = 46485  # cmd
    # TP_MassBalanceModel Initial Values.
    TP_Lake_N[0] = 225  # mg/m3
    TP_Lake_S[0] = 275  # mg/m3
    TP_Lake_Mean[0] = (TP_Lake_N[0] + TP_Lake_S[0]) / 2
    Γ_M_N[0] = 25  # mg/kg
    Γ_S_N[0] = 25  # mg/kg
    Γ_R_N[0] = 25  # mg/kg
    Γ_P_N[0] = 25  # mg/kg
    Γ_M_S[0] = 25  # mg/kg
    Γ_S_S[0] = 25  # mg/kg
    Γ_R_S[0] = 25  # mg/kg
    Γ_P_S[0] = 25  # mg/kg
    DIP_pore_M_N[0] = 700  # 760 #mg/m3
    DIP_pore_S_N[0] = 240  # 205 #mg/m3
    DIP_pore_R_N[0] = 240  # 205 #mg/m3
    DIP_pore_P_N[0] = 160  # 160 #mg/m3
    DIP_pore_M_S[0] = 700  # 760 #mg/m3
    DIP_pore_S_S[0] = 240  # 205 #mg/m3
    DIP_pore_R_S[0] = 240  # 205 #mg/m3
    DIP_pore_P_S[0] = 160  # 160 #mg/m3
    P_sed_M_N[0] = 1100  # mg/kg
    P_sed_S_N[0] = 300  # mg/kg
    P_sed_R_N[0] = 300  # mg/kg
    P_sed_P_N[0] = 200  # mg/kg
    P_sed_M_S[0] = 1100  # mg/kg
    P_sed_S_S[0] = 300  # mg/kg
    P_sed_R_S[0] = 300  # mg/kg
    P_sed_P_S[0] = 200  # mg/kg
    Θ_M = 1 - (
        (TP_Variables.Bulk_density_M / TP_Variables.Particle_density_M)
        * ((100 - TP_Variables.Per_H2O_M) / 100)
    )
    Θ_S = 1 - (
        (TP_Variables.Bulk_density_S / TP_Variables.Particle_density_S)
        * ((100 - TP_Variables.Per_H2O_S) / 100)
    )
    Θ_R = 1 - (
        (TP_Variables.Bulk_density_R / TP_Variables.Particle_density_R)
        * ((100 - TP_Variables.Per_H2O_R) / 100)
    )
    Θ_P = 1 - (
        (TP_Variables.Bulk_density_P / TP_Variables.Particle_density_P)
        * ((100 - TP_Variables.Per_H2O_P) / 100)
    )
    # Mass of sediment in surfacial mix Mud layer in the North Region(kg)
    Mass_sed_M_N = (
        TP_Variables.A_Mud_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_M) / 100)
        * TP_Variables.Bulk_density_M
        * 1000
    )
    # Mass of sediment in surfacial mix Sand layer in the North Region(kg)
    Mass_sed_S_N = (
        TP_Variables.A_Sand_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_S) / 100)
        * TP_Variables.Bulk_density_S
        * 1000
    )
    # Mass of sediment in surfacial mix Rock layer in the North Region(kg)
    Mass_sed_R_N = (
        TP_Variables.A_Rock_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_R) / 100)
        * TP_Variables.Bulk_density_R
        * 1000
    )
    # Mass of sediment in surfacial mix Peat layer in the North Region(kg)
    Mass_sed_P_N = (
        TP_Variables.A_Peat_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_P) / 100)
        * TP_Variables.Bulk_density_P
        * 1000
    )
    # Mass of sediment in surfacial mix Mud layer in the South Region(kg)
    Mass_sed_M_S = (
        TP_Variables.A_Mud_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_M) / 100)
        * TP_Variables.Bulk_density_M
        * 1000
    )
    # Mass of sediment in surfacial mix Sand layer in the South Region(kg)
    Mass_sed_S_S = (
        TP_Variables.A_Sand_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_S) / 100)
        * TP_Variables.Bulk_density_S
        * 1000
    )
    # Mass of sediment in surfacial mix Rock layer in the South Region(kg)
    Mass_sed_R_S = (
        TP_Variables.A_Rock_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_R) / 100)
        * TP_Variables.Bulk_density_R
        * 1000
    )
    # Mass of sediment in surfacial mix Peat layer in the South Region(kg)
    Mass_sed_P_S = (
        TP_Variables.A_Peat_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_P) / 100)
        * TP_Variables.Bulk_density_P
        * 1000
    )

    ######################################################################################################################################
    M_var.Zone_Code[0] = utils.lo_functions.Zone_Code(
        M_var.Lake_Stage[0],
        df_WSMs["A"].iloc[0],
        df_WSMs["B"].iloc[0],
        df_WSMs["C"].iloc[0],
        df_WSMs["D3"].iloc[0],
        df_WSMs["D2"].iloc[0],
        df_WSMs["D1"].iloc[0],
        df_WSMs["D0"].iloc[0],
        df_WSMs["WSM1"].iloc[0],
    )
    M_var.LO_Zone[0] = utils.lo_functions.LO_Zone(M_var.Zone_Code[0])
    for i in range(n_rows - 2):
        M_var.WSM_Zone[i + 2] = utils.lo_functions.WSM_Zone(
            M_var.Lake_Stage[i + 1],
            df_WSMs.at[i + 1, "WSM4"],
            df_WSMs.at[i + 1, "WSM3"],
            df_WSMs.at[i + 1, "WSM2"],
            df_WSMs.at[i + 1, "WSM1"],
        )
        # Calculate Daily Maximum Water Supply
        # Note that in LOSA_dmd we used (i) because this file starts from 1/1/2008 so i at this point =0.
        # Cutbacks are determined based on the WSM Zone.
        M_var.Max_Supply[i + 2] = utils.lo_functions.Max_Supply(
            M_var.WSM_Zone[i + 2],
            Water_dmd.at[i, "Daily_demand"],
            Pre_defined_Variables.Z1_cutback,
            Pre_defined_Variables.Z2_cutback,
            Pre_defined_Variables.Z3_cutback,
            Pre_defined_Variables.Z4_cutback,
        )
        # Actual Daily Water supply
        M_var.LOSA_Supply[i + 2] = utils.lo_functions.LOSA_Supply(
            M_var.WSM_Zone[i + 2],
            LO_Model.at[i + 2, "LOSA_dmd_SFWMM"],
            M_var.Max_Supply[i + 2],
            Pre_defined_Variables.Opt_LOSAws,
        )
        # NetInflow - LOSA Supply
        M_var.NI_Supply[i + 2] = (
            LO_Model.at[i + 2, "Net_Inflow"] - M_var.LOSA_Supply[i + 2]
        )
        # TODO Note: for the pass statement, We will read the Daily Water supply from the SFWMM as an input.
        # Calculate the cutback where Cutback = Demand - Supply
        ctbk = LO_Model.at[i + 2, "LOSA_dmd_SFWMM"] - M_var.LOSA_Supply[i + 2]
        M_var.Cut_back[i + 2] = ctbk
        # Calculate percentage of the demand that is not supplied for each day
        if LO_Model.at[i + 2, "LOSA_dmd_SFWMM"] == 0:
            DNS = 0
        else:
            DNS = (M_var.Cut_back[i + 2] / LO_Model.at[i + 2, "LOSA_dmd_SFWMM"]) * 100
        M_var.Dem_N_Sup[i + 2] = DNS
        # Calculate the Zone Code
        # Note that to calculate the Zone Code in Dec 31 2020 we needed the WSM and breakpoint zones in 1/1/2021!
        # Note Also that i = 0 in Stage indicates Dec 30 1964 while i = 0 in df_WSMs indicates Dec 31 1964!
        M_var.Zone_Code[i + 1] = utils.lo_functions.Zone_Code(
            M_var.Lake_Stage[i + 1],
            df_WSMs.at[i + 1, "A"],
            df_WSMs.at[i + 1, "B"],
            df_WSMs.at[i + 1, "C"],
            df_WSMs.at[i + 1, "D3"],
            df_WSMs.at[i + 1, "D2"],
            df_WSMs.at[i + 1, "D1"],
            df_WSMs.at[i + 1, "D0"],
            df_WSMs.at[i + 1, "WSM1"],
        )
        # Generate the Zone Column based on the corresponding Zone Code.
        M_var.LO_Zone[i + 1] = utils.lo_functions.LO_Zone(M_var.Zone_Code[i + 1])
        M_var.Zone_D_Trib[i] = utils.dec_tree_functions.Zone_D_Trib(
            TC_LONINO_df.at[i, "Tributary_Condition"], Pre_defined_Variables.Opt_NewTree
        )
        M_var.Zone_D_stage[i] = utils.dec_tree_functions.Zone_D_stage(
            M_var.Lake_Stage[i + 1], df_WSMs.at[i, "C-b"]
        )
        M_var.Zone_D_Seas[i] = utils.dec_tree_functions.Zone_D_Seas(
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            M_var.Zone_D_Trib[i],
            Pre_defined_Variables.Opt_NewTree,
        )
        M_var.Zone_D_MSeas[i] = utils.dec_tree_functions.Zone_D_MSeas(
            TC_LONINO_df.at[i, "LONINO_MultiSeasonal_Classes"]
        )
        M_var.Zone_D_Branch_Code[i] = (
            M_var.Zone_D_Trib[i] * 1000
            + M_var.Zone_D_stage[i] * 100
            + M_var.Zone_D_Seas[i] * 10
            + M_var.Zone_D_MSeas[i] * 1
        )
        M_var.Zone_D_Rel_Code[i] = utils.dec_tree_functions.Zone_D_Rel_Code(
            M_var.Zone_D_Branch_Code[i], Pre_defined_Variables.Opt_DecTree
        )
        M_var.Zone_C_Trib[i] = utils.dec_tree_functions.Zone_C_Trib(
            TC_LONINO_df.at[i, "Tributary_Condition"], Pre_defined_Variables.Opt_NewTree
        )
        M_var.Zone_C_Seas[i] = utils.dec_tree_functions.Zone_C_Seas(
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            Pre_defined_Variables.Opt_NewTree,
        )
        M_var.Zone_C_MSeas[i] = utils.dec_tree_functions.Zone_C_MSeas(
            TC_LONINO_df.at[i, "LONINO_MultiSeasonal_Classes"]
        )
        M_var.Zone_C_MetFcast[i] = utils.dec_tree_functions.Zone_C_MetFcast(
            M_var.Zone_C_Seas[i],
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            Pre_defined_Variables.Zone_C_MetFcast_Indicator,
        )
        M_var.Zone_C_Branch_Code[i] = (
            M_var.Zone_C_Trib[i] * 1000
            + M_var.Zone_C_MetFcast[i] * 100
            + M_var.Zone_C_Seas[i] * 10
            + M_var.Zone_C_MSeas[i] * 1
        )
        M_var.Zone_C_Rel_Code[i] = utils.dec_tree_functions.Zone_C_Rel_Code(
            M_var.Zone_C_Branch_Code[i], Pre_defined_Variables.Opt_DecTree
        )
        M_var.Zone_B_Trib[i] = utils.dec_tree_functions.Zone_B_Trib(
            TC_LONINO_df.at[i, "Tributary_Condition"], Pre_defined_Variables.Opt_NewTree
        )
        M_var.Zone_B_Stage[i] = utils.dec_tree_functions.Zone_B_Stage(
            M_var.Lake_Stage[i + 1], Seasons.at[i, "Season"]
        )
        M_var.Zone_B_Seas[i] = utils.dec_tree_functions.Zone_B_Seas(
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"]
        )
        M_var.Zone_B_Branch_Code[i] = (
            M_var.Zone_B_Trib[i] * 1000
            + M_var.Zone_B_Stage[i] * 100
            + DecTree_df.at[i, "Zone_B_MetFcast"] * 10
            + M_var.Zone_B_Seas[i] * 1
        )
        M_var.Zone_B_Rel_Code[i] = utils.dec_tree_functions.Zone_B_Rel_Code(
            M_var.Zone_B_Branch_Code[i], Pre_defined_Variables.Opt_DecTree
        )
        M_var.DecTree_Relslevel[i + 2] = utils.lo_functions.DecTree_Relslevel(
            M_var.Zone_Code[i + 1],
            M_var.Zone_D_Rel_Code[i],
            M_var.Zone_C_Rel_Code[i],
            M_var.Zone_B_Rel_Code[i],
        )
        if i >= 3:
            if (
                startdate.month == LO_Model.at[i, "date"].month
                and startdate.day == LO_Model.at[i, "date"].day
                and (
                    Pre_defined_Variables.CSflag == 0
                    or startdate.year == LO_Model.at[i, "date"].year
                )
            ):
                X2 = "SimDay1"
            else:
                X2 = LO_Model.at[i, "date"].date()
            M_var.DayFlags[i] = X2
        M_var.PlsDay[i + 2] = utils.lo_functions.PlsDay(
            M_var.DayFlags[i + 2],
            M_var.DecTree_Relslevel[i + 2],
            Pre_defined_Variables.PlsDay_Switch,
        )
        M_var.Release_Level[i + 2] = utils.lo_functions.Release_Level(
            M_var.Release_Level[i + 1],
            M_var.Lake_Stage[i + 1],
            TC_LONINO_df.at[i, "Tributary_Condition"],
            M_var.PlsDay[i + 2],
            M_var.Zone_Code[i + 1],
            M_var.DecTree_Relslevel[i + 2],
            Pre_defined_Variables.MaxQstgTrigger,
        )
        if i >= 6:
            dh = M_var.Lake_Stage[i + 1] - M_var.Lake_Stage[i - 6]
            M_var.dh_7days[i + 1] = dh
        M_var.ZoneCodeminus1Code[i + 1] = utils.lo_functions.ZoneCodeminus1Code(
            M_var.Zone_Code[i + 1],
            df_WSMs.at[i + 1, "WSM1"],
            df_WSMs.at[i + 1, "D0"],
            df_WSMs.at[i + 1, "D1"],
            df_WSMs.at[i + 1, "D2"],
            df_WSMs.at[i + 1, "D3"],
            df_WSMs.at[i + 1, "C"],
            df_WSMs.at[i + 1, "B"],
            df_WSMs.at[i + 1, "A"],
        )
        M_var.ZoneCodeCode[i + 1] = utils.lo_functions.ZoneCodeCode(
            M_var.Zone_Code[i + 1],
            df_WSMs.at[i + 1, "WSM1"],
            df_WSMs.at[i + 1, "D0"],
            df_WSMs.at[i + 1, "D1"],
            df_WSMs.at[i + 1, "D2"],
            df_WSMs.at[i + 1, "D3"],
            df_WSMs.at[i + 1, "C"],
            df_WSMs.at[i + 1, "B"],
            df_WSMs.at[i + 1, "A"],
        )
        M_var.Fraction_of_Zone_height[i + 1] = (
            utils.lo_functions.Fraction_of_Zone_height(
                M_var.Zone_Code[i + 1],
                M_var.Lake_Stage[i + 1],
                M_var.ZoneCodeminus1Code[i + 1],
                M_var.ZoneCodeCode[i + 1],
            )
        )
        M_var.ReLevelCode_1[i + 2] = utils.lo_functions.ReLevelCode_1(
            M_var.Release_Level[i + 2],
            Pre_defined_Variables.dstar_D1,
            Pre_defined_Variables.dstar_D2,
            Pre_defined_Variables.dstar_D3,
            Pre_defined_Variables.dstar_C,
            Pre_defined_Variables.dstar_B,
        )
        M_var.ReLevelCode_2[i + 2] = utils.lo_functions.ReLevelCode_2(
            M_var.Release_Level[i + 2],
            Pre_defined_Variables.astar_D1,
            Pre_defined_Variables.astar_D2,
            Pre_defined_Variables.astar_D3,
            Pre_defined_Variables.astar_C,
            Pre_defined_Variables.astar_B,
        )
        M_var.ReLevelCode_3_S80[i + 2] = utils.lo_functions.ReLevelCode_3_S80(
            M_var.Release_Level[i + 2],
            Pre_defined_Variables.bstar_S80_D1,
            Pre_defined_Variables.bstar_S80_D2,
            Pre_defined_Variables.bstar_S80_D3,
            Pre_defined_Variables.bstar_S80_C,
            Pre_defined_Variables.bstar_S80_B,
        )
        M_var.Outlet2DS_Mult[i + 2] = utils.lo_functions.Outlet2DS_Mult(
            Seasons.at[i, "Season"],
            Seasons.at[i, "Month"],
            M_var.dh_7days[i + 1],
            M_var.ReLevelCode_1[i + 2],
            M_var.Fraction_of_Zone_height[i + 1],
            M_var.ReLevelCode_2[i + 2],
            M_var.ReLevelCode_3_S80[i + 2],
            Pre_defined_Variables.Opt_QregMult,
        )
        M_var.Outlet2DS_Mult_2[i + 2] = utils.lo_functions.Outlet2DS_Mult_2(
            LO_Model.at[i + 2, "date"].month,
            LO_Model.at[i + 2, "date"].day,
            M_var.PlsDay[i + 2],
            M_var.Outlet2DS_Mult[i + 2 - M_var.PlsDay[i + 2]],
            M_var.Outlet2DS_Mult[i + 2],
            Pre_defined_Variables.Opt_QregMult,
        )
        M_var.Outlet2DSRS[i + 2] = utils.lo_functions.Outlet2DSRS(
            M_var.Release_Level[i + 2],
            Data.S80_RegRelRates.at[0, "Zone_D1"],
            S80avgL1,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                "S-80_L1_%s" % Pre_defined_Variables.Schedule,
            ],
            M_var.Outlet2DS_Mult_2[i + 2],
            Data.CE_SLE_turns.at[
                LO_Model.at[i + 2, "date"].year - Pre_defined_Variables.startyear,
                "SLEturn",
            ],
            Data.S80_RegRelRates.at[0, "Zone_D2"],
            S80avgL2,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                "S-80_L2_%s" % Pre_defined_Variables.Schedule,
            ],
            Data.S80_RegRelRates.at[0, "Zone_D3"],
            S80avgL3,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                "S-80_L3_%s" % Pre_defined_Variables.Schedule,
            ],
            Data.S80_RegRelRates.at[0, "Zone_C"],
            Data.S80_RegRelRates.at[0, "Zone_B"],
            Data.S80_RegRelRates.at[0, "Zone_A"],
        )
        M_var.Outlet2USRG1[i + 2] = max(
            0, M_var.Outlet2DSRS[i + 2] - LO_Model.at[i + 2, "C44RO"]
        )
        M_var.Sum_Outlet2USRG1[i + 2] = utils.lo_functions.Sum_Outlet2USRG1(
            LO_Model.at[i + 2, "date"].day, M_var.Outlet2USRG1[i + 2]
        )
        M_var.Outlet2DSBS[i + 2] = utils.lo_functions.Outlet2DSBS(
            M_var.Release_Level[i + 2],
            M_var.Sum_Outlet2USRG1[i + 2],
            VLOOKUP1_c[i],
            Outlet2_baseflow,
            Pre_defined_Variables.Option_S80Baseflow,
        )
        M_var.Outlet2USBK[i + 2] = utils.lo_functions.Outlet2USBK(
            M_var.Lake_Stage[i + 1],
            df_WSMs.at[i + 1, "D1"],
            M_var.Outlet2USRG[i + 1],
            LO_Model.at[i + 2, "C44RO"],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S308BK"],
            Pre_defined_Variables.Opt_S308,
            Pre_defined_Variables.S308BK_Const,
            Pre_defined_Variables.S308_BK_Thr,
        )
        M_var.ROeast[i + 2] = LO_Model.at[i + 2, "C44RO"] - M_var.Outlet2USBK[i + 2]
        M_var.Outlet2USBS[i + 2] = utils.lo_functions.Outlet2USBS(
            M_var.Outlet2DSBS[i + 2],
            M_var.Outlet2USRG1[i + 2],
            M_var.ROeast[i + 2],
            Pre_defined_Variables.Option_S80Baseflow,
        )
        M_var.Sum_Outlet2USBK[i + 2] = utils.lo_functions.Sum_Outlet2USBK(
            LO_Model.at[i + 2, "date"].day, M_var.Outlet2USBK[i + 2]
        )
        M_var.Outlet2USRG_Code[i + 2] = utils.lo_functions.Outlet2USRG_Code(
            M_var.Outlet2USRG1[i + 2],
            M_var.Outlet2USBS[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S308RG"],
            Data.SFWMM_Daily_Outputs.at[i + 2, "STEST"],
            Pre_defined_Variables.Option_RegS77S308,
        )
        if Model_Config.Sim_type == 0:
            M_var.Outlet2USRG[i + 2] = utils.lo_functions.Outlet2USRG(
                M_var.Outlet2USRG_Code[i + 2],
                Data.SFWMM_Daily_Outputs.at[i + 2, "S308RG"],
                Data.SFWMM_Daily_Outputs.at[i + 2, "STEST"],
                Pre_defined_Variables.Opt_S308,
                Pre_defined_Variables.S308RG_Const,
            )
        else:
            if M_var.Lake_Stage[i + 1] >= 18:
                M_var.Outlet2USRG[i + 2] = 7200
            elif M_var.Lake_Stage[i + 1] <= 8:
                M_var.Outlet2USRG[i + 2] = 0
            elif (TP_Lake_S[i] <= P_1) and (
                date_rng_6[i + 2].month in [1, 2, 3, 4, 11, 12]
            ):
                M_var.Outlet2USRG[i + 2] = S308_DV[(date_rng_6[i + 2].month) - 1]
            elif (TP_Lake_S[i] <= P_2) and (
                date_rng_6[i + 2].month in [5, 6, 7, 8, 9, 10]
            ):
                M_var.Outlet2USRG[i + 2] = S308_DV[(date_rng_6[i + 2].month) - 1]
            else:
                M_var.Outlet2USRG[i + 2] = 0
        M_var.Outlet2DS[i + 2] = utils.lo_functions.S80(
            M_var.ROeast[i + 2],
            M_var.Outlet2USRG[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S80"],
            Pre_defined_Variables.S80_Const,
        )
        M_var.ReLevelCode_3_S77[i + 2] = utils.lo_functions.ReLevelCode_3_S77(
            M_var.Release_Level[i + 2],
            Pre_defined_Variables.bstar_S77_D1,
            Pre_defined_Variables.bstar_S77_D2,
            Pre_defined_Variables.bstar_S77_D3,
            Pre_defined_Variables.bstar_S77_C,
            Pre_defined_Variables.bstar_S77_B,
        )
        M_var.Outlet1US_Mult[i + 2] = utils.lo_functions.Outlet1US_Mult(
            Seasons.at[i, "Season"],
            Seasons.at[i, "Month"],
            M_var.dh_7days[i + 1],
            M_var.ReLevelCode_1[i + 2],
            M_var.Fraction_of_Zone_height[i + 1],
            M_var.ReLevelCode_2[i + 2],
            M_var.ReLevelCode_3_S77[i + 2],
            Pre_defined_Variables.Opt_QregMult,
        )
        M_var.Outlet1US_Mult_2[i + 2] = utils.lo_functions.Outlet1US_Mult_2(
            LO_Model.at[i + 2, "date"].month,
            LO_Model.at[i + 2, "date"].day,
            M_var.PlsDay[i + 2],
            M_var.Outlet1US_Mult[i + 2 - M_var.PlsDay[i + 2]],
            M_var.Outlet1US_Mult[i + 2],
            Pre_defined_Variables.Opt_QregMult,
        )
        M_var.Outlet1USRS[i + 2] = utils.lo_functions.Outlet1USRS(
            M_var.Release_Level[i + 2],
            Data.S77_RegRelRates.at[0, "Zone_D1"],
            S77avgL1,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                "S-77_L1_%s" % Pre_defined_Variables.Schedule,
            ],
            M_var.Outlet1US_Mult_2[i + 2],
            LO_Model.at[i + 2, "C43RO"],
            Data.CE_SLE_turns.at[
                LO_Model.at[i + 2, "date"].year - Pre_defined_Variables.startyear,
                "CEturn",
            ],
            Data.S77_RegRelRates.at[0, "Zone_D2"],
            S77avgL2,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                "S-77_L2_%s" % Pre_defined_Variables.Schedule,
            ],
            M_var.Zone_Code[i + 1],
            Data.S77_RegRelRates.at[0, "Zone_D3"],
            S77avgL3,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                "S-77_L3_%s" % Pre_defined_Variables.Schedule,
            ],
            Data.S77_RegRelRates.at[0, "Zone_C"],
            Data.S77_RegRelRates.at[0, "Zone_B"],
            Data.S77_RegRelRates.at[0, "Zone_A"],
            Pre_defined_Variables.Opt_Outlet1DSRG,
        )
        M_var.Sum_Outlet1USRS[i + 2] = utils.lo_functions.Sum_Outlet1USRS(
            LO_Model.at[i + 2, "date"].day, M_var.Outlet1USRS[i + 2]
        )
        M_var.Outlet1USBK[i + 2] = utils.lo_functions.Outlet1USBK(
            M_var.Lake_Stage[i + 1],
            M_var.Outlet1USRS[i + 2],
            M_var.Outlet1USBSAP[i + 1],
            M_var.Outlet1USEWS[i + 1],
            LO_Model.at[i + 2, "C43RO"],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S77BK"],
            Pre_defined_Variables.Outlet1USBK_Switch,
            Pre_defined_Variables.Outlet1USBK_Threshold,
        )
        M_var.ROwest[i + 2] = LO_Model.at[i + 2, "C43RO"] - M_var.Outlet1USBK[i + 2]
        M_var.Outlet1DSBS[i + 2] = utils.lo_functions.Outlet1DSBS(
            M_var.Release_Level[i + 2],
            M_var.Sum_Outlet1USRS[i + 2],
            VLOOKUP2_c[i],
            Outlet1_baseflow,
            Pre_defined_Variables.Option_S77Baseflow,
        )
        M_var.Outlet1USBS[i + 2] = utils.lo_functions.Outlet1USBS(
            M_var.Outlet1DSBS[i + 2],
            M_var.Outlet1USRS[i + 2],
            M_var.ROwest[i + 2],
            Pre_defined_Variables.Option_S77Baseflow,
        )
        # Define THC Class Normal or above
        if i < (n_rows - 2):
            M_var.Post_Ap_Baseflow[i] = utils.thc_class.THC_Class(
                i,
                M_var.THC_Class_normal_or_above,
                M_var.Lake_O_Stage_AP,
                M_var.Lake_O_Schedule_Zone,
                M_var.LStgCorres,
                M_var.LowChance_Check,
                M_var.Outlet1USRS_AP,
                M_var.Outlet1USBS_AP,
                M_var.Outlet1USRS_Pre_AP_S77_Baseflow,
                M_var.Forecast_D_Sal,
                M_var.n30d_mavg,
                M_var.n30davgForecast,
                M_var.LORS08_bf_rel,
                M_var.LDS_LC6_1,
                M_var.S_O,
                M_var.All_4,
                M_var.Sabf,
                M_var.Swbf,
                M_var.Swbu,
                M_var.All_4andStage,
                M_var.All_4andStagein,
                M_var.P_AP_BF_Stg,
                M_var.Logic_test_1,
                M_var.Post_Ap_Baseflow,
                M_var.Outlet1USRSplusPreAPS77bsf,
                M_var.AndEstNeedsLakeWater,
                M_var.AndLowChance61stagelessth11,
                M_var.ATHCnora,
                M_var.Choose_PAPEWS_1,
                M_var.Choose_PAPEWS_2,
                M_var.Post_AP_EWS,
                M_var.Post_AP_Baseflow_EWS_cfs,
                AdapProt_df,
                M_var.Lake_Stage,
                M_var.Zone_Code,
                df_WSMs,
                Targ_Stg_df,
                M_var.Outlet1USRS,
                M_var.Outlet1USBS,
                Data.Estuary_needs_water,
                Choose_1,
                M_var.WSM_Zone,
            )["Post_Ap_Baseflow"]
            M_var.Post_AP_EWS[i] = utils.thc_class.THC_Class(
                i,
                M_var.THC_Class_normal_or_above,
                M_var.Lake_O_Stage_AP,
                M_var.Lake_O_Schedule_Zone,
                M_var.LStgCorres,
                M_var.LowChance_Check,
                M_var.Outlet1USRS_AP,
                M_var.Outlet1USBS_AP,
                M_var.Outlet1USRS_Pre_AP_S77_Baseflow,
                M_var.Forecast_D_Sal,
                M_var.n30d_mavg,
                M_var.n30davgForecast,
                M_var.LORS08_bf_rel,
                M_var.LDS_LC6_1,
                M_var.S_O,
                M_var.All_4,
                M_var.Sabf,
                M_var.Swbf,
                M_var.Swbu,
                M_var.All_4andStage,
                M_var.All_4andStagein,
                M_var.P_AP_BF_Stg,
                M_var.Logic_test_1,
                M_var.Post_Ap_Baseflow,
                M_var.Outlet1USRSplusPreAPS77bsf,
                M_var.AndEstNeedsLakeWater,
                M_var.AndLowChance61stagelessth11,
                M_var.ATHCnora,
                M_var.Choose_PAPEWS_1,
                M_var.Choose_PAPEWS_2,
                M_var.Post_AP_EWS,
                M_var.Post_AP_Baseflow_EWS_cfs,
                AdapProt_df,
                M_var.Lake_Stage,
                M_var.Zone_Code,
                df_WSMs,
                Targ_Stg_df,
                M_var.Outlet1USRS,
                M_var.Outlet1USBS,
                Data.Estuary_needs_water,
                Choose_1,
                M_var.WSM_Zone,
            )["Post_AP_EWS"]
        M_var.Outlet1USBSAP[i + 2] = utils.lo_functions.Outlet1USBSAP(
            M_var.Outlet1USBS[i + 2],
            M_var.Post_Ap_Baseflow[i],
            Pre_defined_Variables.Opt_AdapProt,
        )
        M_var.Outlet1USEWS[i + 2] = utils.lo_functions.Outlet1USEWS(
            M_var.Post_AP_EWS[i],
            Data.SFWMM_Daily_Outputs.at[i + 2, "CAEST"],
            Pre_defined_Variables.Outlet1USEWS_Switch,
            Pre_defined_Variables.Opt_AdapProt,
        )
        if Model_Config.Sim_type == 0:
            M_var.Outlet1USREG[i + 2] = utils.lo_functions.Outlet1USREG(
                M_var.Outlet1USRS[i + 2],
                M_var.Outlet1USBSAP[i + 2],
                Data.SFWMM_Daily_Outputs.at[i + 2, "S77RG"],
                Pre_defined_Variables.Outlet1USREG_Switch,
                Pre_defined_Variables.Option_RegS77S308,
            )
        else:
            if M_var.Lake_Stage[i + 1] >= 18:
                M_var.Outlet1USREG[i + 2] = 7800
            elif M_var.Lake_Stage[i + 1] <= 8:
                M_var.Outlet1USREG[i + 2] = 0
            elif (TP_Lake_S[i] <= P_1) and (
                date_rng_6[i + 2].month in [1, 2, 3, 4, 11, 12]
            ):
                M_var.Outlet1USREG[i + 2] = S77_DV[(date_rng_6[i + 2].month) - 1]
            elif (TP_Lake_S[i] <= P_2) and (
                date_rng_6[i + 2].month in [5, 6, 7, 8, 9, 10]
            ):
                M_var.Outlet1USREG[i + 2] = S77_DV[(date_rng_6[i + 2].month) - 1]
            else:
                M_var.Outlet1USREG[i + 2] = 0
        M_var.Outlet1DS[i + 2] = utils.lo_functions.Outlet1DS(
            M_var.Outlet1USREG[i + 2],
            M_var.Outlet1USEWS[i + 2],
            M_var.ROwest[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S79"],
            Pre_defined_Variables.Outlet1DS_Switch,
        )
        M_var.TotRegEW[i + 2] = (
            M_var.Outlet1USREG[i + 2] + M_var.Outlet2USRG[i + 2]
        ) * 1.9835
        M_var.Choose_WCA[i + 2] = utils.lo_functions.Choose_WCA(
            Data.SFWMM_Daily_Outputs.at[i + 2, "RegWCA"],
            Pre_defined_Variables.Option_RegWCA,
            Pre_defined_Variables.Constant_RegWCA,
        )
        M_var.RegWCA[i + 2] = min(
            Pre_defined_Variables.MaxCap_RegWCA,
            Pre_defined_Variables.Multiplier_RegWCA * M_var.Choose_WCA[i + 2],
        )
        M_var.Choose_L8C51[i + 2] = utils.lo_functions.Choose_L8C51(
            Data.SFWMM_Daily_Outputs.at[i + 2, "RegL8C51"],
            Pre_defined_Variables.Option_RegL8C51,
            Pre_defined_Variables.Constant_RegL8C51,
        )
        M_var.RegL8C51[i + 2] = min(
            Pre_defined_Variables.MaxCap_RegL8C51,
            Pre_defined_Variables.Multiplier_RegL8C51 * M_var.Choose_L8C51[i + 2],
        )
        M_var.TotRegSo[i + 2] = (M_var.RegWCA[i + 2] + M_var.RegL8C51[i + 2]) * 1.9835
        M_var.Stage2ar[i + 2] = utils.stg_sto_ar.stg2ar(M_var.Lake_Stage[i + 1], 0)
        M_var.Stage2marsh[i + 2] = utils.stg_sto_ar.stg2mar(M_var.Lake_Stage[i + 1], 0)
        M_var.RF[i + 2] = Data.RF_Vol.at[i + 2, "RF_acft"]
        M_var.ET[i + 2] = utils.lo_functions.ET(
            Data.SFWMM_Daily_Outputs.at[i + 2, "et_dry"],
            M_var.Stage2ar[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "et_litoral"],
            M_var.Stage2marsh[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "et_open"],
            Data.ET_Vol.at[i + 2, "ETVol_acft"],
            Pre_defined_Variables.ET_Switch,
        )
        M_var.Choose_WSA_1[i + 2] = utils.lo_functions.Choose_WSA_1(
            df_WSMs.at[i + 2, "WSM1"],
            Pre_defined_Variables.Opt_WSA,
            Pre_defined_Variables.WSAtrig2,
            Pre_defined_Variables.WSAoff2,
        )
        M_var.Choose_WSA_2[i + 2] = utils.lo_functions.Choose_WSA_2(
            df_WSMs.at[i + 2, "WSM1"],
            Pre_defined_Variables.Opt_WSA,
            Pre_defined_Variables.WSAtrig1,
            Pre_defined_Variables.WSAoff1,
        )
        M_var.WSA_MIA[i + 2] = utils.lo_functions.WSA_MIA(
            WCA_Stages_df.at[i, "Are WCA stages too low?"],
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            M_var.Lake_Stage[i + 1],
            M_var.Choose_WSA_1[i + 2],
            Data.EAA_MIA_RUNOFF.at[i, "MIA"],
            Data.EAA_MIA_RUNOFF.at[i, "S3PMP"],
            M_var.Choose_WSA_2[i + 2],
            Pre_defined_Variables.Opt_WSA,
            Pre_defined_Variables.WSA_THC,
            Pre_defined_Variables.MIAcap2,
            Pre_defined_Variables.MIAcap1,
        )
        M_var.WSA_NNR[i + 2] = utils.lo_functions.WSA_NNR(
            WCA_Stages_df.at[i, "Are WCA stages too low?"],
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            M_var.Lake_Stage[i + 1],
            M_var.Choose_WSA_1[i + 2],
            Data.EAA_MIA_RUNOFF.at[i, "NNR"],
            Data.EAA_MIA_RUNOFF.at[i, "S2PMP"],
            M_var.Choose_WSA_2[i + 2],
            Pre_defined_Variables.Opt_WSA,
            Pre_defined_Variables.WSA_THC,
            Pre_defined_Variables.NNRcap2,
            Pre_defined_Variables.NNRcap1,
        )
        M_var.DSto[i + 2] = (
            M_var.NI_Supply[i + 2]
            + M_var.RF[i + 2]
            - M_var.ET[i + 2]
            + 1.9835
            * (
                M_var.Outlet2USBK[i + 2]
                + M_var.Outlet1USBK[i + 2]
                + M_var.WSA_MIA[i + 2]
                + M_var.WSA_NNR[i + 2]
                - M_var.Outlet1USEWS[i + 2]
            )
            - M_var.TotRegEW[i + 2]
            - M_var.TotRegSo[i + 2]
            + Storage_dev[i + 2]
        )
        M_var.Storage[i + 2] = utils.lo_functions.Storage(
            M_var.DayFlags[i + 2],
            M_var.Storage[i],
            StartStorage,
            M_var.Storage[i + 1],
            M_var.DSto[i + 2],
        )
        M_var.Lake_Stage[i + 2] = utils.lo_functions.Lake_Stage(
            utils.stg_sto_ar.stg2sto(M_var.Storage[i + 2], 1),
            Data.SFWMM_Daily_Outputs.at[i + 2, "EOD Stg(ft,NGVD)"],
            Pre_defined_Variables.Option_Stage,
        )
        # if M_var.Lake_Stage[i+2] >= 18:
        #################################################################################################################################################################
        Q_O[i + 2] = (
            (
                M_var.Outlet1USEWS[i + 2] * 0.028316847
                + ((M_var.TotRegEW[i + 2] + M_var.TotRegSo[i + 2]) / 70.0456)
            )
            * 3600
            * 24
        )

        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48  # m3/d
            Q_O_M[i] = Q_O[i]
            L_ext_M[i] = L_ext[i] + Q_I_M[i] * TP_Lake_N[i]
        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48  # m3/d
            Q_I_M[i] = Q_I[i]
            L_ext_M[i] = L_ext[i]
        Q_N2S[i] = (Q_I_M[i] + Q_O_M[i]) / 2
        M_var.Stage2ar[i + 2] = utils.stg_sto_ar.stg2ar(M_var.Lake_Stage[i + 2], 0)
        LO_WD[i] = M_var.Lake_Stage[i] * 0.3048 - LO_BL
        Lake_O_Storage_N[i] = (
            M_var.Storage[i] * TP_Variables.N_Per * 4046.85642 * 0.305
        )  # m3
        Lake_O_Storage_S[i] = (
            M_var.Storage[i] * TP_Variables.S_Per * 4046.85642 * 0.305
        )  # m3
        Lake_O_A_N[i] = M_var.Stage2ar[i] * TP_Variables.N_Per * 4046.85642  # m2
        Lake_O_A_S[i] = M_var.Stage2ar[i] * TP_Variables.S_Per * 4046.85642  # m2
        Lake_O_A_M_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Mud_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_S_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Sand_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_R_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Rock_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_P_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Peat_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_M_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Mud_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )
        Lake_O_A_S_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Sand_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )
        Lake_O_A_R_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Rock_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )
        Lake_O_A_P_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Peat_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )

        DIP_Lake_N[i] = TP_MBFR.DIP_Lake(TP_Lake_N[i])
        DIP_Lake_S[i] = TP_MBFR.DIP_Lake(TP_Lake_S[i])

        v_settle_N_c[i] = (R * g * d_c**2) / (
            C_1_c * nu_d[i] + (0.75 * C_2_c * R * g * d_c**3) ** 0.5
        )
        v_settle_N_s[i] = (R * g * d_s**2) / (
            C_1_s * nu_d[i] + (0.75 * C_2_s * R * g * d_s**3) ** 0.5
        )
        v_settle_N[i] = v_settle_N_c[i] * (
            (TP_Variables.A_Mud_N + TP_Variables.A_Peat_N) / TP_Variables.A_N
        ) + v_settle_N_s[i] * (
            (TP_Variables.A_Sand_N + TP_Variables.A_Rock_N) / TP_Variables.A_N
        )

        v_settle_S_c[i] = (R * g * d_c**2) / (
            C_1_c * nu_d[i] + (0.75 * C_2_c * R * g * d_c**3) ** 0.5
        )
        v_settle_S_s[i] = (R * g * d_s**2) / (
            C_1_s * nu_d[i] + (0.75 * C_2_s * R * g * d_s**3) ** 0.5
        )
        v_settle_S[i] = v_settle_S_c[i] * (
            (TP_Variables.A_Mud_S + TP_Variables.A_Peat_S) / TP_Variables.A_S
        ) + v_settle_S_s[i] * (
            (TP_Variables.A_Sand_S + TP_Variables.A_Rock_S) / TP_Variables.A_S
        )

        J_des_M_N[i] = TP_MBFR.Des_flux(Γ_M_N[i], Mass_sed_M_N, TP_Variables.K_des_M)
        J_des_S_N[i] = TP_MBFR.Des_flux(Γ_S_N[i], Mass_sed_S_N, TP_Variables.K_des_S)
        J_des_R_N[i] = TP_MBFR.Des_flux(Γ_R_N[i], Mass_sed_R_N, TP_Variables.K_des_R)
        J_des_P_N[i] = TP_MBFR.Des_flux(Γ_P_N[i], Mass_sed_P_N, TP_Variables.K_des_P)
        J_des_M_S[i] = TP_MBFR.Des_flux(Γ_M_S[i], Mass_sed_M_S, TP_Variables.K_des_M)
        J_des_S_S[i] = TP_MBFR.Des_flux(Γ_S_S[i], Mass_sed_S_S, TP_Variables.K_des_S)
        J_des_R_S[i] = TP_MBFR.Des_flux(Γ_R_S[i], Mass_sed_R_S, TP_Variables.K_des_R)
        J_des_P_S[i] = TP_MBFR.Des_flux(Γ_P_S[i], Mass_sed_P_S, TP_Variables.K_des_P)

        J_ads_M_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_M_N[i],
            Γ_M_N[i],
            Mass_sed_M_N,
            TP_Variables.K_ads_M,
            TP_Variables.Γ_inf,
        )
        J_ads_S_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_S_N[i],
            Γ_S_N[i],
            Mass_sed_S_N,
            TP_Variables.K_ads_S,
            TP_Variables.Γ_inf,
        )
        J_ads_R_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_R_N[i],
            Γ_R_N[i],
            Mass_sed_R_N,
            TP_Variables.K_ads_R,
            TP_Variables.Γ_inf,
        )
        J_ads_P_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_P_N[i],
            Γ_P_N[i],
            Mass_sed_P_N,
            TP_Variables.K_ads_P,
            TP_Variables.Γ_inf,
        )
        J_ads_M_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_M_S[i],
            Γ_M_S[i],
            Mass_sed_M_S,
            TP_Variables.K_ads_M,
            TP_Variables.Γ_inf,
        )
        J_ads_S_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_S_S[i],
            Γ_S_S[i],
            Mass_sed_S_S,
            TP_Variables.K_ads_S,
            TP_Variables.Γ_inf,
        )
        J_ads_R_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_R_S[i],
            Γ_R_S[i],
            Mass_sed_R_S,
            TP_Variables.K_ads_R,
            TP_Variables.Γ_inf,
        )
        J_ads_P_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_P_S[i],
            Γ_P_S[i],
            Mass_sed_P_S,
            TP_Variables.K_ads_P,
            TP_Variables.Γ_inf,
        )

        J_sedburial_M_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_M_N[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_N,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_sedburial_S_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_S_N[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_N,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_sedburial_R_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_R_N[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_N,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_sedburial_P_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_P_N[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_N,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )
        J_sedburial_M_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_M_S[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_S,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_sedburial_S_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_S_S[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_S,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_sedburial_R_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_R_S[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_S,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_sedburial_P_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_P_S[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_S,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )

        Sed_Resusp_M_N[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_M_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_S_N[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_S_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_R_N[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_R_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_P_N[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_P_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_M_S[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_M_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_S_S[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_S_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_R_S[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_R_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_P_S[i] = (
            ((E_0 / Td**E_1) * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2)
            * 10
            / LO_WD[i]
            * P_sed_P_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )

        P_sed_M_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_M_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.K_decomp_M,
                v_settle_N[i],
            )
            - Sed_Resusp_M_N[i] * Lake_O_Storage_N[i] / Mass_sed_M_N
            if TP_MBFR.P_sed(
                Lake_O_A_M_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.K_decomp_M,
                v_settle_N[i],
            )
            - Sed_Resusp_M_N[i] * Lake_O_Storage_N[i] / Mass_sed_M_N
            > 0
            else 0
        )
        P_sed_S_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_S_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.K_decomp_S,
                v_settle_N[i],
            )
            - Sed_Resusp_S_N[i] * Lake_O_Storage_N[i] / Mass_sed_S_N
            if TP_MBFR.P_sed(
                Lake_O_A_S_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.K_decomp_S,
                v_settle_N[i],
            )
            - Sed_Resusp_S_N[i] * Lake_O_Storage_N[i] / Mass_sed_S_N
            > 0
            else 0
        )
        P_sed_R_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_R_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.K_decomp_R,
                v_settle_N[i],
            )
            - Sed_Resusp_R_N[i] * Lake_O_Storage_N[i] / Mass_sed_R_N
            if TP_MBFR.P_sed(
                Lake_O_A_R_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.K_decomp_R,
                v_settle_N[i],
            )
            - Sed_Resusp_R_N[i] * Lake_O_Storage_N[i] / Mass_sed_R_N
            > 0
            else 0
        )
        P_sed_P_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_P_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.K_decomp_P,
                v_settle_N[i],
            )
            - Sed_Resusp_P_N[i] * Lake_O_Storage_N[i] / Mass_sed_P_N
            if TP_MBFR.P_sed(
                Lake_O_A_P_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.K_decomp_P,
                v_settle_N[i],
            )
            - Sed_Resusp_P_N[i] * Lake_O_Storage_N[i] / Mass_sed_P_N
            > 0
            else 0
        )
        P_sed_M_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_M_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.K_decomp_M,
                v_settle_S[i],
            )
            - Sed_Resusp_M_S[i] * Lake_O_Storage_S[i] / Mass_sed_M_S
            if TP_MBFR.P_sed(
                Lake_O_A_M_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.K_decomp_M,
                v_settle_S[i],
            )
            - Sed_Resusp_M_S[i] * Lake_O_Storage_S[i] / Mass_sed_M_S
            > 0
            else 0
        )
        P_sed_S_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_S_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.K_decomp_S,
                v_settle_S[i],
            )
            - Sed_Resusp_S_S[i] * Lake_O_Storage_S[i] / Mass_sed_S_S
            if TP_MBFR.P_sed(
                Lake_O_A_S_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.K_decomp_S,
                v_settle_S[i],
            )
            - Sed_Resusp_S_S[i] * Lake_O_Storage_S[i] / Mass_sed_S_S
            > 0
            else 0
        )
        P_sed_R_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_R_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.K_decomp_R,
                v_settle_S[i],
            )
            - Sed_Resusp_R_S[i] * Lake_O_Storage_S[i] / Mass_sed_R_S
            if TP_MBFR.P_sed(
                Lake_O_A_R_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.K_decomp_R,
                v_settle_S[i],
            )
            - Sed_Resusp_R_S[i] * Lake_O_Storage_S[i] / Mass_sed_R_S
            > 0
            else 0
        )
        P_sed_P_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_P_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.K_decomp_P,
                v_settle_S[i],
            )
            - Sed_Resusp_P_S[i] * Lake_O_Storage_S[i] / Mass_sed_P_S
            if TP_MBFR.P_sed(
                Lake_O_A_P_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.K_decomp_P,
                v_settle_S[i],
            )
            - Sed_Resusp_P_S[i] * Lake_O_Storage_S[i] / Mass_sed_P_S
            > 0
            else 0
        )

        J_Γburial_M_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_M_N[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_N,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_Γburial_S_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_S_N[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_N,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_Γburial_R_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_R_N[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_N,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_Γburial_P_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_P_N[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_N,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )
        J_Γburial_M_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_M_S[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_S,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_Γburial_S_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_S_S[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_S,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_Γburial_R_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_R_S[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_S,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_Γburial_P_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_P_S[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_S,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )

        Γ_M_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_M_N[i], J_des_M_N[i], J_Γburial_M_N[i], Γ_M_N[i], Mass_sed_M_N
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_M_N[i], J_des_M_N[i], J_Γburial_M_N[i], Γ_M_N[i], Mass_sed_M_N
            )
            > 0
            else 0
        )
        Γ_S_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_S_N[i], J_des_S_N[i], J_Γburial_S_N[i], Γ_S_N[i], Mass_sed_S_N
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_S_N[i], J_des_S_N[i], J_Γburial_S_N[i], Γ_S_N[i], Mass_sed_S_N
            )
            > 0
            else 0
        )
        Γ_R_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_R_N[i], J_des_R_N[i], J_Γburial_R_N[i], Γ_R_N[i], Mass_sed_R_N
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_R_N[i], J_des_R_N[i], J_Γburial_R_N[i], Γ_R_N[i], Mass_sed_R_N
            )
            > 0
            else 0
        )
        Γ_P_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_P_N[i], J_des_P_N[i], J_Γburial_P_N[i], Γ_P_N[i], Mass_sed_P_N
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_P_N[i], J_des_P_N[i], J_Γburial_P_N[i], Γ_P_N[i], Mass_sed_P_N
            )
            > 0
            else 0
        )
        Γ_M_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_M_S[i], J_des_M_S[i], J_Γburial_M_S[i], Γ_M_S[i], Mass_sed_M_S
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_M_S[i], J_des_M_S[i], J_Γburial_M_S[i], Γ_M_S[i], Mass_sed_M_S
            )
            > 0
            else 0
        )
        Γ_S_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_S_S[i], J_des_S_S[i], J_Γburial_S_S[i], Γ_S_S[i], Mass_sed_S_S
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_S_S[i], J_des_S_S[i], J_Γburial_S_S[i], Γ_S_S[i], Mass_sed_S_S
            )
            > 0
            else 0
        )
        Γ_R_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_R_S[i], J_des_R_S[i], J_Γburial_R_S[i], Γ_R_S[i], Mass_sed_R_S
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_R_S[i], J_des_R_S[i], J_Γburial_R_S[i], Γ_R_S[i], Mass_sed_R_S
            )
            > 0
            else 0
        )
        Γ_P_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_P_S[i], J_des_P_S[i], J_Γburial_P_S[i], Γ_P_S[i], Mass_sed_P_S
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_P_S[i], J_des_P_S[i], J_Γburial_P_S[i], Γ_P_S[i], Mass_sed_P_S
            )
            > 0
            else 0
        )

        J_decomp_M_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_M, P_sed_M_N[i], Mass_sed_M_N
        )
        J_decomp_S_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_S, P_sed_S_N[i], Mass_sed_S_N
        )
        J_decomp_R_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_R, P_sed_R_N[i], Mass_sed_R_N
        )
        J_decomp_P_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_P, P_sed_P_N[i], Mass_sed_P_N
        )
        J_decomp_M_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_M, P_sed_M_S[i], Mass_sed_M_S
        )
        J_decomp_S_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_S, P_sed_S_S[i], Mass_sed_S_S
        )
        J_decomp_R_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_R, P_sed_R_S[i], Mass_sed_R_S
        )
        J_decomp_P_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_P, P_sed_P_S[i], Mass_sed_P_S
        )

        DIP_pore_M_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_N[i],
                DIP_Lake_N[i],
                J_des_M_N[i],
                J_ads_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_N,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            if TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_N[i],
                DIP_Lake_N[i],
                J_des_M_N[i],
                J_ads_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_N,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            > 0
            else 0
        )
        DIP_pore_S_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_N[i],
                DIP_Lake_N[i],
                J_des_S_N[i],
                J_ads_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_N,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            if TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_N[i],
                DIP_Lake_N[i],
                J_des_S_N[i],
                J_ads_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_N,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            > 0
            else 0
        )
        DIP_pore_R_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_N[i],
                DIP_Lake_N[i],
                J_des_R_N[i],
                J_ads_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_N,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            if TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_N[i],
                DIP_Lake_N[i],
                J_des_R_N[i],
                J_ads_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_N,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            > 0
            else 0
        )
        DIP_pore_P_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                J_des_P_N[i],
                J_ads_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_N,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            if TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                J_des_P_N[i],
                J_ads_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_N,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            > 0
            else 0
        )
        DIP_pore_M_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_S[i],
                DIP_Lake_S[i],
                J_des_M_S[i],
                J_ads_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_S,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            if TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_S[i],
                DIP_Lake_S[i],
                J_des_M_S[i],
                J_ads_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_S,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            > 0
            else 0
        )
        DIP_pore_S_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_S[i],
                DIP_Lake_S[i],
                J_des_S_S[i],
                J_ads_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_S,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            if TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_S[i],
                DIP_Lake_S[i],
                J_des_S_S[i],
                J_ads_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_S,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            > 0
            else 0
        )
        DIP_pore_R_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_S[i],
                DIP_Lake_S[i],
                J_des_R_S[i],
                J_ads_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_S,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            if TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_S[i],
                DIP_Lake_S[i],
                J_des_R_S[i],
                J_ads_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_S,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            > 0
            else 0
        )
        DIP_pore_P_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                J_des_P_S[i],
                J_ads_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_S,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            if TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                J_des_P_S[i],
                J_ads_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_S,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            > 0
            else 0
        )

        Settling_P_N[i] = TP_MBFR.Sett_P(
            TP_Lake_N[i],
            DIP_Lake_N[i],
            Lake_O_A_N[i],
            Lake_O_Storage_N[i],
            v_settle_N[i],
        )
        Settling_P_S[i] = TP_MBFR.Sett_P(
            TP_Lake_S[i],
            DIP_Lake_S[i],
            Lake_O_A_S[i],
            Lake_O_Storage_S[i],
            v_settle_S[i],
        )

        P_diff_M_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_M,
            DIP_pore_M_N[i],
            DIP_Lake_N[i],
            Θ_M,
            TP_Variables.A_Mud_N,
            Lake_O_Storage_N[i],
        )
        P_diff_S_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_S,
            DIP_pore_S_N[i],
            DIP_Lake_N[i],
            Θ_S,
            TP_Variables.A_Sand_N,
            Lake_O_Storage_N[i],
        )
        P_diff_R_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_R,
            DIP_pore_R_N[i],
            DIP_Lake_N[i],
            Θ_R,
            TP_Variables.A_Rock_N,
            Lake_O_Storage_N[i],
        )
        P_diff_P_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_P,
            DIP_pore_P_N[i],
            DIP_Lake_N[i],
            Θ_P,
            TP_Variables.A_Peat_N,
            Lake_O_Storage_N[i],
        )
        P_diff_M_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_M,
            DIP_pore_M_S[i],
            DIP_Lake_S[i],
            Θ_M,
            TP_Variables.A_Mud_S,
            Lake_O_Storage_S[i],
        )
        P_diff_S_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_S,
            DIP_pore_S_S[i],
            DIP_Lake_S[i],
            Θ_S,
            TP_Variables.A_Sand_S,
            Lake_O_Storage_S[i],
        )
        P_diff_R_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_R,
            DIP_pore_R_S[i],
            DIP_Lake_S[i],
            Θ_R,
            TP_Variables.A_Rock_S,
            Lake_O_Storage_S[i],
        )
        P_diff_P_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_P,
            DIP_pore_P_S[i],
            DIP_Lake_S[i],
            Θ_P,
            TP_Variables.A_Peat_S,
            Lake_O_Storage_S[i],
        )

        # TP_N_to_S[i] = TP_MBFR.P_N_to_S(Q_N2S[i], TP_Lake_N[i], Lake_O_Storage_N[i])
        # TP_Out[i] = TP_MBFR.P_Out(Q_O_M[i], TP_Lake_S[i], Lake_O_Storage_S[i])

        TP_Lake_N[i + 1] = (
            TP_MBFR.TP_Lake_N(
                L_ext_M[i],
                Atm_Dep_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_N[i],
                DIP_pore_S_N[i],
                DIP_pore_R_N[i],
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                Q_N2S[i],
                Lake_O_A_N[i],
                TP_Lake_N[i],
                Lake_O_Storage_N[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_N[i],
            )
            + (
                Sed_Resusp_M_N[i]
                + Sed_Resusp_S_N[i]
                + Sed_Resusp_R_N[i]
                + Sed_Resusp_P_N[i]
            )
            if TP_MBFR.TP_Lake_N(
                L_ext_M[i],
                Atm_Dep_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_N[i],
                DIP_pore_S_N[i],
                DIP_pore_R_N[i],
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                Q_N2S[i],
                Lake_O_A_N[i],
                TP_Lake_N[i],
                Lake_O_Storage_N[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_N[i],
            )
            + (
                Sed_Resusp_M_N[i]
                + Sed_Resusp_S_N[i]
                + Sed_Resusp_R_N[i]
                + Sed_Resusp_P_N[i]
            )
            > 0
            else 0
        )
        TP_Lake_S[i + 1] = (
            TP_MBFR.TP_Lake_S(
                Atm_Dep_S[i],
                Q_N2S[i],
                TP_Lake_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_S[i],
                DIP_pore_S_S[i],
                DIP_pore_R_S[i],
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                Q_O_M[i],
                Lake_O_A_S[i],
                TP_Lake_S[i],
                Lake_O_Storage_S[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_S[i],
            )
            + (
                Sed_Resusp_M_S[i]
                + Sed_Resusp_S_S[i]
                + Sed_Resusp_R_S[i]
                + Sed_Resusp_P_S[i]
            )
            if TP_MBFR.TP_Lake_S(
                Atm_Dep_S[i],
                Q_N2S[i],
                TP_Lake_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_S[i],
                DIP_pore_S_S[i],
                DIP_pore_R_S[i],
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                Q_O_M[i],
                Lake_O_A_S[i],
                TP_Lake_S[i],
                Lake_O_Storage_S[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_S[i],
            )
            + (
                Sed_Resusp_M_S[i]
                + Sed_Resusp_S_S[i]
                + Sed_Resusp_R_S[i]
                + Sed_Resusp_P_S[i]
            )
            > 0
            else 0
        )
        TP_Lake_Mean[i + 1] = (TP_Lake_N[i + 1] + TP_Lake_S[i + 1]) / 2

        P_Load_Cal[i] = (
            M_var.Outlet1USREG[i] * 0.028316847 * 3600 * 24 * TP_Lake_S[i]
        )  # mg/d P
        P_Load_StL[i] = (
            M_var.Outlet2USRG[i] * 0.028316847 * 3600 * 24 * TP_Lake_S[i]
        )  # mg/d P
        P_Load_South[i] = M_var.TotRegSo[i] * 1233.48 * TP_Lake_S[i]  # mg/d P

    Output_df = pd.DataFrame(date_rng_2, columns=["Date"])  # 1/1/2008-12/31/2018

    Output_df["Stage_LO"] = M_var.Lake_Stage[2:]
    Output_df["S308_Q"] = M_var.Outlet2USRG[2:]
    Output_df["S77_Q"] = M_var.Outlet1USREG[2:]
    Output_df["Storage"] = M_var.Storage[2:]
    Output_df["Cut_back"] = M_var.Cut_back[2:]
    Output_df["P_Lake"] = TP_Lake_Mean
    # Output_df['P_Lake_N'] = TP_Lake_N
    # Output_df['P_Lake_S'] = TP_Lake_S
    # Output_df['DIP_pore_M_N'] = DIP_pore_M_N
    # Output_df['Q_N2S'] = Q_N2S
    # Output_df['Lake_O_A_N'] = Lake_O_A_N
    # Output_df['Lake_O_Storage_N'] = Lake_O_Storage_N
    # Output_df['Sed_Resusp_M_N'] = Sed_Resusp_M_N
    Output_df["P_Load_Cal"] = P_Load_Cal / 1e9  # tons
    Output_df["P_Load_StL"] = P_Load_StL / 1e9  # tons
    Output_df["P_Load_South"] = P_Load_South / 1e9  # tons

    return Output_df


Exported_File = LOONE_HydNut()
Exported_File.drop(index=Exported_File.index[-1], axis=0, inplace=True)
Exported_File.drop(index=Exported_File.index[-1], axis=0, inplace=True)
Exported_File["Stage_LO"] = Exported_File["Stage_LO"].astype(float)
Exported_File["Storage"] = Exported_File["Storage"].astype(float)
Exported_File["S308_Q"] = Exported_File["S308_Q"].astype(float)
Exported_File["S77_Q"] = Exported_File["S77_Q"].astype(float)
Exported_File["Cut_back"] = Exported_File["Cut_back"].astype(float)
Exported_File["P_Lake"] = pd.to_numeric(Exported_File["P_Lake"])
# Exported_File['P_Lake_N']=pd.to_numeric(Exported_File['P_Lake_N'])
# Exported_File['P_Lake_S']=pd.to_numeric(Exported_File['P_Lake_S'])
# Exported_File['DIP_pore_M_N']=pd.to_numeric(Exported_File['DIP_pore_M_N'])
# Exported_File['Q_N2S']=pd.to_numeric(Exported_File['Q_N2S'])
# Exported_File['Lake_O_A_N']=pd.to_numeric(Exported_File['Lake_O_A_N'])
# Exported_File['Lake_O_Storage_N']=pd.to_numeric(Exported_File['Lake_O_Storage_N'])
# Exported_File['Sed_Resusp_M_N']=pd.to_numeric(Exported_File['Sed_Resusp_M_N'])
Exported_File["P_Load_Cal"] = pd.to_numeric(Exported_File["P_Load_Cal"])
Exported_File["P_Load_StL"] = pd.to_numeric(Exported_File["P_Load_StL"])
Exported_File["P_Load_South"] = pd.to_numeric(Exported_File["P_Load_South"])

Exported_File = Exported_File.set_index("Date")
Exported_File.index = pd.to_datetime(Exported_File.index, unit="ns")
Exported_File_Mean = Exported_File.resample("M").mean()
Exported_File_Sum = Exported_File.resample("M").sum()

# Exported_File.to_csv('./Outputs/Daily.csv')
Exported_File_Mean.to_csv(
    "./Outputs/Exported_File_Opt_0809_Mean_%s.csv" % Pre_defined_Variables.Schedule
)
Exported_File_Sum.to_csv(
    "./Outputs/Exported_File_Opt_0809_Sum_%s.csv" % Pre_defined_Variables.Schedule
)
