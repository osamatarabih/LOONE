import os
import argparse
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from calendar import monthrange
from loone.utils import (
    load_config,
    replicate,
    df_wsms,
    stg_sto_ar,
    thc_class,
    lo_functions,
    dec_tree_functions,
    trib_hc,
)
from loone.data.model_variables import M_var as MVarClass
from loone.utils.wca_stages_class import WCA_Stages_Cls
from loone.data import Data as DClass


def LOONE_Q(workspace, p1, p2, s77_dv, s308_dv, tp_lake_s):
    """This function runs the LOONE Q module.

    Args:
        workspace (str): The path to the workspace directory.
        p1 (float): Parameter 1.
        p2 (float): Parameter 2.
        s77_dv (float): s77_dv value.
        s308_dv (float): s308_dv value.
        tp_lake_s (float): tp_lake_s value.

    Returns:
        None
    """
    os.chdir(workspace)
    config = load_config(workspace)

    Data = DClass(workspace)
    M_var = MVarClass(config)
    print("LOONE Q Module is Running!")
    # Based on the defined Start and End year, month, and day on the
    # Pre_defined_Variables File, Startdate and enddate are defined.
    year, month, day = map(int, config["start_date_entry"])
    startdate = datetime(year, month, day).date()
    year, month, day = map(int, config["start_date_entry"])
    begdateCS = datetime(year, month, day).date()
    year, month, day = map(int, config["end_date_entry"])
    enddate = datetime(year, month, day).date()

    ###################################################################
    if config["sim_type"] in [0, 1]:
        df_wsms.WSMs(workspace)

    df_WSMs = pd.read_csv("df_WSMs.csv")

    # The Following Code interpolates daily LOSA demand from weekly
    # data for 6 differnet datasets where the user defines the LOSA
    # demand that will be used based on a Code (1:6).
    # Set time frame for model run
    date_rng_2 = pd.date_range(start=startdate, end=enddate, freq="D")
    # Create a data frame with a date column
    Water_dmd = pd.DataFrame(date_rng_2, columns=["date"])

    N = []
    Wk = []
    # Generate a count list
    for i in Water_dmd["date"]:
        if i.month == startdate.month and i.day == startdate.day:
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
        D = ((Data.Weekly_dmd[f'C{config["code"]}'].iloc[i - 1]) / 7) * (
            config["multiplier"] / 100
        )
        dd.append(D)
    Water_dmd["Daily_demand"] = dd
    ###################################################################
    # Determine Tributary Hydrologic Conditions
    TC_LONINO_df = trib_hc.Trib_HC(workspace)
    # Determine WCA Stages
    WCA_Stages_df = WCA_Stages_Cls(workspace, TC_LONINO_df)
    # A dataframe to determine eachday's season (Months 11,12,1,2 are
    # Season 1, Months 3,4,5 are season 2, Months 6,7 are season 3,
    # Months 8,9,10 are season 4 )
    date_rng_5 = pd.date_range(start=startdate, end=enddate, freq="D")
    Seasons = pd.DataFrame(date_rng_5, columns=["date"])
    Seas_Count = len(Seasons.index)
    for i in range(Seas_Count):
        if (
            Seasons["date"].iloc[i].month > 2
            and Seasons["date"].iloc[i].month < 6
        ):
            S = 2
        elif (
            Seasons["date"].iloc[i].month > 5
            and Seasons["date"].iloc[i].month < 8
        ):
            S = 3
        elif (
            Seasons["date"].iloc[i].month > 7
            and Seasons["date"].iloc[i].month < 11
        ):
            S = 4
        else:
            S = 1
        M_var.Daily_Seasons[i] = S
        M_var.Mon[i] = Seasons["date"].iloc[i].month
    Seasons["Season"] = M_var.Daily_Seasons
    Seasons["Month"] = M_var.Mon
    ###################################################################
    # This following Script runs the main model daily simulations.
    date_rng_6 = pd.date_range(
        start=startdate - timedelta(days=1),
        end=enddate,
        freq="D",
    )
    LO_Model = pd.DataFrame(date_rng_6, columns=["date"])
    LO_Model["Net_Inflow"] = Data.NetInf_Input["Netflows_acft"]
    n_rows = len(LO_Model.index)
    LO_Model["LOSA_dmd_SFWMM"] = Data.SFWMM_W_dmd["LOSA_dmd"] * (
        config["mult_losa"] / 100
    )
    LO_Model["C44RO"] = Data.C44_Runoff["C44RO"]
    ##################################
    DecTree_df = pd.DataFrame(date_rng_5, columns=["Date"])
    DecTree_df["Zone_B_MetFcast"] = TC_LONINO_df["LONINO_Seasonal_Classes"]
    # Create a dataframe that includes Monthly Mean Basin Runoff &
    # BaseFlow-Runoff & Runoff-Baseflow (cfs)
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
            max(0, Data.C44RO["C44RO"].iloc[i] - Outlet2_baseflow)
            * Num_days[i]
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
    S80avgL1 = Data.Pulses["S-80_L1_LORS20082023"].mean()
    S80avgL2 = Data.Pulses["S-80_L2_LORS20082023"].mean()
    S80avgL3 = Data.Pulses["S-80_L3_LORS20082023"].mean()
    S77avgL1 = Data.Pulses["S-77_L1_LORS20082023"].mean()
    S77avgL2 = Data.Pulses["S-77_L2_LORS20082023"].mean()
    S77avgL3 = Data.Pulses["S-77_L3_LORS20082023"].mean()
    Basin_RO = Basin_RO.set_index(["date"])
    Basin_RO.index = pd.to_datetime(Basin_RO.index)
    Basin_RO_Daily = Basin_RO.reindex(date_rng_11d, method="ffill")
    Basin_RO = Basin_RO.reset_index()
    VLOOKUP1 = Basin_RO_Daily["BS-C44RO"]
    VLOOKUP1_c = [x for x in VLOOKUP1 if ~np.isnan(x)]
    ###################################################################
    # This following script contains the logic and calculations for
    # the proposed Lake Okeechobee Adaptive Protocol.
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
    AdapProt_df["Tributary Hydrologic Condition"] = TC_LONINO_df[
        "Tributary_Condition"
    ]
    # Define "Low Chance" 6/1 stg<11'
    if config["opt_date_targ_stg"] == 1:
        Targ_Stg = Data.Targ_Stg_June_1st
    else:
        Targ_Stg = Data.Targ_Stg_May_1st

    Targ_Stg_df = pd.DataFrame(date_rng_5, columns=["dates"])
    for i in range(len(Targ_Stg_df)):
        M_var.V10per[i] = replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            10,
            Targ_Stg,
        )
        M_var.V20per[i] = replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            20,
            Targ_Stg,
        )
        M_var.V25per[i] = replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            25,
            Targ_Stg,
        )
        M_var.V30per[i] = replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            30,
            Targ_Stg,
        )
        M_var.V40per[i] = replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            40,
            Targ_Stg,
        )
        M_var.V45per[i] = replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            45,
            Targ_Stg,
        )
        M_var.V50per[i] = replicate(
            Targ_Stg_df["dates"].iloc[i].year,
            Targ_Stg_df["dates"].iloc[i].timetuple().tm_yday,
            50,
            Targ_Stg,
        )
        M_var.V60per[i] = replicate(
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
    ###################################################################
    M_var.Lake_Stage[0] = config["beg_stage_cs"]
    M_var.Lake_Stage[1] = config["beg_stage_cs"]
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
    StartStorage = stg_sto_ar.stg2sto(config["start_stage"], 0)
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
    ###################################################################
    M_var.Zone_Code[0] = lo_functions.Zone_Code(
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
    M_var.LO_Zone[0] = lo_functions.LO_Zone(M_var.Zone_Code[0])
    for i in range(n_rows - 2):
        M_var.WSM_Zone[i + 2] = lo_functions.WSM_Zone(
            M_var.Lake_Stage[i + 1],
            df_WSMs.at[i + 1, "WSM4"],
            df_WSMs.at[i + 1, "WSM3"],
            df_WSMs.at[i + 1, "WSM2"],
            df_WSMs.at[i + 1, "WSM1"],
        )
        # Calculate Daily Maximum Water Supply
        # Note that in LOSA_dmd we used (i) because this file starts
        # from 1/1/2008 so i at this point =0.
        # Cutbacks are determined based on the WSM Zone.
        M_var.Max_Supply[i + 2] = lo_functions.Max_Supply(
            M_var.WSM_Zone[i + 2],
            Water_dmd.at[i, "Daily_demand"],
            config["z1_cutback"],
            config["z2_cutback"],
            config["z3_cutback"],
            config["z4_cutback"],
        )
        # Actual Daily Water supply
        M_var.LOSA_Supply[i + 2] = lo_functions.LOSA_Supply(
            M_var.WSM_Zone[i + 2],
            LO_Model.at[i + 2, "LOSA_dmd_SFWMM"],
            M_var.Max_Supply[i + 2],
            config["opt_losa_ws"],
        )
        # NetInflow - LOSA Supply
        M_var.NI_Supply[i + 2] = (
            LO_Model.at[i + 2, "Net_Inflow"] - M_var.LOSA_Supply[i + 2]
        )
        # TODO Note: for the pass statement, We will read the Daily
        # Water supply from the SFWMM as an input.
        # Calculate the cutback where Cutback = Demand - Supply
        ctbk = LO_Model.at[i + 2, "LOSA_dmd_SFWMM"] - M_var.LOSA_Supply[i + 2]
        M_var.Cut_back[i + 2] = ctbk
        # Calculate percentage of the demand that is not supplied for
        # each day
        if LO_Model.at[i + 2, "LOSA_dmd_SFWMM"] == 0:
            DNS = 0
        else:
            DNS = (
                M_var.Cut_back[i + 2] / LO_Model.at[i + 2, "LOSA_dmd_SFWMM"]
            ) * 100
        M_var.Dem_N_Sup[i + 2] = DNS
        # Calculate the Zone Code
        # Note that to calculate the Zone Code in Dec 31 2020 we
        # needed the WSM and breakpoint zones in 1/1/2021!
        # Note Also that i = 0 in Stage indicates Dec 30 1964 while
        # i = 0 in df_WSMs indicates Dec 31 1964!
        M_var.Zone_Code[i + 1] = lo_functions.Zone_Code(
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
        # Generate the Zone Column based on the corresponding Zone Code
        M_var.LO_Zone[i + 1] = lo_functions.LO_Zone(M_var.Zone_Code[i + 1])
        M_var.Zone_D_Trib[i] = dec_tree_functions.Zone_D_Trib(
            TC_LONINO_df.at[i, "Tributary_Condition"],
            config["opt_new_tree"],
        )
        M_var.Zone_D_stage[i] = dec_tree_functions.Zone_D_stage(
            M_var.Lake_Stage[i + 1], df_WSMs.at[i, "C-b"]
        )
        M_var.Zone_D_Seas[i] = dec_tree_functions.Zone_D_Seas(
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            M_var.Zone_D_Trib[i],
            config["opt_new_tree"],
        )
        M_var.Zone_D_MSeas[i] = dec_tree_functions.Zone_D_MSeas(
            TC_LONINO_df.at[i, "LONINO_MultiSeasonal_Classes"]
        )
        M_var.Zone_D_Branch_Code[i] = (
            M_var.Zone_D_Trib[i] * 1000
            + M_var.Zone_D_stage[i] * 100
            + M_var.Zone_D_Seas[i] * 10
            + M_var.Zone_D_MSeas[i] * 1
        )
        M_var.Zone_D_Rel_Code[i] = dec_tree_functions.Zone_D_Rel_Code(
            M_var.Zone_D_Branch_Code[i], config["opt_dec_tree"]
        )
        M_var.Zone_C_Trib[i] = dec_tree_functions.Zone_C_Trib(
            TC_LONINO_df.at[i, "Tributary_Condition"],
            config["opt_new_tree"],
        )
        M_var.Zone_C_Seas[i] = dec_tree_functions.Zone_C_Seas(
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            config["opt_new_tree"],
        )
        M_var.Zone_C_MSeas[i] = dec_tree_functions.Zone_C_MSeas(
            TC_LONINO_df.at[i, "LONINO_MultiSeasonal_Classes"]
        )
        M_var.Zone_C_MetFcast[i] = dec_tree_functions.Zone_C_MetFcast(
            M_var.Zone_C_Seas[i],
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            config["zone_c_met_fcast_indicator"],
        )
        M_var.Zone_C_Branch_Code[i] = (
            M_var.Zone_C_Trib[i] * 1000
            + M_var.Zone_C_MetFcast[i] * 100
            + M_var.Zone_C_Seas[i] * 10
            + M_var.Zone_C_MSeas[i] * 1
        )
        M_var.Zone_C_Rel_Code[i] = dec_tree_functions.Zone_C_Rel_Code(
            M_var.Zone_C_Branch_Code[i], config["opt_dec_tree"]
        )
        M_var.Zone_B_Trib[i] = dec_tree_functions.Zone_B_Trib(
            TC_LONINO_df.at[i, "Tributary_Condition"],
            config["opt_new_tree"],
        )
        M_var.Zone_B_Stage[i] = dec_tree_functions.Zone_B_Stage(
            M_var.Lake_Stage[i + 1], Seasons.at[i, "Season"]
        )
        M_var.Zone_B_Seas[i] = dec_tree_functions.Zone_B_Seas(
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"]
        )
        M_var.Zone_B_Branch_Code[i] = (
            M_var.Zone_B_Trib[i] * 1000
            + M_var.Zone_B_Stage[i] * 100
            + DecTree_df.at[i, "Zone_B_MetFcast"] * 10
            + M_var.Zone_B_Seas[i] * 1
        )
        M_var.Zone_B_Rel_Code[i] = dec_tree_functions.Zone_B_Rel_Code(
            M_var.Zone_B_Branch_Code[i], config["opt_dec_tree"]
        )
        M_var.DecTree_Relslevel[i + 2] = lo_functions.DecTree_Relslevel(
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
                    config["cs_flag"] == 0
                    or startdate.year == LO_Model.at[i, "date"].year
                )
            ):
                X2 = "SimDay1"
            else:
                X2 = LO_Model.at[i, "date"].date()
            M_var.DayFlags[i] = X2
        M_var.PlsDay[i + 2] = lo_functions.PlsDay(
            M_var.DayFlags[i + 2],
            M_var.DecTree_Relslevel[i + 2],
            config["pls_day_switch"],
        )
        M_var.Release_Level[i + 2] = lo_functions.Release_Level(
            M_var.Release_Level[i + 1],
            M_var.Lake_Stage[i + 1],
            TC_LONINO_df.at[i, "Tributary_Condition"],
            M_var.PlsDay[i + 2],
            M_var.Zone_Code[i + 1],
            M_var.DecTree_Relslevel[i + 2],
            config["max_qstg_trigger"],
        )
        if i >= 6:
            dh = M_var.Lake_Stage[i + 1] - M_var.Lake_Stage[i - 6]
            M_var.dh_7days[i + 1] = dh
        M_var.ZoneCodeminus1Code[i + 1] = lo_functions.ZoneCodeminus1Code(
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
        M_var.ZoneCodeCode[i + 1] = lo_functions.ZoneCodeCode(
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
            lo_functions.Fraction_of_Zone_height(
                M_var.Zone_Code[i + 1],
                M_var.Lake_Stage[i + 1],
                M_var.ZoneCodeminus1Code[i + 1],
                M_var.ZoneCodeCode[i + 1],
            )
        )
        M_var.ReLevelCode_1[i + 2] = lo_functions.ReLevelCode_1(
            M_var.Release_Level[i + 2],
            config["dstar_d1"],
            config["dstar_d2"],
            config["dstar_d3"],
            config["dstar_c"],
            config["dstar_b"],
        )
        M_var.ReLevelCode_2[i + 2] = lo_functions.ReLevelCode_2(
            M_var.Release_Level[i + 2],
            config["astar_d1"],
            config["astar_d2"],
            config["astar_d3"],
            config["astar_c"],
            config["astar_b"],
        )
        M_var.ReLevelCode_3_S80[i + 2] = lo_functions.ReLevelCode_3_S80(
            M_var.Release_Level[i + 2],
            config["bstar_s80_d1"],
            config["bstar_s80_d2"],
            config["bstar_s80_d3"],
            config["bstar_s80_c"],
            config["bstar_s80_b"],
        )
        M_var.Outlet2DS_Mult[i + 2] = lo_functions.Outlet2DS_Mult(
            Seasons.at[i, "Season"],
            Seasons.at[i, "Month"],
            M_var.dh_7days[i + 1],
            M_var.ReLevelCode_1[i + 2],
            M_var.Fraction_of_Zone_height[i + 1],
            M_var.ReLevelCode_2[i + 2],
            M_var.ReLevelCode_3_S80[i + 2],
            config["opt_qreg_mult"],
        )
        M_var.Outlet2DS_Mult_2[i + 2] = lo_functions.Outlet2DS_Mult_2(
            LO_Model.at[i + 2, "date"].month,
            LO_Model.at[i + 2, "date"].day,
            M_var.PlsDay[i + 2],
            M_var.Outlet2DS_Mult[i + 2 - M_var.PlsDay[i + 2]],
            M_var.Outlet2DS_Mult[i + 2],
            config["opt_qreg_mult"],
        )
        M_var.Outlet2DSRS[i + 2] = lo_functions.Outlet2DSRS(
            M_var.Release_Level[i + 2],
            Data.S80_RegRelRates.at[0, "Zone_D1"],
            S80avgL1,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                f'S-80_L1_{config["schedule"]}',
            ],
            M_var.Outlet2DS_Mult_2[i + 2],
            Data.CE_SLE_turns.at[
                LO_Model.at[i + 2, "date"].year - config["start_year"],
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
                f'S-80_L2_{config["schedule"]}',
            ],
            Data.S80_RegRelRates.at[0, "Zone_D3"],
            S80avgL3,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                f'S-80_L3_{config["schedule"]}',
            ],
            Data.S80_RegRelRates.at[0, "Zone_C"],
            Data.S80_RegRelRates.at[0, "Zone_B"],
            Data.S80_RegRelRates.at[0, "Zone_A"],
        )
        diff = M_var.Outlet2DSRS[i + 2] - LO_Model.at[i + 2, "C44RO"]
        M_var.Outlet2USRG1[i + 2] = max(
            0, diff.values[0] if isinstance(diff, pd.Series) else diff
        )
        M_var.Sum_Outlet2USRG1[i + 2] = lo_functions.Sum_Outlet2USRG1(
            LO_Model.at[i + 2, "date"].day, M_var.Outlet2USRG1[i + 2]
        )
        M_var.Outlet2DSBS[i + 2] = lo_functions.Outlet2DSBS(
            M_var.Release_Level[i + 2],
            M_var.Sum_Outlet2USRG1[i + 2],
            VLOOKUP1_c[i],
            Outlet2_baseflow,
            config["option_s80_baseflow"],
        )
        M_var.Outlet2USBK[i + 2] = lo_functions.Outlet2USBK(
            M_var.Lake_Stage[i + 1],
            df_WSMs.at[i + 1, "D1"],
            M_var.Outlet2USRG[i + 1],
            LO_Model.at[i + 2, "C44RO"],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S308BK"],
            config["opt_s308"],
            config["s308_bk_const"],
            config["s308_bk_thr"],
        )
        M_var.ROeast[i + 2] = (
            LO_Model.at[i + 2, "C44RO"] - M_var.Outlet2USBK[i + 2]
        )
        M_var.Outlet2USBS[i + 2] = lo_functions.Outlet2USBS(
            M_var.Outlet2DSBS[i + 2],
            M_var.Outlet2USRG1[i + 2],
            M_var.ROeast[i + 2],
            config["option_s80_baseflow"],
        )
        M_var.Sum_Outlet2USBK[i + 2] = lo_functions.Sum_Outlet2USBK(
            LO_Model.at[i + 2, "date"].day, M_var.Outlet2USBK[i + 2]
        )
        M_var.Outlet2USRG_Code[i + 2] = lo_functions.Outlet2USRG_Code(
            M_var.Outlet2USRG1[i + 2],
            M_var.Outlet2USBS[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S308RG"],
            Data.SFWMM_Daily_Outputs.at[i + 2, "STEST"],
            config["option_reg_s77_s308"],
        )
        if config["sim_type"] == 0:
            M_var.Outlet2USRG[i + 2] = lo_functions.Outlet2USRG(
                M_var.Outlet2USRG_Code[i + 2],
                Data.SFWMM_Daily_Outputs.at[i + 2, "S308RG"],
                Data.SFWMM_Daily_Outputs.at[i + 2, "STEST"],
                config["opt_s308"],
                config["s308_rg_const"],
            )
        else:
            if M_var.Lake_Stage[i + 1] >= 18:
                M_var.Outlet2USRG[i + 2] = 7200
            elif M_var.Lake_Stage[i + 1] <= 8:
                M_var.Outlet2USRG[i + 2] = 0
            elif (tp_lake_s[i] <= p1) and (
                date_rng_6[i + 2].month in [1, 2, 3, 4, 11, 12]
            ):
                M_var.Outlet2USRG[i + 2] = s308_dv[
                    (date_rng_6[i + 2].month) - 1
                ]
            elif (tp_lake_s[i] <= p2) and (
                date_rng_6[i + 2].month in [5, 6, 7, 8, 9, 10]
            ):
                M_var.Outlet2USRG[i + 2] = s308_dv[
                    (date_rng_6[i + 2].month) - 1
                ]
            else:
                M_var.Outlet2USRG[i + 2] = 0
        M_var.Outlet2DS[i + 2] = lo_functions.S80(
            M_var.ROeast[i + 2],
            M_var.Outlet2USRG[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S80"],
            config["s80_const"],
        )
        M_var.ReLevelCode_3_S77[i + 2] = lo_functions.ReLevelCode_3_S77(
            M_var.Release_Level[i + 2],
            config["bstar_s77_d1"],
            config["bstar_s77_d2"],
            config["bstar_s77_d3"],
            config["bstar_s77_c"],
            config["bstar_s77_b"],
        )
        M_var.Outlet1US_Mult[i + 2] = lo_functions.Outlet1US_Mult(
            Seasons.at[i, "Season"],
            Seasons.at[i, "Month"],
            M_var.dh_7days[i + 1],
            M_var.ReLevelCode_1[i + 2],
            M_var.Fraction_of_Zone_height[i + 1],
            M_var.ReLevelCode_2[i + 2],
            M_var.ReLevelCode_3_S77[i + 2],
            config["opt_qreg_mult"],
        )
        M_var.Outlet1US_Mult_2[i + 2] = lo_functions.Outlet1US_Mult_2(
            LO_Model.at[i + 2, "date"].month,
            LO_Model.at[i + 2, "date"].day,
            M_var.PlsDay[i + 2],
            M_var.Outlet1US_Mult[i + 2 - M_var.PlsDay[i + 2]],
            M_var.Outlet1US_Mult[i + 2],
            config["opt_qreg_mult"],
        )
        M_var.Outlet1USRS[i + 2] = lo_functions.Outlet1USRS(
            M_var.Release_Level[i + 2],
            Data.S77_RegRelRates.at[0, "Zone_D1"],
            S77avgL1,
            Data.Pulses.at[
                (
                    M_var.PlsDay[i + 2] - 1
                    if M_var.PlsDay[i + 2] - 1 >= 0
                    else len(Data.Pulses) - 1
                ),
                f'S-77_L1_{config["schedule"]}',
            ],
            M_var.Outlet1US_Mult_2[i + 2],
            LO_Model.at[i + 2, "C43RO"],
            Data.CE_SLE_turns.at[
                LO_Model.at[i + 2, "date"].year - config["start_year"],
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
                f'S-77_L2_{config["schedule"]}',
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
                f'S-77_L3_{config["schedule"]}',
            ],
            Data.S77_RegRelRates.at[0, "Zone_C"],
            Data.S77_RegRelRates.at[0, "Zone_B"],
            Data.S77_RegRelRates.at[0, "Zone_A"],
            config["opt_outlet1_dsrg"],
        )
        M_var.Sum_Outlet1USRS[i + 2] = lo_functions.Sum_Outlet1USRS(
            LO_Model.at[i + 2, "date"].day, M_var.Outlet1USRS[i + 2]
        )
        M_var.Outlet1USBK[i + 2] = lo_functions.Outlet1USBK(
            M_var.Lake_Stage[i + 1],
            M_var.Outlet1USRS[i + 2],
            M_var.Outlet1USBSAP[i + 1],
            M_var.Outlet1USEWS[i + 1],
            LO_Model.at[i + 2, "C43RO"],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S77BK"],
            config["outlet1_usbk_switch"],
            config["outlet1_usbk_threshold"],
        )
        M_var.ROwest[i + 2] = (
            LO_Model.at[i + 2, "C43RO"] - M_var.Outlet1USBK[i + 2]
        )
        M_var.Outlet1DSBS[i + 2] = lo_functions.Outlet1DSBS(
            M_var.Release_Level[i + 2],
            M_var.Sum_Outlet1USRS[i + 2],
            VLOOKUP2_c[i],
            Outlet1_baseflow,
            config["option_s77_baseflow"],
        )
        M_var.Outlet1USBS[i + 2] = lo_functions.Outlet1USBS(
            M_var.Outlet1DSBS[i + 2],
            M_var.Outlet1USRS[i + 2],
            M_var.ROwest[i + 2],
            config["option_s77_baseflow"],
        )
        # Define THC Class Normal or above
        if i < (n_rows - 2):
            M_var.Post_Ap_Baseflow[i] = thc_class.THC_Class(
                config,
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
            M_var.Post_AP_EWS[i] = thc_class.THC_Class(
                config,
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
        M_var.Outlet1USBSAP[i + 2] = lo_functions.Outlet1USBSAP(
            M_var.Outlet1USBS[i + 2],
            M_var.Post_Ap_Baseflow[i],
            config["opt_adap_prot"],
        )
        M_var.Outlet1USEWS[i + 2] = lo_functions.Outlet1USEWS(
            M_var.Post_AP_EWS[i],
            Data.SFWMM_Daily_Outputs.at[i + 2, "CAEST"],
            config["outlet1_usews_switch"],
            config["opt_adap_prot"],
        )
        if config["sim_type"] == 0:
            M_var.Outlet1USREG[i + 2] = lo_functions.Outlet1USREG(
                M_var.Outlet1USRS[i + 2],
                M_var.Outlet1USBSAP[i + 2],
                Data.SFWMM_Daily_Outputs.at[i + 2, "S77RG"],
                config["outlet1_usreg_switch"],
                config["option_reg_s77_s308"],
            )
        else:
            if M_var.Lake_Stage[i + 1] >= 18:
                M_var.Outlet1USREG[i + 2] = 7800
            elif M_var.Lake_Stage[i + 1] <= 8:
                M_var.Outlet1USREG[i + 2] = 0
            elif (tp_lake_s[i] <= p1) and (
                date_rng_6[i + 2].month in [1, 2, 3, 4, 11, 12]
            ):
                M_var.Outlet1USREG[i + 2] = s77_dv[
                    (date_rng_6[i + 2].month) - 1
                ]
            elif (tp_lake_s[i] <= p2) and (
                date_rng_6[i + 2].month in [5, 6, 7, 8, 9, 10]
            ):
                M_var.Outlet1USREG[i + 2] = s77_dv[
                    (date_rng_6[i + 2].month) - 1
                ]
            else:
                M_var.Outlet1USREG[i + 2] = 0
        M_var.Outlet1DS[i + 2] = lo_functions.Outlet1DS(
            M_var.Outlet1USREG[i + 2],
            M_var.Outlet1USEWS[i + 2],
            M_var.ROwest[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "S79"],
            config["outlet1_ds_switch"],
        )
        M_var.TotRegEW[i + 2] = (
            M_var.Outlet1USREG[i + 2] + M_var.Outlet2USRG[i + 2]
        ) * 1.9835
        M_var.Choose_WCA[i + 2] = lo_functions.Choose_WCA(
            Data.SFWMM_Daily_Outputs.at[i + 2, "RegWCA"],
            config["option_reg_wca"],
            config["constant_reg_wca"],
        )
        M_var.RegWCA[i + 2] = min(
            config["max_cap_reg_wca"],
            config["multiplier_reg_wca"] * M_var.Choose_WCA[i + 2],
        )
        M_var.Choose_L8C51[i + 2] = lo_functions.Choose_L8C51(
            Data.SFWMM_Daily_Outputs.at[i + 2, "RegL8C51"],
            config["option_reg_l8_c51"],
            config["constant_reg_l8_c51"],
        )
        M_var.RegL8C51[i + 2] = min(
            config["max_cap_reg_l8_c51"],
            config["multiplier_reg_l8_c51"] * M_var.Choose_L8C51[i + 2],
        )
        M_var.TotRegSo[i + 2] = (
            M_var.RegWCA[i + 2] + M_var.RegL8C51[i + 2]
        ) * 1.9835
        M_var.Stage2ar[i + 2] = stg_sto_ar.stg2ar(M_var.Lake_Stage[i + 1], 0)
        M_var.Stage2marsh[i + 2] = stg_sto_ar.stg2mar(
            M_var.Lake_Stage[i + 1], 0
        )
        M_var.RF[i + 2] = Data.RF_Vol.at[i + 2, "RFVol_acft"]
        M_var.ET[i + 2] = lo_functions.ET(
            Data.SFWMM_Daily_Outputs.at[i + 2, "et_dry"],
            M_var.Stage2ar[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "et_litoral"],
            M_var.Stage2marsh[i + 2],
            Data.SFWMM_Daily_Outputs.at[i + 2, "et_open"],
            Data.ET_Vol.at[i + 2, "ETVol_acft"],
            config["et_switch"],
        )
        M_var.Choose_WSA_1[i + 2] = lo_functions.Choose_WSA_1(
            df_WSMs.at[i + 2, "WSM1"],
            config["opt_wsa"],
            config["wsa_trig2"],
            config["wsa_off2"],
        )
        M_var.Choose_WSA_2[i + 2] = lo_functions.Choose_WSA_2(
            df_WSMs.at[i + 2, "WSM1"],
            config["opt_wsa"],
            config["wsa_trig1"],
            config["wsa_off1"],
        )
        M_var.WSA_MIA[i + 2] = lo_functions.WSA_MIA(
            WCA_Stages_df.at[i, "Are WCA stages too low?"],
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            M_var.Lake_Stage[i + 1],
            M_var.Choose_WSA_1[i + 2],
            Data.EAA_MIA_RUNOFF.at[i, "MIA"],
            Data.EAA_MIA_RUNOFF.at[i, "S3PMP"],
            M_var.Choose_WSA_2[i + 2],
            config["opt_wsa"],
            config["wsa_thc"],
            config["mia_cap2"],
            config["mia_cap1"],
        )
        M_var.WSA_NNR[i + 2] = lo_functions.WSA_NNR(
            WCA_Stages_df.at[i, "Are WCA stages too low?"],
            TC_LONINO_df.at[i, "LONINO_Seasonal_Classes"],
            M_var.Lake_Stage[i + 1],
            M_var.Choose_WSA_1[i + 2],
            Data.EAA_MIA_RUNOFF.at[i, "NNR"],
            Data.EAA_MIA_RUNOFF.at[i, "S2PMP"],
            M_var.Choose_WSA_2[i + 2],
            config["opt_wsa"],
            config["wsa_thc"],
            config["nnr_cap2"],
            config["nnr_cap1"],
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
        M_var.Storage[i + 2] = lo_functions.Storage(
            M_var.DayFlags[i + 2],
            M_var.Storage[i],
            StartStorage,
            M_var.Storage[i + 1],
            M_var.DSto[i + 2],
        )
        M_var.Lake_Stage[i + 2] = lo_functions.Lake_Stage(
            stg_sto_ar.stg2sto(M_var.Storage[i + 2], 1),
            Data.SFWMM_Daily_Outputs.at[i + 2, "EOD Stg(ft,NGVD)"],
            config["option_stage"],
        )
        # if M_var.Lake_Stage[i+2] >= 18:
        #     Flood[i+2] = 1
        # else:
        #     Flood[i+2] = 0
        LO_Model["Cutback"] = M_var.Cut_back
        LO_Model["Stage"] = M_var.Lake_Stage
        LO_Model["Storage"] = M_var.Storage
        LO_Model["S77_Q"] = M_var.Outlet1USREG
        LO_Model["S308_Q"] = M_var.Outlet2USRG
        LO_Model["S77EW"] = M_var.Outlet1USEWS
        LO_Model["TotRegEW"] = M_var.TotRegEW
        LO_Model["TotRegSo"] = M_var.TotRegSo

    # Write out the results to a file - Needed because
    # execute_loone.py calls this script as a subprocess and can't get
    # this data.
    LO_Model.to_csv("LOONE_Q_Outputs.csv")

    return [LO_Model]


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "workspace",
        nargs=1,
        help="The path to the working directory.",
    )
    argparser.add_argument("--p1", nargs=1)
    argparser.add_argument("--p2", nargs=1)
    argparser.add_argument("--s77_dv", nargs=1)
    argparser.add_argument("--s308_dv", nargs=1)
    argparser.add_argument("--tp_lake_s", nargs=1)
    args = argparser.parse_args()
    workspace = args.workspace[0]
    p1 = args.p1[0] if args.p1 else 0
    p2 = args.p2[0] if args.p2 else 0
    s77_dv = args.s77_dv[0] if args.s77_dv else 0
    s308_dv = args.s308_dv[0] if args.s308_dv else 0
    tp_lake_s = args.tp_lake_s[0] if args.tp_lake_s else 0
    LOONE_Q(workspace, p1, p2, s77_dv, s308_dv, tp_lake_s)
