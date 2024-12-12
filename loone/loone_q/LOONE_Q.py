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

    data = DClass(workspace)
    model_variables = MVarClass(config)
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
    date_range_2 = pd.date_range(start=startdate, end=enddate, freq="D")
    # Create a data frame with a date column
    water_demand = pd.DataFrame(date_range_2, columns=["date"])

    N = []
    week_number = []
    # Generate a count list
    for i in water_demand["date"]:
        if i.month == startdate.month and i.day == startdate.day:
            n = 0
        else:
            n = n + 1
        N.append(n)
    water_demand["count"] = N
    # Calculate the week number for all rows in the data frame
    for i in water_demand["count"]:
        if i > 363:
            J = 52
        else:
            J = int(i / 7) + 1
        week_number.append(J)
    water_demand["Week_num"] = week_number
    daily_demand = []  # daily demand
    # Calculate daily water demand
    for i in water_demand["Week_num"]:
        demand = ((data.Weekly_dmd[f'C{config["code"]}'].iloc[i - 1]) / 7) * (
            config["multiplier"] / 100
        )
        daily_demand.append(demand)
    water_demand["Daily_demand"] = daily_demand
    ###################################################################
    # Determine Tributary Hydrologic Conditions
    tc_lonino_df = trib_hc.Trib_HC(workspace)
    # Determine WCA Stages
    wca_stages_df = WCA_Stages_Cls(workspace, tc_lonino_df)
    # A dataframe to determine eachday's season (Months 11,12,1,2 are
    # Season 1, Months 3,4,5 are season 2, Months 6,7 are season 3,
    # Months 8,9,10 are season 4 )
    date_range_5 = pd.date_range(start=startdate, end=enddate, freq="D")
    seasons = pd.DataFrame(date_range_5, columns=["date"])
    seasons_count = len(seasons.index)
    for i in range(seasons_count):
        if (
            seasons["date"].iloc[i].month > 2
            and seasons["date"].iloc[i].month < 6
        ):
            S = 2
        elif (
            seasons["date"].iloc[i].month > 5
            and seasons["date"].iloc[i].month < 8
        ):
            S = 3
        elif (
            seasons["date"].iloc[i].month > 7
            and seasons["date"].iloc[i].month < 11
        ):
            S = 4
        else:
            S = 1
        model_variables.Daily_Seasons[i] = S
        model_variables.Mon[i] = seasons["date"].iloc[i].month
    seasons["Season"] = model_variables.Daily_Seasons
    seasons["Month"] = model_variables.Mon
    ###################################################################
    # This following Script runs the main model daily simulations.
    date_range_6 = pd.date_range(
        start=startdate - timedelta(days=1),
        end=enddate,
        freq="D",
    )
    lo_model = pd.DataFrame(date_range_6, columns=["date"])
    lo_model["Net_Inflow"] = data.NetInf_Input["Netflows_acft"]
    n_rows = len(lo_model.index)
    lo_model["LOSA_dmd_SFWMM"] = data.SFWMM_W_dmd["LOSA_dmd"] * (
        config["mult_losa"] / 100
    )
    lo_model["C44RO"] = data.C44_Runoff["C44RO"]
    ##################################
    dec_tree_df = pd.DataFrame(date_range_5, columns=["Date"])
    dec_tree_df["Zone_B_MetFcast"] = tc_lonino_df["LONINO_Seasonal_Classes"]
    # Create a dataframe that includes Monthly Mean Basin Runoff &
    # BaseFlow-Runoff & Runoff-Baseflow (cfs)
    date_range_11 = pd.date_range(start=startdate, end=enddate, freq="MS")
    date_range_11d = pd.date_range(start=startdate, end=enddate, freq="D")
    date_range_11d.name = "Date"
    basin_ro = pd.DataFrame(date_range_11, columns=["date"])
    # Baseflows
    outlet1_baseflow = data.S77_RegRelRates["Zone_D0"].iloc[0]
    outlet2_baseflow = data.S80_RegRelRates["Zone_D0"].iloc[0]
    # Calculta number of months in the timeseries data.
    num_B_R = len(basin_ro.index)
    bs_c43ro = np.zeros(num_B_R)
    bs_c44ro = np.zeros(num_B_R)
    c44ro_sltrib = np.zeros(num_B_R)
    c44ro_bs = np.zeros(num_B_R)
    num_days = np.zeros(num_B_R)
    for i in range(num_B_R):
        num_days[i] = monthrange(
            basin_ro["date"].iloc[i].year, basin_ro["date"].iloc[i].month
        )[
            1
        ]  # no. of days in each time step month.
        bs_c43ro[i] = max(0, (outlet1_baseflow - data.C43RO["C43RO"].iloc[i]))
        bs_c44ro[i] = max(0, (outlet2_baseflow - data.C44RO["C44RO"].iloc[i]))
        c44ro_sltrib[i] = bs_c44ro[i] + data.SLTRIB["SLTRIB_cfs"].iloc[i]
        c44ro_bs[i] = (
            max(0, data.C44RO["C44RO"].iloc[i] - outlet2_baseflow)
            * num_days[i]
        )
    basin_ro["Ndays"] = num_days
    basin_ro["C43RO"] = data.C43RO["C43RO"]
    basin_ro["BS-C43RO"] = bs_c43ro
    basin_ro["C44RO"] = data.C44RO["C44RO"]
    basin_ro["BS-C44RO"] = bs_c44ro
    basin_ro["SLTRIB"] = data.SLTRIB["SLTRIB_cfs"]
    basin_ro["C44RO_SLTRIB"] = c44ro_sltrib
    basin_ro["C44RO-BS"] = c44ro_bs
    lo_model["C43RO"] = data.C43RO_Daily["C43RO"]
    s80avg_l1 = data.Pulses["S-80_L1_LORS20082023"].mean()
    s80avg_l2 = data.Pulses["S-80_L2_LORS20082023"].mean()
    s80avg_l3 = data.Pulses["S-80_L3_LORS20082023"].mean()
    s77avg_l1 = data.Pulses["S-77_L1_LORS20082023"].mean()
    s77avg_l2 = data.Pulses["S-77_L2_LORS20082023"].mean()
    s77avg_l3 = data.Pulses["S-77_L3_LORS20082023"].mean()
    basin_ro = basin_ro.set_index(["date"])
    basin_ro.index = pd.to_datetime(basin_ro.index)
    basin_ro_daily = basin_ro.reindex(date_range_11d, method="ffill")
    basin_ro = basin_ro.reset_index()
    vlookup1 = basin_ro_daily["BS-C44RO"]
    vlookup1_c = [x for x in vlookup1 if ~np.isnan(x)]
    ###################################################################
    # This following script contains the logic and calculations for
    # the proposed Lake Okeechobee Adaptive Protocol.
    adaptive_protocol_df = pd.DataFrame(date_range_5, columns=["date"])
    # Calculate Late Dry Season (Apr-May) logic.
    late_dry_season = []
    for i in adaptive_protocol_df["date"]:
        if i.month > 3 and i.month < 6:
            L = True
        else:
            L = False
        late_dry_season.append(L)
    adaptive_protocol_df["Late_Dry_Season"] = late_dry_season
    adaptive_protocol_df["Tributary Hydrologic Condition"] = tc_lonino_df[
        "Tributary_Condition"
    ]
    # Define "Low Chance" 6/1 stg<11'
    if config["opt_date_targ_stg"] == 1:
        targ_stg = data.Targ_Stg_June_1st
    else:
        targ_stg = data.Targ_Stg_May_1st

    targ_stg_df = pd.DataFrame(date_range_5, columns=["dates"])
    for i in range(len(targ_stg_df)):
        model_variables.V10per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            10,
            targ_stg,
        )
        model_variables.V20per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            20,
            targ_stg,
        )
        model_variables.V25per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            25,
            targ_stg,
        )
        model_variables.V30per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            30,
            targ_stg,
        )
        model_variables.V40per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            40,
            targ_stg,
        )
        model_variables.V45per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            45,
            targ_stg,
        )
        model_variables.V50per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            50,
            targ_stg,
        )
        model_variables.V60per[i] = replicate(
            targ_stg_df["dates"].iloc[i].year,
            targ_stg_df["dates"].iloc[i].timetuple().tm_yday,
            60,
            targ_stg,
        )

    targ_stg_df["10%"] = [x for x in model_variables.V10per if ~np.isnan(x)]
    targ_stg_df["20%"] = [x for x in model_variables.V20per if ~np.isnan(x)]
    targ_stg_df["25%"] = [x for x in model_variables.V25per if ~np.isnan(x)]
    targ_stg_df["30%"] = [x for x in model_variables.V30per if ~np.isnan(x)]
    targ_stg_df["40%"] = [x for x in model_variables.V40per if ~np.isnan(x)]
    targ_stg_df["45%"] = [x for x in model_variables.V45per if ~np.isnan(x)]
    targ_stg_df["50%"] = [x for x in model_variables.V50per if ~np.isnan(x)]
    targ_stg_df["60%"] = [x for x in model_variables.V60per if ~np.isnan(x)]

    # Outlet1_baseflow = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    outlet1_baseflow = 450  # cfs
    vlookup2 = basin_ro_daily["BS-C43RO"]
    vlookup2_c = [x for x in vlookup2 if ~np.isnan(x)]
    ###################################################################
    model_variables.Lake_Stage[0] = config["beg_stage_cs"]
    model_variables.Lake_Stage[1] = config["beg_stage_cs"]
    model_variables.DecTree_Relslevel[0] = np.nan
    model_variables.DecTree_Relslevel[1] = np.nan
    if (
        startdate.month == lo_model["date"].iloc[2].month
        and startdate.day == lo_model["date"].iloc[2].day
    ):
        x1 = "SimDay1"
    elif (
        begdateCS.year == lo_model["date"].iloc[2].year
        and begdateCS.month == lo_model["date"].iloc[2].month
        and begdateCS.day == lo_model["date"].iloc[2].day
    ):
        x1 = "CS start date"
    else:
        x1 = lo_model["date"].iloc[2]
    model_variables.DayFlags[2] = x1
    start_storage = stg_sto_ar.stg2sto(config["start_stage"], 0)
    model_variables.Storage[0] = start_storage
    model_variables.Storage[1] = start_storage
    # Flood = np.zeros(n_rows, dtype = object)
    ##Here, I will insert the Storage Deviaiton Values as Input!
    storage_deviation = data.Storage_dev_df["DS_dev"]
    # Create a Choose Function for AP Post Baseflow
    # if Pre_defined_Variables.Opt_AdapProt == 0:
    #     C = 450
    # elif Pre_defined_Variables.Opt_AdapProt == 1:
    #     C = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    # Choose_1 = C
    choose_1 = 450  # cfs
    ###################################################################
    model_variables.Zone_Code[0] = lo_functions.Zone_Code(
        model_variables.Lake_Stage[0],
        df_WSMs["A"].iloc[0],
        df_WSMs["B"].iloc[0],
        df_WSMs["C"].iloc[0],
        df_WSMs["D3"].iloc[0],
        df_WSMs["D2"].iloc[0],
        df_WSMs["D1"].iloc[0],
        df_WSMs["D0"].iloc[0],
        df_WSMs["WSM1"].iloc[0],
    )
    model_variables.LO_Zone[0] = lo_functions.LO_Zone(model_variables.Zone_Code[0])
    for i in range(n_rows - 2):
        model_variables.WSM_Zone[i + 2] = lo_functions.WSM_Zone(
            model_variables.Lake_Stage[i + 1],
            df_WSMs.at[i + 1, "WSM4"],
            df_WSMs.at[i + 1, "WSM3"],
            df_WSMs.at[i + 1, "WSM2"],
            df_WSMs.at[i + 1, "WSM1"],
        )
        # Calculate Daily Maximum Water Supply
        # Note that in LOSA_dmd we used (i) because this file starts
        # from 1/1/2008 so i at this point =0.
        # Cutbacks are determined based on the WSM Zone.
        model_variables.Max_Supply[i + 2] = lo_functions.Max_Supply(
            model_variables.WSM_Zone[i + 2],
            water_demand.at[i, "Daily_demand"],
            config["z1_cutback"],
            config["z2_cutback"],
            config["z3_cutback"],
            config["z4_cutback"],
        )
        # Actual Daily Water supply
        model_variables.LOSA_Supply[i + 2] = lo_functions.LOSA_Supply(
            model_variables.WSM_Zone[i + 2],
            lo_model.at[i + 2, "LOSA_dmd_SFWMM"],
            model_variables.Max_Supply[i + 2],
            config["opt_losa_ws"],
        )
        # NetInflow - LOSA Supply
        model_variables.NI_Supply[i + 2] = (
            lo_model.at[i + 2, "Net_Inflow"] - model_variables.LOSA_Supply[i + 2]
        )
        # TODO Note: for the pass statement, We will read the Daily
        # Water supply from the SFWMM as an input.
        # Calculate the cutback where Cutback = Demand - Supply
        cutback = lo_model.at[i + 2, "LOSA_dmd_SFWMM"] - model_variables.LOSA_Supply[i + 2]
        model_variables.Cut_back[i + 2] = cutback
        # Calculate percentage of the demand that is not supplied for
        # each day
        if lo_model.at[i + 2, "LOSA_dmd_SFWMM"] == 0:
            demand_not_supplied = 0
        else:
            demand_not_supplied = (
                model_variables.Cut_back[i + 2] / lo_model.at[i + 2, "LOSA_dmd_SFWMM"]
            ) * 100
        model_variables.Dem_N_Sup[i + 2] = demand_not_supplied
        # Calculate the Zone Code
        # Note that to calculate the Zone Code in Dec 31 2020 we
        # needed the WSM and breakpoint zones in 1/1/2021!
        # Note Also that i = 0 in Stage indicates Dec 30 1964 while
        # i = 0 in df_WSMs indicates Dec 31 1964!
        model_variables.Zone_Code[i + 1] = lo_functions.Zone_Code(
            model_variables.Lake_Stage[i + 1],
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
        model_variables.LO_Zone[i + 1] = lo_functions.LO_Zone(model_variables.Zone_Code[i + 1])
        model_variables.Zone_D_Trib[i] = dec_tree_functions.Zone_D_Trib(
            tc_lonino_df.at[i, "Tributary_Condition"],
            config["opt_new_tree"],
        )
        model_variables.Zone_D_stage[i] = dec_tree_functions.Zone_D_stage(
            model_variables.Lake_Stage[i + 1], df_WSMs.at[i, "C-b"]
        )
        model_variables.Zone_D_Seas[i] = dec_tree_functions.Zone_D_Seas(
            tc_lonino_df.at[i, "LONINO_Seasonal_Classes"],
            model_variables.Zone_D_Trib[i],
            config["opt_new_tree"],
        )
        model_variables.Zone_D_MSeas[i] = dec_tree_functions.Zone_D_MSeas(
            tc_lonino_df.at[i, "LONINO_MultiSeasonal_Classes"]
        )
        model_variables.Zone_D_Branch_Code[i] = (
            model_variables.Zone_D_Trib[i] * 1000
            + model_variables.Zone_D_stage[i] * 100
            + model_variables.Zone_D_Seas[i] * 10
            + model_variables.Zone_D_MSeas[i] * 1
        )
        model_variables.Zone_D_Rel_Code[i] = dec_tree_functions.Zone_D_Rel_Code(
            model_variables.Zone_D_Branch_Code[i], config["opt_dec_tree"]
        )
        model_variables.Zone_C_Trib[i] = dec_tree_functions.Zone_C_Trib(
            tc_lonino_df.at[i, "Tributary_Condition"],
            config["opt_new_tree"],
        )
        model_variables.Zone_C_Seas[i] = dec_tree_functions.Zone_C_Seas(
            tc_lonino_df.at[i, "LONINO_Seasonal_Classes"],
            config["opt_new_tree"],
        )
        model_variables.Zone_C_MSeas[i] = dec_tree_functions.Zone_C_MSeas(
            tc_lonino_df.at[i, "LONINO_MultiSeasonal_Classes"]
        )
        model_variables.Zone_C_MetFcast[i] = dec_tree_functions.Zone_C_MetFcast(
            model_variables.Zone_C_Seas[i],
            tc_lonino_df.at[i, "LONINO_Seasonal_Classes"],
            config["zone_c_met_fcast_indicator"],
        )
        model_variables.Zone_C_Branch_Code[i] = (
            model_variables.Zone_C_Trib[i] * 1000
            + model_variables.Zone_C_MetFcast[i] * 100
            + model_variables.Zone_C_Seas[i] * 10
            + model_variables.Zone_C_MSeas[i] * 1
        )
        model_variables.Zone_C_Rel_Code[i] = dec_tree_functions.Zone_C_Rel_Code(
            model_variables.Zone_C_Branch_Code[i], config["opt_dec_tree"]
        )
        model_variables.Zone_B_Trib[i] = dec_tree_functions.Zone_B_Trib(
            tc_lonino_df.at[i, "Tributary_Condition"],
            config["opt_new_tree"],
        )
        model_variables.Zone_B_Stage[i] = dec_tree_functions.Zone_B_Stage(
            model_variables.Lake_Stage[i + 1], seasons.at[i, "Season"]
        )
        model_variables.Zone_B_Seas[i] = dec_tree_functions.Zone_B_Seas(
            tc_lonino_df.at[i, "LONINO_Seasonal_Classes"]
        )
        model_variables.Zone_B_Branch_Code[i] = (
            model_variables.Zone_B_Trib[i] * 1000
            + model_variables.Zone_B_Stage[i] * 100
            + dec_tree_df.at[i, "Zone_B_MetFcast"] * 10
            + model_variables.Zone_B_Seas[i] * 1
        )
        model_variables.Zone_B_Rel_Code[i] = dec_tree_functions.Zone_B_Rel_Code(
            model_variables.Zone_B_Branch_Code[i], config["opt_dec_tree"]
        )
        model_variables.DecTree_Relslevel[i + 2] = lo_functions.DecTree_Relslevel(
            model_variables.Zone_Code[i + 1],
            model_variables.Zone_D_Rel_Code[i],
            model_variables.Zone_C_Rel_Code[i],
            model_variables.Zone_B_Rel_Code[i],
        )
        if i >= 3:
            if (
                startdate.month == lo_model.at[i, "date"].month
                and startdate.day == lo_model.at[i, "date"].day
                and (
                    config["cs_flag"] == 0
                    or startdate.year == lo_model.at[i, "date"].year
                )
            ):
                X2 = "SimDay1"
            else:
                X2 = lo_model.at[i, "date"].date()
            model_variables.DayFlags[i] = X2
        model_variables.PlsDay[i + 2] = lo_functions.PlsDay(
            model_variables.DayFlags[i + 2],
            model_variables.DecTree_Relslevel[i + 2],
            config["pls_day_switch"],
        )
        model_variables.Release_Level[i + 2] = lo_functions.release_level(
            model_variables.Release_Level[i + 1],
            model_variables.Lake_Stage[i + 1],
            tc_lonino_df.at[i, "Tributary_Condition"],
            model_variables.PlsDay[i + 2],
            model_variables.Zone_Code[i + 1],
            model_variables.DecTree_Relslevel[i + 2],
            config["max_qstg_trigger"],
        )
        if i >= 6:
            dh = model_variables.Lake_Stage[i + 1] - model_variables.Lake_Stage[i - 6]
            model_variables.dh_7days[i + 1] = dh
        model_variables.ZoneCodeminus1Code[i + 1] = lo_functions.ZoneCodeminus1Code(
            model_variables.Zone_Code[i + 1],
            df_WSMs.at[i + 1, "WSM1"],
            df_WSMs.at[i + 1, "D0"],
            df_WSMs.at[i + 1, "D1"],
            df_WSMs.at[i + 1, "D2"],
            df_WSMs.at[i + 1, "D3"],
            df_WSMs.at[i + 1, "C"],
            df_WSMs.at[i + 1, "B"],
            df_WSMs.at[i + 1, "A"],
        )
        model_variables.ZoneCodeCode[i + 1] = lo_functions.ZoneCodeCode(
            model_variables.Zone_Code[i + 1],
            df_WSMs.at[i + 1, "WSM1"],
            df_WSMs.at[i + 1, "D0"],
            df_WSMs.at[i + 1, "D1"],
            df_WSMs.at[i + 1, "D2"],
            df_WSMs.at[i + 1, "D3"],
            df_WSMs.at[i + 1, "C"],
            df_WSMs.at[i + 1, "B"],
            df_WSMs.at[i + 1, "A"],
        )
        model_variables.Fraction_of_Zone_height[i + 1] = (
            lo_functions.Fraction_of_Zone_height(
                model_variables.Zone_Code[i + 1],
                model_variables.Lake_Stage[i + 1],
                model_variables.ZoneCodeminus1Code[i + 1],
                model_variables.ZoneCodeCode[i + 1],
            )
        )
        model_variables.ReLevelCode_1[i + 2] = lo_functions.ReLevelCode_1(
            model_variables.Release_Level[i + 2],
            config["dstar_d1"],
            config["dstar_d2"],
            config["dstar_d3"],
            config["dstar_c"],
            config["dstar_b"],
        )
        model_variables.ReLevelCode_2[i + 2] = lo_functions.ReLevelCode_2(
            model_variables.Release_Level[i + 2],
            config["astar_d1"],
            config["astar_d2"],
            config["astar_d3"],
            config["astar_c"],
            config["astar_b"],
        )
        model_variables.ReLevelCode_3_S80[i + 2] = lo_functions.ReLevelCode_3_S80(
            model_variables.Release_Level[i + 2],
            config["bstar_s80_d1"],
            config["bstar_s80_d2"],
            config["bstar_s80_d3"],
            config["bstar_s80_c"],
            config["bstar_s80_b"],
        )
        model_variables.Outlet2DS_Mult[i + 2] = lo_functions.Outlet2DS_Mult(
            seasons.at[i, "Season"],
            seasons.at[i, "Month"],
            model_variables.dh_7days[i + 1],
            model_variables.ReLevelCode_1[i + 2],
            model_variables.Fraction_of_Zone_height[i + 1],
            model_variables.ReLevelCode_2[i + 2],
            model_variables.ReLevelCode_3_S80[i + 2],
            config["opt_qreg_mult"],
        )
        model_variables.Outlet2DS_Mult_2[i + 2] = lo_functions.Outlet2DS_Mult_2(
            lo_model.at[i + 2, "date"].month,
            lo_model.at[i + 2, "date"].day,
            model_variables.PlsDay[i + 2],
            model_variables.Outlet2DS_Mult[i + 2 - model_variables.PlsDay[i + 2]],
            model_variables.Outlet2DS_Mult[i + 2],
            config["opt_qreg_mult"],
        )
        model_variables.Outlet2DSRS[i + 2] = lo_functions.Outlet2DSRS(
            model_variables.Release_Level[i + 2],
            data.S80_RegRelRates.at[0, "Zone_D1"],
            s80avg_l1,
            data.Pulses.at[
                (
                    model_variables.PlsDay[i + 2] - 1
                    if model_variables.PlsDay[i + 2] - 1 >= 0
                    else len(data.Pulses) - 1
                ),
                f'S-80_L1_{config["schedule"]}',
            ],
            model_variables.Outlet2DS_Mult_2[i + 2],
            data.CE_SLE_turns.at[
                lo_model.at[i + 2, "date"].year - config["start_year"],
                "SLEturn",
            ],
            data.S80_RegRelRates.at[0, "Zone_D2"],
            s80avg_l2,
            data.Pulses.at[
                (
                    model_variables.PlsDay[i + 2] - 1
                    if model_variables.PlsDay[i + 2] - 1 >= 0
                    else len(data.Pulses) - 1
                ),
                f'S-80_L2_{config["schedule"]}',
            ],
            data.S80_RegRelRates.at[0, "Zone_D3"],
            s80avg_l3,
            data.Pulses.at[
                (
                    model_variables.PlsDay[i + 2] - 1
                    if model_variables.PlsDay[i + 2] - 1 >= 0
                    else len(data.Pulses) - 1
                ),
                f'S-80_L3_{config["schedule"]}',
            ],
            data.S80_RegRelRates.at[0, "Zone_C"],
            data.S80_RegRelRates.at[0, "Zone_B"],
            data.S80_RegRelRates.at[0, "Zone_A"],
        )
        diff = model_variables.Outlet2DSRS[i + 2] - lo_model.at[i + 2, "C44RO"]
        model_variables.Outlet2USRG1[i + 2] = max(
            0, diff.values[0] if isinstance(diff, pd.Series) else diff
        )
        model_variables.Sum_Outlet2USRG1[i + 2] = lo_functions.Sum_Outlet2USRG1(
            lo_model.at[i + 2, "date"].day, model_variables.Outlet2USRG1[i + 2]
        )
        model_variables.Outlet2DSBS[i + 2] = lo_functions.Outlet2DSBS(
            model_variables.Release_Level[i + 2],
            model_variables.Sum_Outlet2USRG1[i + 2],
            vlookup1_c[i],
            outlet2_baseflow,
            config["option_s80_baseflow"],
        )
        model_variables.Outlet2USBK[i + 2] = lo_functions.Outlet2USBK(
            model_variables.Lake_Stage[i + 1],
            df_WSMs.at[i + 1, "D1"],
            model_variables.Outlet2USRG[i + 1],
            lo_model.at[i + 2, "C44RO"],
            data.SFWMM_Daily_Outputs.at[i + 2, "S308BK"],
            config["opt_s308"],
            config["s308_bk_const"],
            config["s308_bk_thr"],
        )
        model_variables.ROeast[i + 2] = (
            lo_model.at[i + 2, "C44RO"] - model_variables.Outlet2USBK[i + 2]
        )
        model_variables.Outlet2USBS[i + 2] = lo_functions.Outlet2USBS(
            model_variables.Outlet2DSBS[i + 2],
            model_variables.Outlet2USRG1[i + 2],
            model_variables.ROeast[i + 2],
            config["option_s80_baseflow"],
        )
        model_variables.Sum_Outlet2USBK[i + 2] = lo_functions.Sum_Outlet2USBK(
            lo_model.at[i + 2, "date"].day, model_variables.Outlet2USBK[i + 2]
        )
        model_variables.Outlet2USRG_Code[i + 2] = lo_functions.Outlet2USRG_Code(
            model_variables.Outlet2USRG1[i + 2],
            model_variables.Outlet2USBS[i + 2],
            data.SFWMM_Daily_Outputs.at[i + 2, "S308RG"],
            data.SFWMM_Daily_Outputs.at[i + 2, "STEST"],
            config["option_reg_s77_s308"],
        )
        if config["sim_type"] == 0:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet2USRG(
                model_variables.Outlet2USRG_Code[i + 2],
                data.SFWMM_Daily_Outputs.at[i + 2, "S308RG"],
                data.SFWMM_Daily_Outputs.at[i + 2, "STEST"],
                config["opt_s308"],
                config["s308_rg_const"],
            )
        else:
            if model_variables.Lake_Stage[i + 1] >= 18:
                model_variables.Outlet2USRG[i + 2] = 7200
            elif model_variables.Lake_Stage[i + 1] <= 8:
                model_variables.Outlet2USRG[i + 2] = 0
            elif (tp_lake_s[i] <= p1) and (
                date_range_6[i + 2].month in [1, 2, 3, 4, 11, 12]
            ):
                model_variables.Outlet2USRG[i + 2] = s308_dv[
                    (date_range_6[i + 2].month) - 1
                ]
            elif (tp_lake_s[i] <= p2) and (
                date_range_6[i + 2].month in [5, 6, 7, 8, 9, 10]
            ):
                model_variables.Outlet2USRG[i + 2] = s308_dv[
                    (date_range_6[i + 2].month) - 1
                ]
            else:
                model_variables.Outlet2USRG[i + 2] = 0
        model_variables.Outlet2DS[i + 2] = lo_functions.S80(
            model_variables.ROeast[i + 2],
            model_variables.Outlet2USRG[i + 2],
            data.SFWMM_Daily_Outputs.at[i + 2, "S80"],
            config["s80_const"],
        )
        model_variables.ReLevelCode_3_S77[i + 2] = lo_functions.ReLevelCode_3_S77(
            model_variables.Release_Level[i + 2],
            config["bstar_s77_d1"],
            config["bstar_s77_d2"],
            config["bstar_s77_d3"],
            config["bstar_s77_c"],
            config["bstar_s77_b"],
        )
        model_variables.Outlet1US_Mult[i + 2] = lo_functions.Outlet1US_Mult(
            seasons.at[i, "Season"],
            seasons.at[i, "Month"],
            model_variables.dh_7days[i + 1],
            model_variables.ReLevelCode_1[i + 2],
            model_variables.Fraction_of_Zone_height[i + 1],
            model_variables.ReLevelCode_2[i + 2],
            model_variables.ReLevelCode_3_S77[i + 2],
            config["opt_qreg_mult"],
        )
        model_variables.Outlet1US_Mult_2[i + 2] = lo_functions.Outlet1US_Mult_2(
            lo_model.at[i + 2, "date"].month,
            lo_model.at[i + 2, "date"].day,
            model_variables.PlsDay[i + 2],
            model_variables.Outlet1US_Mult[i + 2 - model_variables.PlsDay[i + 2]],
            model_variables.Outlet1US_Mult[i + 2],
            config["opt_qreg_mult"],
        )
        model_variables.Outlet1USRS[i + 2] = lo_functions.Outlet1USRS(
            model_variables.Release_Level[i + 2],
            data.S77_RegRelRates.at[0, "Zone_D1"],
            s77avg_l1,
            data.Pulses.at[
                (
                    model_variables.PlsDay[i + 2] - 1
                    if model_variables.PlsDay[i + 2] - 1 >= 0
                    else len(data.Pulses) - 1
                ),
                f'S-77_L1_{config["schedule"]}',
            ],
            model_variables.Outlet1US_Mult_2[i + 2],
            lo_model.at[i + 2, "C43RO"],
            data.CE_SLE_turns.at[
                lo_model.at[i + 2, "date"].year - config["start_year"],
                "CEturn",
            ],
            data.S77_RegRelRates.at[0, "Zone_D2"],
            s77avg_l2,
            data.Pulses.at[
                (
                    model_variables.PlsDay[i + 2] - 1
                    if model_variables.PlsDay[i + 2] - 1 >= 0
                    else len(data.Pulses) - 1
                ),
                f'S-77_L2_{config["schedule"]}',
            ],
            model_variables.Zone_Code[i + 1],
            data.S77_RegRelRates.at[0, "Zone_D3"],
            s77avg_l3,
            data.Pulses.at[
                (
                    model_variables.PlsDay[i + 2] - 1
                    if model_variables.PlsDay[i + 2] - 1 >= 0
                    else len(data.Pulses) - 1
                ),
                f'S-77_L3_{config["schedule"]}',
            ],
            data.S77_RegRelRates.at[0, "Zone_C"],
            data.S77_RegRelRates.at[0, "Zone_B"],
            data.S77_RegRelRates.at[0, "Zone_A"],
            config["opt_outlet1_dsrg"],
        )
        model_variables.Sum_Outlet1USRS[i + 2] = lo_functions.Sum_Outlet1USRS(
            lo_model.at[i + 2, "date"].day, model_variables.Outlet1USRS[i + 2]
        )
        model_variables.Outlet1USBK[i + 2] = lo_functions.Outlet1USBK(
            model_variables.Lake_Stage[i + 1],
            model_variables.Outlet1USRS[i + 2],
            model_variables.Outlet1USBSAP[i + 1],
            model_variables.Outlet1USEWS[i + 1],
            lo_model.at[i + 2, "C43RO"],
            data.SFWMM_Daily_Outputs.at[i + 2, "S77BK"],
            config["outlet1_usbk_switch"],
            config["outlet1_usbk_threshold"],
        )
        model_variables.ROwest[i + 2] = (
            lo_model.at[i + 2, "C43RO"] - model_variables.Outlet1USBK[i + 2]
        )
        model_variables.Outlet1DSBS[i + 2] = lo_functions.Outlet1DSBS(
            model_variables.Release_Level[i + 2],
            model_variables.Sum_Outlet1USRS[i + 2],
            vlookup2_c[i],
            outlet1_baseflow,
            config["option_s77_baseflow"],
        )
        model_variables.Outlet1USBS[i + 2] = lo_functions.Outlet1USBS(
            model_variables.Outlet1DSBS[i + 2],
            model_variables.Outlet1USRS[i + 2],
            model_variables.ROwest[i + 2],
            config["option_s77_baseflow"],
        )
        # Define THC Class Normal or above
        if i < (n_rows - 2):
            model_variables.Post_Ap_Baseflow[i] = thc_class.THC_Class(
                config,
                i,
                model_variables.THC_Class_normal_or_above,
                model_variables.Lake_O_Stage_AP,
                model_variables.Lake_O_Schedule_Zone,
                model_variables.LStgCorres,
                model_variables.LowChance_Check,
                model_variables.Outlet1USRS_AP,
                model_variables.Outlet1USBS_AP,
                model_variables.Outlet1USRS_Pre_AP_S77_Baseflow,
                model_variables.Forecast_D_Sal,
                model_variables.n30d_mavg,
                model_variables.n30davgForecast,
                model_variables.LORS08_bf_rel,
                model_variables.LDS_LC6_1,
                model_variables.S_O,
                model_variables.All_4,
                model_variables.Sabf,
                model_variables.Swbf,
                model_variables.Swbu,
                model_variables.All_4andStage,
                model_variables.All_4andStagein,
                model_variables.P_AP_BF_Stg,
                model_variables.Logic_test_1,
                model_variables.Post_Ap_Baseflow,
                model_variables.Outlet1USRSplusPreAPS77bsf,
                model_variables.AndEstNeedsLakeWater,
                model_variables.AndLowChance61stagelessth11,
                model_variables.ATHCnora,
                model_variables.Choose_PAPEWS_1,
                model_variables.Choose_PAPEWS_2,
                model_variables.Post_AP_EWS,
                model_variables.Post_AP_Baseflow_EWS_cfs,
                adaptive_protocol_df,
                model_variables.Lake_Stage,
                model_variables.Zone_Code,
                df_WSMs,
                targ_stg_df,
                model_variables.Outlet1USRS,
                model_variables.Outlet1USBS,
                data.Estuary_needs_water,
                choose_1,
                model_variables.WSM_Zone,
            )["Post_Ap_Baseflow"]
            model_variables.Post_AP_EWS[i] = thc_class.THC_Class(
                config,
                i,
                model_variables.THC_Class_normal_or_above,
                model_variables.Lake_O_Stage_AP,
                model_variables.Lake_O_Schedule_Zone,
                model_variables.LStgCorres,
                model_variables.LowChance_Check,
                model_variables.Outlet1USRS_AP,
                model_variables.Outlet1USBS_AP,
                model_variables.Outlet1USRS_Pre_AP_S77_Baseflow,
                model_variables.Forecast_D_Sal,
                model_variables.n30d_mavg,
                model_variables.n30davgForecast,
                model_variables.LORS08_bf_rel,
                model_variables.LDS_LC6_1,
                model_variables.S_O,
                model_variables.All_4,
                model_variables.Sabf,
                model_variables.Swbf,
                model_variables.Swbu,
                model_variables.All_4andStage,
                model_variables.All_4andStagein,
                model_variables.P_AP_BF_Stg,
                model_variables.Logic_test_1,
                model_variables.Post_Ap_Baseflow,
                model_variables.Outlet1USRSplusPreAPS77bsf,
                model_variables.AndEstNeedsLakeWater,
                model_variables.AndLowChance61stagelessth11,
                model_variables.ATHCnora,
                model_variables.Choose_PAPEWS_1,
                model_variables.Choose_PAPEWS_2,
                model_variables.Post_AP_EWS,
                model_variables.Post_AP_Baseflow_EWS_cfs,
                adaptive_protocol_df,
                model_variables.Lake_Stage,
                model_variables.Zone_Code,
                df_WSMs,
                targ_stg_df,
                model_variables.Outlet1USRS,
                model_variables.Outlet1USBS,
                data.Estuary_needs_water,
                choose_1,
                model_variables.WSM_Zone,
            )["Post_AP_EWS"]
        model_variables.Outlet1USBSAP[i + 2] = lo_functions.Outlet1USBSAP(
            model_variables.Outlet1USBS[i + 2],
            model_variables.Post_Ap_Baseflow[i],
            config["opt_adap_prot"],
        )
        model_variables.Outlet1USEWS[i + 2] = lo_functions.Outlet1USEWS(
            model_variables.Post_AP_EWS[i],
            data.SFWMM_Daily_Outputs.at[i + 2, "CAEST"],
            config["outlet1_usews_switch"],
            config["opt_adap_prot"],
        )
        if config["sim_type"] == 0:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet1USREG(
                model_variables.Outlet1USRS[i + 2],
                model_variables.Outlet1USBSAP[i + 2],
                data.SFWMM_Daily_Outputs.at[i + 2, "S77RG"],
                config["outlet1_usreg_switch"],
                config["option_reg_s77_s308"],
            )
        else:
            if model_variables.Lake_Stage[i + 1] >= 18:
                model_variables.Outlet1USREG[i + 2] = 7800
            elif model_variables.Lake_Stage[i + 1] <= 8:
                model_variables.Outlet1USREG[i + 2] = 0
            elif (tp_lake_s[i] <= p1) and (
                date_range_6[i + 2].month in [1, 2, 3, 4, 11, 12]
            ):
                model_variables.Outlet1USREG[i + 2] = s77_dv[
                    (date_range_6[i + 2].month) - 1
                ]
            elif (tp_lake_s[i] <= p2) and (
                date_range_6[i + 2].month in [5, 6, 7, 8, 9, 10]
            ):
                model_variables.Outlet1USREG[i + 2] = s77_dv[
                    (date_range_6[i + 2].month) - 1
                ]
            else:
                model_variables.Outlet1USREG[i + 2] = 0
        model_variables.Outlet1DS[i + 2] = lo_functions.Outlet1DS(
            model_variables.Outlet1USREG[i + 2],
            model_variables.Outlet1USEWS[i + 2],
            model_variables.ROwest[i + 2],
            data.SFWMM_Daily_Outputs.at[i + 2, "S79"],
            config["outlet1_ds_switch"],
        )
        model_variables.TotRegEW[i + 2] = (
            model_variables.Outlet1USREG[i + 2] + model_variables.Outlet2USRG[i + 2]
        ) * 1.9835
        model_variables.Choose_WCA[i + 2] = lo_functions.Choose_WCA(
            data.SFWMM_Daily_Outputs.at[i + 2, "RegWCA"],
            config["option_reg_wca"],
            config["constant_reg_wca"],
        )
        model_variables.RegWCA[i + 2] = min(
            config["max_cap_reg_wca"],
            config["multiplier_reg_wca"] * model_variables.Choose_WCA[i + 2],
        )
        model_variables.Choose_L8C51[i + 2] = lo_functions.Choose_L8C51(
            data.SFWMM_Daily_Outputs.at[i + 2, "RegL8C51"],
            config["option_reg_l8_c51"],
            config["constant_reg_l8_c51"],
        )
        model_variables.RegL8C51[i + 2] = min(
            config["max_cap_reg_l8_c51"],
            config["multiplier_reg_l8_c51"] * model_variables.Choose_L8C51[i + 2],
        )
        model_variables.TotRegSo[i + 2] = (
            model_variables.RegWCA[i + 2] + model_variables.RegL8C51[i + 2]
        ) * 1.9835
        model_variables.Stage2ar[i + 2] = stg_sto_ar.stg2ar(model_variables.Lake_Stage[i + 1], 0)
        model_variables.Stage2marsh[i + 2] = stg_sto_ar.stg2mar(
            model_variables.Lake_Stage[i + 1], 0
        )
        model_variables.RF[i + 2] = data.RF_Vol.at[i + 2, "RFVol_acft"]
        model_variables.ET[i + 2] = lo_functions.ET(
            data.SFWMM_Daily_Outputs.at[i + 2, "et_dry"],
            model_variables.Stage2ar[i + 2],
            data.SFWMM_Daily_Outputs.at[i + 2, "et_litoral"],
            model_variables.Stage2marsh[i + 2],
            data.SFWMM_Daily_Outputs.at[i + 2, "et_open"],
            data.ET_Vol.at[i + 2, "ETVol_acft"],
            config["et_switch"],
        )
        model_variables.Choose_WSA_1[i + 2] = lo_functions.Choose_WSA_1(
            df_WSMs.at[i + 2, "WSM1"],
            config["opt_wsa"],
            config["wsa_trig2"],
            config["wsa_off2"],
        )
        model_variables.Choose_WSA_2[i + 2] = lo_functions.Choose_WSA_2(
            df_WSMs.at[i + 2, "WSM1"],
            config["opt_wsa"],
            config["wsa_trig1"],
            config["wsa_off1"],
        )
        model_variables.WSA_MIA[i + 2] = lo_functions.WSA_MIA(
            wca_stages_df.at[i, "Are WCA stages too low?"],
            tc_lonino_df.at[i, "LONINO_Seasonal_Classes"],
            model_variables.Lake_Stage[i + 1],
            model_variables.Choose_WSA_1[i + 2],
            data.EAA_MIA_RUNOFF.at[i, "MIA"],
            data.EAA_MIA_RUNOFF.at[i, "S3PMP"],
            model_variables.Choose_WSA_2[i + 2],
            config["opt_wsa"],
            config["wsa_thc"],
            config["mia_cap2"],
            config["mia_cap1"],
        )
        model_variables.WSA_NNR[i + 2] = lo_functions.WSA_NNR(
            wca_stages_df.at[i, "Are WCA stages too low?"],
            tc_lonino_df.at[i, "LONINO_Seasonal_Classes"],
            model_variables.Lake_Stage[i + 1],
            model_variables.Choose_WSA_1[i + 2],
            data.EAA_MIA_RUNOFF.at[i, "NNR"],
            data.EAA_MIA_RUNOFF.at[i, "S2PMP"],
            model_variables.Choose_WSA_2[i + 2],
            config["opt_wsa"],
            config["wsa_thc"],
            config["nnr_cap2"],
            config["nnr_cap1"],
        )
        model_variables.DSto[i + 2] = (
            model_variables.NI_Supply[i + 2]
            + model_variables.RF[i + 2]
            - model_variables.ET[i + 2]
            + 1.9835
            * (
                model_variables.Outlet2USBK[i + 2]
                + model_variables.Outlet1USBK[i + 2]
                + model_variables.WSA_MIA[i + 2]
                + model_variables.WSA_NNR[i + 2]
                - model_variables.Outlet1USEWS[i + 2]
            )
            - model_variables.TotRegEW[i + 2]
            - model_variables.TotRegSo[i + 2]
            + storage_deviation[i + 2]
        )
        model_variables.Storage[i + 2] = lo_functions.Storage(
            model_variables.DayFlags[i + 2],
            model_variables.Storage[i],
            start_storage,
            model_variables.Storage[i + 1],
            model_variables.DSto[i + 2],
        )
        model_variables.Lake_Stage[i + 2] = lo_functions.Lake_Stage(
            stg_sto_ar.stg2sto(model_variables.Storage[i + 2], 1),
            data.SFWMM_Daily_Outputs.at[i + 2, "EOD Stg(ft,NGVD)"],
            config["option_stage"],
        )
        # if M_var.Lake_Stage[i+2] >= 18:
        #     Flood[i+2] = 1
        # else:
        #     Flood[i+2] = 0
        lo_model["Cutback"] = model_variables.Cut_back
        lo_model["Stage"] = model_variables.Lake_Stage
        lo_model["Storage"] = model_variables.Storage
        lo_model["S77_Q"] = model_variables.Outlet1USREG
        lo_model["S308_Q"] = model_variables.Outlet2USRG
        lo_model["S77EW"] = model_variables.Outlet1USEWS
        lo_model["TotRegEW"] = model_variables.TotRegEW
        lo_model["TotRegSo"] = model_variables.TotRegSo

    # Write out the results to a file - Needed because
    # execute_loone.py calls this script as a subprocess and can't get
    # this data.
    lo_model.to_csv("LOONE_Q_Outputs.csv")

    return [lo_model]


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
