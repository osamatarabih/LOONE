import os
import argparse
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from calendar import monthrange
from typing import List
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


def _calculate_wsm_zone(
    i: int, model_variables: object, lo_functions: object, df_WSMs: pd.DataFrame
) -> None:
    """
    Calculate the WSM zone for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.

    Returns:
        None
    """
    model_variables.WSM_Zone[i + 2] = lo_functions.WSM_Zone(
        model_variables.Lake_Stage[i + 1],
        df_WSMs.at[i + 1, "WSM4"],
        df_WSMs.at[i + 1, "WSM3"],
        df_WSMs.at[i + 1, "WSM2"],
        df_WSMs.at[i + 1, "WSM1"],
    )


def _calculate_max_supply(
    i: int,
    model_variables: object,
    lo_functions: object,
    water_demand: pd.DataFrame,
    config: dict,
) -> None:
    """
    Calculate the maximum supply for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        water_demand (pd.DataFrame): The water demand DataFrame.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Max_Supply[i + 2] = lo_functions.Max_Supply(
        model_variables.WSM_Zone[i + 2],
        water_demand.at[i, "Daily_demand"],
        config["z1_cutback"],
        config["z2_cutback"],
        config["z3_cutback"],
        config["z4_cutback"],
    )


def _calculate_losa_supply(
    i: int,
    model_variables: object,
    lo_functions: object,
    lo_model: pd.DataFrame,
    config: dict,
) -> None:
    """
    Calculate the LOSA supply for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        lo_model (pd.DataFrame): The model DataFrame.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.LOSA_Supply[i + 2] = lo_functions.LOSA_Supply(
        model_variables.WSM_Zone[i + 2],
        lo_model.at[i + 2, "LOSA_dmd_SFWMM"],
        model_variables.Max_Supply[i + 2],
        config["opt_losa_ws"],
    )


def _calculate_ni_supply(
    i: int, model_variables: object, lo_model: pd.DataFrame
) -> None:
    """
    Calculate the NI supply for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_model (pd.DataFrame): The model DataFrame.

    Returns:
        None
    """
    model_variables.NI_Supply[i + 2] = (
        lo_model.at[i + 2, "Net_Inflow"] - model_variables.LOSA_Supply[i + 2]
    )


def _calculate_cutback(i: int, model_variables: object, lo_model: pd.DataFrame) -> None:
    """
    Calculate the cutback for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_model (pd.DataFrame): The model DataFrame.

    Returns:
        None
    """
    cutback = lo_model.at[i + 2, "LOSA_dmd_SFWMM"] - model_variables.LOSA_Supply[i + 2]
    model_variables.Cut_back[i + 2] = cutback


def _calculate_demand_not_supplied(
    i: int, model_variables: object, lo_model: pd.DataFrame
) -> None:
    """
    Calculate the demand not supplied for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_model (pd.DataFrame): The model DataFrame.

    Returns:
        None
    """
    if lo_model.at[i + 2, "LOSA_dmd_SFWMM"] == 0:
        demand_not_supplied = 0
    else:
        demand_not_supplied = (
            model_variables.Cut_back[i + 2] / lo_model.at[i + 2, "LOSA_dmd_SFWMM"]
        ) * 100
    model_variables.Dem_N_Sup[i + 2] = demand_not_supplied


def _calculate_initial_zone_code_and_lo_zone(
    model_variables: object, lo_functions: object, df_WSMs: pd.DataFrame
) -> None:
    """
    Calculate the Zone_Code and LO_Zone for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.

    Returns:
        None
    """
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


def _calculate_zone_code(
    i: int, model_variables: object, lo_functions: object, df_WSMs: pd.DataFrame
) -> None:
    """
    Calculate the Zone_Code for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.

    Returns:
        None
    """
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


def _generate_zone_column(
    i: int, model_variables: object, lo_functions: object
) -> None:
    """
    Generate the Zone columns for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.

    Returns:
        None
    """
    model_variables.LO_Zone[i + 1] = lo_functions.LO_Zone(
        model_variables.Zone_Code[i + 1]
    )


def _calculate_decision_tree_variables(
    i: int,
    model_variables: object,
    dec_tree_functions: object,
    tc_lonino_df: pd.DataFrame,
    config: dict,
    seasons: pd.DataFrame,
    dec_tree_df: pd.DataFrame,
    df_WSMs: pd.DataFrame,
) -> None:
    """
    Calculate the decision tree variables for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        dec_tree_functions (object): The decision tree functions object.
        tc_lonino_df (pd.DataFrame): The TC LONINO DataFrame.
        config (dict): The configuration dictionary.
        seasons (pd.DataFrame): The seasons DataFrame.
        dec_tree_df (pd.DataFrame): The decision tree DataFrame.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.

    Returns:
        None
    """
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


def _calculate_release_levels(
    i: int,
    startdate: datetime,
    lo_model: pd.DataFrame,
    model_variables: object,
    lo_functions: object,
    config: dict,
    tc_lonino_df: pd.DataFrame,
) -> None:
    """
    Calculate the release levels for the given index.

    Args:
        i (int): The current index.
        startdate (datetime): The start date.
        lo_model (pd.DataFrame): The model DataFrame.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        config (dict): The configuration dictionary.
        tc_lonino_df (pd.DataFrame): The TC LONINO DataFrame.

    Returns:
        None
    """
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
                config["cs_flag"] == 0 or startdate.year == lo_model.at[i, "date"].year
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


def _calculate_zone_codes(
    i: int, model_variables: object, lo_functions: object, df_WSMs: pd.DataFrame
) -> None:
    """
    Calculate the zone codes for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.

    Returns:
        None
    """
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


def _calculate_release_level_codes(
    i: int, model_variables: object, lo_functions: object, config: dict
) -> None:
    """
    Calculate the release level codes for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
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


def _calculate_relevel_code_3_s77(
    i: int, model_variables: object, lo_functions: object, config: dict
) -> None:
    """
    Calculate the ReLevelCode_3_S77 for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.ReLevelCode_3_S77[i + 2] = lo_functions.ReLevelCode_3_S77(
        model_variables.Release_Level[i + 2],
        config["bstar_s77_d1"],
        config["bstar_s77_d2"],
        config["bstar_s77_d3"],
        config["bstar_s77_c"],
        config["bstar_s77_b"],
    )


def _calculate_outlet2ds_mult(
    i: int,
    lo_model: pd.DataFrame,
    model_variables: object,
    lo_functions: object,
    seasons: pd.DataFrame,
    config: dict,
) -> None:
    """
    Calculate the Outlet2DS multiplier for the given index.

    Args:
        i (int): The current index.
        lo_model (pd.DataFrame): The model DataFrame.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        seasons (pd.DataFrame): The seasons DataFrame.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
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


def _calculate_outlet2dsrs(
    i: int,
    lo_model: pd.DataFrame,
    model_variables: object,
    lo_functions: object,
    data: object,
    s80avg_l1: float,
    s80avg_l2: float,
    s80avg_l3: float,
    config: dict,
) -> None:
    """
    Calculate the Outlet2DSRS for the given index.

    Args:
        i (int): The current index.
        lo_model (pd.DataFrame): The model DataFrame.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        s80avg_l1 (float): The S80 average level 1.
        s80avg_l2 (float): The S80 average level 2.
        s80avg_l3 (float): The S80 average level 3.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    year_index = lo_model.at[i + 2, "date"].year - config["start_year"]

    if year_index in data.CE_SLE_turns.index:
        sle_turn = data.CE_SLE_turns.at[year_index, "SLEturn"]
    else:
        #TODO not sure that this should be hard coded like this, but this is because of the forecast things
        sle_turn = 1
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
            f'S-80_L1',
        ],
        model_variables.Outlet2DS_Mult_2[i + 2],
        sle_turn,
        data.S80_RegRelRates.at[0, "Zone_D2"],
        s80avg_l2,
        data.Pulses.at[
            (
                model_variables.PlsDay[i + 2] - 1
                if model_variables.PlsDay[i + 2] - 1 >= 0
                else len(data.Pulses) - 1
            ),
            f'S-80_L2',
        ],
        data.S80_RegRelRates.at[0, "Zone_D3"],
        s80avg_l3,
        data.Pulses.at[
            (
                model_variables.PlsDay[i + 2] - 1
                if model_variables.PlsDay[i + 2] - 1 >= 0
                else len(data.Pulses) - 1
            ),
            f'S-80_L3',
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


def _calculate_outlet2dsbs(
    i: int,
    model_variables: object,
    lo_functions: object,
    vlookup1_c: list,
    outlet2_baseflow: float,
    config: dict,
) -> None:
    """
    Calculate the Outlet2DSBS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        vlookup1_c (list): The vlookup1_c list.
        outlet2_baseflow (float): The outlet2 baseflow.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet2DSBS[i + 2] = lo_functions.Outlet2DSBS(
        model_variables.Release_Level[i + 2],
        model_variables.Sum_Outlet2USRG1[i + 2],
        vlookup1_c[i],
        outlet2_baseflow,
        config["option_s80_baseflow"],
    )


def _calculate_outlet2usbk(
    i: int,
    model_variables: object,
    lo_functions: object,
    df_WSMs: pd.DataFrame,
    lo_model: pd.DataFrame,
    data: object,
    config: dict,
) -> None:
    """
    Calculate the Outlet2USBK for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.
        lo_model (pd.DataFrame): The model DataFrame.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
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


def _calculate_outlet2usrg_code(
    i: int,
    model_variables: object,
    lo_functions: object,
    data: object,
    config: dict,
    p1: float,
    p2: float,
    s308_dv: float,
    tp_lake_s: float,
    lo_model: pd.DataFrame,
    date_range_one_day: pd.date_range,
) -> None:
    """
    Calculate the Outlet2USRG code for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.
        p1 (float): The P1 value.
        p2 (float): The P2 value.
        s308_dv (float): The S308 DV value.
        tp_lake_s (float): The TP lake S value.
        date_range_one_day (pd.date_range): The date range.

    Returns:
        None
    """
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
    elif config["sim_type"] == 1:
        if model_variables.Lake_Stage[i + 1] >= 18:
            model_variables.Outlet2USRG[i + 2] = 7200
        elif model_variables.Lake_Stage[i + 1] <= 8:
            model_variables.Outlet2USRG[i + 2] = 0
        elif (tp_lake_s[i] <= p1) and (
            date_range_one_day[i + 2].month in [1, 2, 3, 4, 11, 12]
        ):
            model_variables.Outlet2USRG[i + 2] = s308_dv[
                (date_range_one_day[i + 2].month) - 1
            ]
        elif (tp_lake_s[i] <= p2) and (
            date_range_one_day[i + 2].month in [5, 6, 7, 8, 9, 10]
        ):
            model_variables.Outlet2USRG[i + 2] = s308_dv[
                (date_range_one_day[i + 2].month) - 1
            ]
        else:
            model_variables.Outlet2USRG[i + 2] = 0
    elif config["sim_type"] == 3:
        if model_variables.Lake_Stage[i + 1] >= 18:
            model_variables.Outlet2USRG[i + 2] = 7200
        elif model_variables.Lake_Stage[i + 1] <= 8:
            model_variables.Outlet2USRG[i + 2] = 0
        elif lo_model.at[i + 2, "date"].month == 1:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[0]
            )
        elif lo_model.at[i + 2, "date"].month == 2:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[1]
            )
        elif lo_model.at[i + 2, "date"].month == 3:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[2]
            )
        elif lo_model.at[i + 2, "date"].month == 4:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[3]
            )
        elif lo_model.at[i + 2, "date"].month == 5:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[4]
            )
        elif lo_model.at[i + 2, "date"].month == 6:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[5]
            )

        elif lo_model.at[i + 2, "date"].month == 7:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[6]
            )
        elif lo_model.at[i + 2, "date"].month == 8:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[7]
            )
        elif lo_model.at[i + 2, "date"].month == 9:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[8]
            )
        elif lo_model.at[i + 2, "date"].month == 10:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[9]
            )
        elif lo_model.at[i + 2, "date"].month == 11:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[10]
            )
        elif lo_model.at[i + 2, "date"].month == 12:
            model_variables.Outlet2USRG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[11]
            )

def _calculate_outlet2ds(
    i: int, model_variables: object, lo_functions: object, data: object, config: dict
) -> None:
    """
    Calculate the Outlet2DS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet2DS[i + 2] = lo_functions.S80(
        model_variables.ROeast[i + 2],
        model_variables.Outlet2USRG[i + 2],
        data.SFWMM_Daily_Outputs.at[i + 2, "S80"],
        config["s80_const"],
    )


def _calculate_outlet1us_mult(
    i: int,
    model_variables: object,
    lo_functions: object,
    seasons: pd.DataFrame,
    config: dict,
) -> None:
    """
    Calculate the Outlet1US_Mult for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        seasons (pd.DataFrame): The seasons DataFrame.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
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


def _calculate_outlet1us_mult_2(
    i: int,
    model_variables: object,
    lo_functions: object,
    lo_model: pd.DataFrame,
    config: dict,
) -> None:
    """
    Calculate the Outlet1US_Mult_2 for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        lo_model (pd.DataFrame): The model DataFrame.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet1US_Mult_2[i + 2] = lo_functions.Outlet1US_Mult_2(
        lo_model.at[i + 2, "date"].month,
        lo_model.at[i + 2, "date"].day,
        model_variables.PlsDay[i + 2],
        model_variables.Outlet1US_Mult[i + 2 - model_variables.PlsDay[i + 2]],
        model_variables.Outlet1US_Mult[i + 2],
        config["opt_qreg_mult"],
    )


def _calculate_outlet1usrs(
    i: int,
    model_variables: object,
    lo_functions: object,
    lo_model: pd.DataFrame,
    data: object,
    config: dict,
    s77avg_l1: float,
    s77avg_l2: float,
    s77avg_l3: float,
) -> None:
    """
    Calculate the Outlet1USRS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        lo_model (pd.DataFrame): The model DataFrame.
        data (object): The data object.
        config (dict): The configuration dictionary.
        s77avg_l1 (float): The S77 average level 1.
        s77avg_l2 (float): The S77 average level 2.
        s77avg_l3 (float): The S77 average level 3.

    Returns:
        None
    """

    target_year = lo_model.at[i + 2, "date"].year

    available_years = data.CE_SLE_turns["Year"].values
    closest_year = available_years[np.argmin(np.abs(available_years - target_year))]

    CEturn_value = data.CE_SLE_turns.loc[data.CE_SLE_turns["Year"] == closest_year, "CEturn"].values[0]

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
            f'S-77_L1',
        ],
        model_variables.Outlet1US_Mult_2[i + 2],
        lo_model.at[i + 2, "C43RO"],
        CEturn_value,
        data.S77_RegRelRates.at[0, "Zone_D2"],
        s77avg_l2,
        data.Pulses.at[
            (
                model_variables.PlsDay[i + 2] - 1
                if model_variables.PlsDay[i + 2] - 1 >= 0
                else len(data.Pulses) - 1
            ),
            f'S-77_L2',
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
            f'S-77_L3',
        ],
        data.S77_RegRelRates.at[0, "Zone_C"],
        data.S77_RegRelRates.at[0, "Zone_B"],
        data.S77_RegRelRates.at[0, "Zone_A"],
        config["opt_outlet1_dsrg"],
)


def _calculate_sum_outlet1usrs(
    i: int, model_variables: object, lo_functions: object, lo_model: pd.DataFrame
) -> None:
    """
    Calculate the Sum_Outlet1USRS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        lo_model (pd.DataFrame): The model DataFrame.

    Returns:
        None
    """
    model_variables.Sum_Outlet1USRS[i + 2] = lo_functions.Sum_Outlet1USRS(
        lo_model.at[i + 2, "date"].day, model_variables.Outlet1USRS[i + 2]
    )


def _calculate_outlet1usbk(
    i: int,
    model_variables: object,
    lo_functions: object,
    lo_model: pd.DataFrame,
    data: object,
    config: dict,
) -> None:
    """
    Calculate the Outlet1USBK for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        lo_model (pd.DataFrame): The model DataFrame.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
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


def _calculate_rowest(i: int, model_variables: object, lo_model: pd.DataFrame) -> None:
    """
    Calculate the ROwest for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_model (pd.DataFrame): The model DataFrame.

    Returns:
        None
    """
    model_variables.ROwest[i + 2] = (
        lo_model.at[i + 2, "C43RO"] - model_variables.Outlet1USBK[i + 2]
    )


def _calculate_outlet1dsbs(
    i: int,
    model_variables: object,
    lo_functions: object,
    vlookup2_c: list,
    outlet1_baseflow: float,
    config: dict,
) -> None:
    """
    Calculate the Outlet1DSBS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        vlookup2_c (list): The vlookup2_c list.
        outlet1_baseflow (float): The outlet1 baseflow.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet1DSBS[i + 2] = lo_functions.Outlet1DSBS(
        model_variables.Release_Level[i + 2],
        model_variables.Sum_Outlet1USRS[i + 2],
        vlookup2_c[i],
        outlet1_baseflow,
        config["option_s77_baseflow"],
    )


def _calculate_outlet1usbs(
    i: int, model_variables: object, lo_functions: object, config: dict
) -> None:
    """
    Calculate the Outlet1USBS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet1USBS[i + 2] = lo_functions.Outlet1USBS(
        model_variables.Outlet1DSBS[i + 2],
        model_variables.Outlet1USRS[i + 2],
        model_variables.ROwest[i + 2],
        config["option_s77_baseflow"],
    )


def _calculate_outlet1usbsap(
    i: int, model_variables: object, lo_functions: object, config: dict
) -> None:
    """
    Calculate the Outlet1USBSAP for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet1USBSAP[i + 2] = lo_functions.Outlet1USBSAP(
        model_variables.Outlet1USBS[i + 2],
        model_variables.Post_Ap_Baseflow[i],
        config["opt_adap_prot"],
    )


def _calculate_outlet1usews(
    i: int, model_variables: object, lo_functions: object, data: object, config: dict
) -> None:
    """
    Calculate the Outlet1USEWS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet1USEWS[i + 2] = lo_functions.Outlet1USEWS(
        model_variables.Post_AP_EWS[i],
        data.SFWMM_Daily_Outputs.at[i + 2, "CAEST"],
        config["outlet1_usews_switch"],
        config["opt_adap_prot"],
    )


def _define_thc_class_normal_or_above(
    i: int,
    n_rows: int,
    model_variables: object,
    config: dict,
    adaptive_protocol_df: pd.DataFrame,
    df_WSMs: pd.DataFrame,
    targ_stg_df: pd.DataFrame,
    data: object,
    choose_1: float,
) -> None:
    """
    Define THC Class Normal or above for the given index.

    Args:
        i (int): The current index.
        n_rows (int): The number of rows.
        model_variables (object): The model variables object.
        config (dict): The configuration dictionary.
        adaptive_protocol_df (pd.DataFrame): The adaptive protocol DataFrame.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.
        targ_stg_df (pd.DataFrame): The target stage DataFrame.
        data (object): The data object.
        choose_1 (float): The choose_1 value.

    Returns:
        None
    """
    if i < (n_rows - 2):
        thc_class_result = thc_class.THC_Class(
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
        )
        model_variables.Post_Ap_Baseflow[i] = thc_class_result["Post_Ap_Baseflow"]
        model_variables.Post_AP_EWS[i] = thc_class_result["Post_AP_EWS"]


def _calculate_outlet1usreg(
    i: int,
    model_variables: object,
    lo_functions: object,
    data: object,
    config: dict,
    p1: float,
    p2: float,
    s77_dv: float,
    tp_lake_s: float,
    date_range_one_day: pd.DatetimeIndex,
    lo_model: pd.DataFrame,
    sensitivity_analysis_params,
) -> None:
    """
    Calculate the Outlet1USREG for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.
        p1 (float): The p1 value.
        p2 (float): The p2 value.
        s77_dv (float): The s77_dv value.
        tp_lake_s (float): The tp_lake_s value.
        date_range_one_day (pd.DatetimeIndex): The date range.
        lo_model (pd.DataFrame): The model DataFrame.
        sensitivity_analysis_params (list): Holds a release rate for each month of the year.

    Returns:
        None
    """
    if config["sim_type"] == 0:
        model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet1USREG(
            model_variables.Outlet1USRS[i + 2],
            model_variables.Outlet1USBSAP[i + 2],
            data.SFWMM_Daily_Outputs.at[i + 2, "S77RG"],
            config["outlet1_usreg_switch"],
            config["option_reg_s77_s308"],
        )
    elif config["sim_type"] == 1:
        if model_variables.Lake_Stage[i + 1] >= 18:
            model_variables.Outlet1USREG[i + 2] = 7800
        elif model_variables.Lake_Stage[i + 1] <= 8:
            model_variables.Outlet1USREG[i + 2] = 0
        elif (tp_lake_s[i] <= p1) and (
            date_range_one_day[i + 2].month in [1, 2, 3, 4, 11, 12]
        ):
            model_variables.Outlet1USREG[i + 2] = s77_dv[
                (date_range_one_day[i + 2].month) - 1
            ]
        elif (tp_lake_s[i] <= p2) and (
            date_range_one_day[i + 2].month in [5, 6, 7, 8, 9, 10]
        ):
            model_variables.Outlet1USREG[i + 2] = s77_dv[
                (date_range_one_day[i + 2].month) - 1
            ]
        else:
            model_variables.Outlet1USREG[i + 2] = 0
    elif config["sim_type"] == 3:
        if model_variables.Lake_Stage[i + 1] >= 18:
            model_variables.Outlet1USREG[i + 2] = 7800
        elif model_variables.Lake_Stage[i + 1] <= 8:
            model_variables.Outlet1USREG[i + 2] = 0
        elif lo_model.at[i + 2, "date"].month == 1:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[12]
            )
        elif lo_model.at[i + 2, "date"].month == 2:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[13]
            )
        elif lo_model.at[i + 2, "date"].month == 3:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[14]
            )
        elif lo_model.at[i + 2, "date"].month == 4:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[15]
            )
        elif lo_model.at[i + 2, "date"].month == 5:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[16]
            )
        elif lo_model.at[i + 2, "date"].month == 6:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[17]
            )
        elif lo_model.at[i + 2, "date"].month == 7:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[18]
            )
        elif lo_model.at[i + 2, "date"].month == 8:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[19]
            )
        elif lo_model.at[i + 2, "date"].month == 9:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[20]
            )
        elif lo_model.at[i + 2, "date"].month == 10:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[21]
            )
        elif lo_model.at[i + 2, "date"].month == 11:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[22]
            )
        elif lo_model.at[i + 2, "date"].month == 12:
            model_variables.Outlet1USREG[i + 2] = lo_functions.Outlet_Rel_Sim(
                model_variables.Release_Level[i + 2], sensitivity_analysis_params[23]
            )


def _calculate_outlet1ds(
    i: int, model_variables: object, lo_functions: object, data: object, config: dict
) -> None:
    """
    Calculate the Outlet1DS for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Outlet1DS[i + 2] = lo_functions.Outlet1DS(
        model_variables.Outlet1USREG[i + 2],
        model_variables.Outlet1USEWS[i + 2],
        model_variables.ROwest[i + 2],
        data.SFWMM_Daily_Outputs.at[i + 2, "S79"],
        config["outlet1_ds_switch"],
    )


def _calculate_tot_reg_ew(i: int, model_variables: object) -> None:
    """
    Calculate the TotRegEW for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.

    Returns:
        None
    """
    model_variables.TotRegEW[i + 2] = (
        model_variables.Outlet1USREG[i + 2] + model_variables.Outlet2USRG[i + 2]
    ) * 1.9835


def _calculate_choose_wca(
    i: int, model_variables: object, lo_functions: object, data: object, config: dict
) -> None:
    """
    Calculate the Choose_WCA and RegWCA for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Choose_WCA[i + 2] = lo_functions.Choose_WCA(
        data.SFWMM_Daily_Outputs.at[i + 2, "RegWCA"],
        config["option_reg_wca"],
        config["constant_reg_wca"],
    )
    model_variables.RegWCA[i + 2] = min(
        config["max_cap_reg_wca"],
        config["multiplier_reg_wca"] * model_variables.Choose_WCA[i + 2],
    )


def _calculate_choose_l8c51(
    i: int, model_variables: object, lo_functions: object, data: object, config: dict
) -> None:
    """
    Calculate the Choose_L8C51 and RegL8C51 for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Choose_L8C51[i + 2] = lo_functions.Choose_L8C51(
        data.SFWMM_Daily_Outputs.at[i + 2, "RegL8C51"],
        config["option_reg_l8_c51"],
        config["constant_reg_l8_c51"],
    )
    model_variables.RegL8C51[i + 2] = min(
        config["max_cap_reg_l8_c51"],
        config["multiplier_reg_l8_c51"] * model_variables.Choose_L8C51[i + 2],
    )


def _calculate_tot_reg_so(i: int, model_variables: object) -> None:
    """
    Calculate the TotRegSo for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.

    Returns:
        None
    """
    model_variables.TotRegSo[i + 2] = (
        model_variables.RegWCA[i + 2] + model_variables.RegL8C51[i + 2]
    ) * 1.9835


def _calculate_stage2ar_and_stage2marsh(
    i: int, model_variables: object, stg_sto_ar: object
) -> None:
    """
    Calculate the Stage2ar and Stage2marsh for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        stg_sto_ar (object): The stage storage area object.

    Returns:
        None
    """
    model_variables.Stage2ar[i + 2] = stg_sto_ar.stg2ar(
        model_variables.Lake_Stage[i + 1], 0
    )
    model_variables.Stage2marsh[i + 2] = stg_sto_ar.stg2mar(
        model_variables.Lake_Stage[i + 1], 0
    )


def _calculate_rf_and_et(
    i: int, model_variables: object, lo_functions: object, data: object, config: dict
) -> None:
    """
    Calculate the RF and ET for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
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


def _calculate_wsa(
    i: int,
    model_variables: object,
    lo_functions: object,
    df_WSMs: pd.DataFrame,
    wca_stages_df: pd.DataFrame,
    tc_lonino_df: pd.DataFrame,
    data: object,
    config: dict,
) -> None:
    """
    Calculate the WSA for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        df_WSMs (pd.DataFrame): The WSMs DataFrame.
        wca_stages_df (pd.DataFrame): The WCA stages DataFrame.
        tc_lonino_df (pd.DataFrame): The TC Lonino DataFrame.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
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


def _calculate_dsto_and_storage(
    i: int,
    model_variables: object,
    lo_functions: object,
    storage_deviation: list,
    start_storage: float,
) -> None:
    """
    Calculate the DSto and Storage for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        storage_deviation (list): The storage deviation list.
        start_storage (float): The start storage value.

    Returns:
        None
    """
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


def _calculate_lake_stage(
    i: int,
    model_variables: object,
    lo_functions: object,
    stg_sto_ar: object,
    data: object,
    config: dict,
) -> None:
    """
    Calculate the Lake_Stage for the given index.

    Args:
        i (int): The current index.
        model_variables (object): The model variables object.
        lo_functions (object): The functions object.
        stg_sto_ar (object): The stage storage area object.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        None
    """
    model_variables.Lake_Stage[i + 2] = lo_functions.Lake_Stage(
        stg_sto_ar.stg2sto(model_variables.Storage[i + 2], 1),
        data.SFWMM_Daily_Outputs.at[i + 2, "EOD Stg(ft,NGVD)"],
        config["option_stage"],
    )


def _update_lo_model(lo_model: pd.DataFrame, model_variables: object) -> None:
    """
    Update the model DataFrame with the calculated model variables.

    Args:
        lo_model (pd.DataFrame): The model DataFrame.
        model_variables (object): The model variables object.

    Returns:
        None
    """
    lo_model["Cutback"] = model_variables.Cut_back
    lo_model["Stage"] = model_variables.Lake_Stage
    lo_model["Storage"] = model_variables.Storage
    lo_model["S77_Q"] = model_variables.Outlet1USREG
    lo_model["S308_Q"] = model_variables.Outlet2USRG
    lo_model["S77EW"] = model_variables.Outlet1USEWS
    lo_model["TotRegEW"] = model_variables.TotRegEW
    lo_model["TotRegSo"] = model_variables.TotRegSo


def _set_starting_storage(model_variables, start_storage):
    """
    Set the starting storage for the model variables.

    Args:
        model_variables (object): The model variables object.
        start_storage (float): The starting storage value.

    Returns:
        None
    """
    model_variables.Storage[0] = start_storage
    model_variables.Storage[1] = start_storage


def _determine_day_flags(
    startdate: pd.Timestamp, lo_model: pd.DataFrame, begdateCS: pd.Timestamp
) -> str:
    """
    Determine the day flag based on the start date and model.

    Args:
        startdate (pd.Timestamp): The start date.
        lo_model (pd.DataFrame): The model DataFrame.
        begdateCS (pd.Timestamp): The beginning date for CS.

    Returns:
        str: The simulation day.
    """
    if (
        startdate.month == lo_model["date"].iloc[2].month
        and startdate.day == lo_model["date"].iloc[2].day
    ):
        return "SimDay1"
    elif (
        begdateCS.year == lo_model["date"].iloc[2].year
        and begdateCS.month == lo_model["date"].iloc[2].month
        and begdateCS.day == lo_model["date"].iloc[2].day
    ):
        return "CS start date"
    else:
        return lo_model["date"].iloc[2]


def _calculate_targ_stg_df(
    config: dict,
    data: object,
    daily_date_range: pd.DatetimeIndex,
    model_variables: object,
) -> pd.DataFrame:
    """
    Calculate the target stage DataFrame.

    Args:
        config (dict): The configuration dictionary.
        data (object): The data object.
        daily_date_range (pd.DatetimeIndex): The date range.
        model_variables (object): The model variables object.

    Returns:
        pd.DataFrame: The target stage DataFrame.
    """
    targ_stg = (
        data.Targ_Stg_June_1st
        if config["opt_date_targ_stg"] == 1
        else data.Targ_Stg_May_1st
    )

    targ_stg_df = pd.DataFrame(daily_date_range, columns=["dates"])

    # Extract year and day of year (vectorized)
    targ_stg_df["year"] = targ_stg_df["dates"].dt.year
    targ_stg_df["day_of_year"] = targ_stg_df["dates"].dt.dayofyear

    # Compute percentiles in a single vectorized operation
    for perc in [10, 20, 25, 30, 40, 45, 50, 60]:
        targ_stg_df[f"{perc}%"] = targ_stg_df.apply(
            lambda row: replicate(row["year"], row["day_of_year"], perc, targ_stg),
            axis=1,
        )

    return targ_stg_df


def _calculate_late_dry_season(
    adaptive_protocol_df: pd.DataFrame, tc_lonino_df: pd.DataFrame
) -> None:
    """
    Calculate the Late Dry Season (Apr-May) logic and update the adaptive protocol DataFrame.

    Args:
        adaptive_protocol_df (pd.DataFrame): The adaptive protocol DataFrame.
        tc_lonino_df (pd.DataFrame): The TC Lonino DataFrame.

    Returns:
        None
    """
    late_dry_season = [3 < date.month < 6 for date in adaptive_protocol_df["date"]]
    adaptive_protocol_df["Late_Dry_Season"] = late_dry_season
    adaptive_protocol_df["Tributary Hydrologic Condition"] = tc_lonino_df[
        "Tributary_Condition"
    ]


def _calculate_basin_runoff(
    basin_ro: pd.DataFrame,
    data: object,
    outlet1_baseflow: float,
    outlet2_baseflow: float,
    forecast: bool = False,
) -> None:
    """
    Calculate the basin runoff and update the basin runoff DataFrame.

    Args:
        basin_ro (pd.DataFrame): The basin runoff DataFrame.
        data (object): The data object.
        outlet1_baseflow (float): The outlet1 baseflow.
        outlet2_baseflow (float): The outlet2 baseflow.

    Returns:
        None
    """
    # Number of rows
    num_B_R = len(basin_ro)

    # Compute number of days in each month (vectorized)
    basin_ro["Ndays"] = basin_ro["date"].apply(lambda d: monthrange(d.year, d.month)[1])

    # Compute baseflow shortfalls (vectorized) 
    #TODO: keep the date in the monthly flow so that we can compare and make sure it is correct
    basin_ro["BS-C43RO"] = np.maximum(0, outlet1_baseflow - data.C43RO["C43RO"])
    basin_ro["BS-C44RO"] = np.maximum(0, outlet2_baseflow - data.C44RO["C44RO"])

    # Compute C44RO_SLTrib (vectorized)
    # Only keeping SLTRIB if forecast is false
    if not forecast:
        basin_ro["C44RO_SLTRIB"] = basin_ro["BS-C44RO"] + data.SLTRIB["SLTRIB_cfs"]

    # Compute baseflow contribution for C44 (vectorized)
    basin_ro["C44RO-BS"] = (
        np.maximum(0, data.C44RO["C44RO"] - outlet2_baseflow) * basin_ro["Ndays"]
    )

    # Assign direct column values from data
    basin_ro["C43RO"] = data.C43RO["C43RO"]
    basin_ro["C44RO"] = data.C44RO["C44RO"]
    if not forecast:
        basin_ro["SLTRIB"] = data.SLTRIB["SLTRIB_cfs"]


def _initialize_lo_model(
    date_range: pd.DatetimeIndex, data: object, config: dict
) -> pd.DataFrame:
    """
    Initialize the model DataFrame with the given date range, data, and configuration.

    Args:
        date_range (pd.DatetimeIndex): The date range.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        pd.DataFrame: The initialized LOONE model DataFrame.
    """
    lo_model = pd.DataFrame(date_range, columns=["date"])

    lo_model["Net_Inflow"] = data.NetInf_Input["Netflows_acft"]
    lo_model["LOSA_dmd_SFWMM"] = data.SFWMM_W_dmd["LOSA_dmd"] * (
        config["mult_losa"] / 100
    )
    lo_model["C44RO"] = data.C44_Runoff["C44RO"]
    return lo_model


def _determine_seasons(
    daily_date_range: pd.DatetimeIndex, model_variables: object
) -> pd.DataFrame:
    """
    Determine the seasons for each date in the given date range and update the model variables.

    Args:
        daily_date_range (pd.DatetimeIndex): The date range.
        model_variables (object): The model variables object.

    Returns:
        pd.DataFrame: The seasons DataFrame.
    """
    seasons = pd.DataFrame(daily_date_range, columns=["date"])
    seasons_count = len(seasons.index)
    for i in range(seasons_count):
        month = seasons["date"].iloc[i].month
        if 2 < month < 6:
            S = 2
        elif 5 < month < 8:
            S = 3
        elif 7 < month < 11:
            S = 4
        else:
            S = 1
        model_variables.Daily_Seasons[i] = S
        model_variables.Mon[i] = month
    seasons["Season"] = model_variables.Daily_Seasons
    seasons["Month"] = model_variables.Mon
    return seasons


def _calculate_daily_water_demand(
    date_range: pd.DatetimeIndex, startdate: pd.Timestamp, data: object, config: dict
) -> pd.DataFrame:
    """
    Calculate the daily water demand based on the given date range, start date, data, and configuration.

    Args:
        date_range (pd.DatetimeIndex): The date range.
        startdate (pd.Timestamp): The start date.
        data (object): The data object.
        config (dict): The configuration dictionary.

    Returns:
        pd.DataFrame: The water demand DataFrame with daily demand.
    """
    water_demand = pd.DataFrame(date_range, columns=["date"])

    # Compute 'count' using a cumulative sum
    water_demand["count"] = (water_demand["date"].dt.month != date_range[0].month) | (
        water_demand["date"].dt.day != date_range[0].day
    )
    water_demand["count"] = water_demand["count"].cumsum()

    # Compute 'Week_num' directly
    water_demand["Week_num"] = water_demand["count"].floordiv(7).add(1)
    water_demand.loc[water_demand["count"] > 363, "Week_num"] = 52  # Cap at 52

    # Compute 'Daily_demand' using vectorized operations
    weekly_demand = data.Weekly_dmd[f'C{config["code"]}']
    water_demand["Daily_demand"] = (
        weekly_demand.iloc[water_demand["Week_num"] - 1].values / 7
    ) * (config["multiplier"] / 100)

    return water_demand


def _define_start_and_end_dates(config: dict) -> tuple:
    """
    Define the start and end dates based on the configuration.

    Args:
        config (dict): The configuration dictionary.

    Returns:
        tuple: A tuple containing the start date, begdateCS, and end date.
    """
    year, month, day = map(int, config["start_date_entry"])
    startdate = datetime(year, month, day).date()
    year, month, day = map(int, config["start_date_entry"])
    begdateCS = datetime(year, month, day).date()
    year, month, day = map(int, config["end_date_entry"])
    enddate = datetime(year, month, day).date()

    return startdate, begdateCS, enddate


def _create_dectree_df(
    daily_date_range: pd.DatetimeIndex, tc_lonino_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Create the dectree DataFrame with the given date range and TC Lonino DataFrame.

    Args:
        daily_date_range (pd.DatetimeIndex): The date range.
        tc_lonino_df (pd.DataFrame): The TC Lonino DataFrame.

    Returns:
        pd.DataFrame: The decision tree DataFrame.
    """
    dec_tree_df = pd.DataFrame(daily_date_range, columns=["Date"])
    dec_tree_df["Zone_B_MetFcast"] = tc_lonino_df["LONINO_Seasonal_Classes"]
    return dec_tree_df


def _calculate_pulse_averages(data: object, config: dict) -> dict:
    """
    Calculate the pulse averages for S-80 and S-77.

    Args:
        data (object): The data object containing pulse information.
        config (dict): The configuration dictionary.

    Returns:
        dict: A dictionary containing the pulse averages for S-80 and S-77.
    """
    pulse_averages = {
        "s80avg_l1": data.Pulses[f"S-80_L1"].mean(),
        "s80avg_l2": data.Pulses[f"S-80_L2"].mean(),
        "s80avg_l3": data.Pulses[f"S-80_L3"].mean(),
        "s77avg_l1": data.Pulses[f"S-77_L1"].mean(),
        "s77avg_l2": data.Pulses[f"S-77_L2"].mean(),
        "s77avg_l3": data.Pulses[f"S-77_L3"].mean(),
    }
    return pulse_averages


def _initialize_model_variables_stage_levels_flags(
    model_variables: object,
    config: dict,
    startdate: pd.Timestamp,
    lo_model: pd.DataFrame,
    begdateCS: pd.Timestamp,
    forecast: bool = False,
    workspace: str = "",
) -> None:
    """
    Initialize the model variables with the given configuration and start date.

    Args:
        model_variables (object): The model variables object.
        config (dict): The configuration dictionary.
        startdate (pd.Timestamp): The start date.
        lo_model (pd.DataFrame): The LOONE model DataFrame.
        begdateCS (pd.Timestamp): The beginning date for CS.
        forecast (bool, optional): Whether to run in forecast mode. Defaults to False.
        workspace (str, optional): The path to the workspace directory. Needed in forecast mode Defaults to "".
    Returns:
        None
    """
    # TODO - is this where the stage should be fixed?
    if forecast:
        stage = pd.read_csv(os.path.join(workspace, "LO_Stage.csv"))
        model_variables.Lake_Stage[0] = stage["Average_Stage"].iloc[-1]
        model_variables.Lake_Stage[1] = stage["Average_Stage"].iloc[-1]
    else:
        model_variables.Lake_Stage[0] = config["beg_stage_cs"]
        model_variables.Lake_Stage[1] = config["beg_stage_cs"]
    model_variables.DecTree_Relslevel[0] = np.nan
    model_variables.DecTree_Relslevel[1] = np.nan
    model_variables.DayFlags[2] = _determine_day_flags(startdate, lo_model, begdateCS)


def LOONE_Q(
    workspace: str,
    sensitivity_analysis_params: List[float] = [],
    optimization_params: dict[float] = None,
    forecast: bool = False,
    ensemble: int = None,
) -> None:
    """This function runs the LOONE Q module.

    Args:
        workspace (str): The path to the workspace directory.
        sensitivity_analysis_params (List[float], optional): Holds a release rate for each month of the year.
            Required when sim_type is 3. Defaults to [].
        optimization_params (dict[float], optional): Holds the optimization parameters,
            including p1, p2, s77_dv, s308_dv, and tp_lake_s. Required when sim_type is 2. Defaults to None.
        forecast (bool, optional): Whether to run in forecast mode. Defaults to False.
        ensemble (int, optional): The ensemble number. Defaults to None. Only necesary for forecast mode.

    Returns:
        None
    """
    os.chdir(workspace)
    config = load_config(workspace)

    if optimization_params:
        p1 = optimization_params["p1"]
        p2 = optimization_params["p2"]
        s77_dv = optimization_params["s77_dv"]
        s308_dv = optimization_params["s308_dv"]
        tp_lake_s = optimization_params["tp_lake_s"]
    else:
        p1 = p2 = s77_dv = s308_dv = tp_lake_s = 0

    data = DClass(workspace, forecast, ensemble)
    model_variables = MVarClass(config, forecast)
    print("LOONE Q Module is Running!")
    # Based on the defined Start and End year, month, and day on the
    # Pre_defined_Variables File, Startdate and enddate are defined.
    startdate, begdateCS, enddate = _define_start_and_end_dates(config)

    ###################################################################
    if config["sim_type"] in [0, 1, 3]:
        df_wsms.WSMs(workspace, forecast, ensemble)

    df_WSMs = pd.read_csv("df_WSMs.csv")

    # The Following Code interpolates daily LOSA demand from weekly
    # data for 6 differnet datasets where the user defines the LOSA
    # demand that will be used based on a Code (1:6).
    # Set time frame for model run
    if forecast:
        today_date = datetime.today()
        future_date = today_date + timedelta(days=15)
        daily_date_range = pd.date_range(start=today_date, end=future_date, freq="D")
    else:
        daily_date_range = pd.date_range(start=startdate, end=enddate, freq="D")

    # date_range_2 = pd.date_range(start=startdate, end=enddate, freq="D")
    # water_demand = _calculate_daily_water_demand(date_range_2, startdate, data, config)
    water_demand = _calculate_daily_water_demand(
        daily_date_range, startdate, data, config
    )

    ###################################################################
    # Determine Tributary Hydrologic Conditions
    tc_lonino_df = trib_hc.Trib_HC(workspace, forecast, ensemble)
    # Determine WCA Stages
    wca_stages_df = WCA_Stages_Cls(workspace, tc_lonino_df, forecast)
    # A dataframe to determine eachday's season (Months 11,12,1,2 are
    # Season 1, Months 3,4,5 are season 2, Months 6,7 are season 3,
    # Months 8,9,10 are season 4 )
    # date_range_5 = pd.date_range(start=startdate, end=enddate, freq="D")
    # seasons = _determine_seasons(date_range_5, model_variables)
    seasons = _determine_seasons(daily_date_range, model_variables)

    ###################################################################
    # This following Script runs the main model daily simulations.
    if forecast:
        date_range_one_day = pd.date_range(
            start=today_date - timedelta(days=1),
            end=future_date,
            freq="D",
        )
    else:
        date_range_one_day = pd.date_range(
            start=startdate - timedelta(days=1),
            end=enddate,
            freq="D",
        )
    lo_model = _initialize_lo_model(date_range_one_day, data, config)
    n_rows = len(lo_model.index)

    ##################################
    # dec_tree_df = _create_dectree_df(date_range_5, tc_lonino_df)
    dec_tree_df = _create_dectree_df(daily_date_range, tc_lonino_df)

    # Create a dataframe that includes Monthly Mean Basin Runoff &
    ##  BaseFlow-Runoff & Runoff-Baseflow (cfs)
    if forecast:
        new_startdate = today_date
        if new_startdate.day > 1:
            new_startdate = new_startdate.replace(day=1)
        monthly_date_range = pd.date_range(
            start=new_startdate, end=future_date, freq="MS"
        )
    else:
        monthly_date_range = pd.date_range(start=startdate, end=enddate, freq="MS")
    if monthly_date_range.empty:
        monthly_date_range = pd.DatetimeIndex(
            [pd.to_datetime(startdate).replace(day=1)]
        )
    # date_range_11 = pd.date_range(start=startdate, end=enddate, freq="MS")
    # date_range_11d = pd.date_range(start=startdate, end=enddate, freq="D")
    # date_range_11d.name = "Date"
    basin_ro = pd.DataFrame(monthly_date_range, columns=["date"])
    # Baseflows
    outlet1_baseflow = data.S77_RegRelRates["Zone_D0"].iloc[0]
    outlet2_baseflow = data.S80_RegRelRates["Zone_D0"].iloc[0]
    _calculate_basin_runoff(basin_ro, data, outlet1_baseflow, outlet2_baseflow, forecast)
    lo_model["C43RO"] = data.C43RO_Daily["C43RO"]
    pulse_averages = _calculate_pulse_averages(data, config)
    s80avg_l1 = pulse_averages["s80avg_l1"]
    s80avg_l2 = pulse_averages["s80avg_l2"]
    s80avg_l3 = pulse_averages["s80avg_l3"]
    s77avg_l1 = pulse_averages["s77avg_l1"]
    s77avg_l2 = pulse_averages["s77avg_l2"]
    s77avg_l3 = pulse_averages["s77avg_l3"]
    basin_ro = basin_ro.set_index(["date"])
    basin_ro.index = pd.to_datetime(basin_ro.index)
    # basin_ro_daily = basin_ro.reindex(date_range_11d, method="ffill")
    basin_ro_daily = basin_ro.reindex(daily_date_range, method="ffill")
    basin_ro = basin_ro.reset_index()
    vlookup1 = basin_ro_daily["BS-C44RO"]
    vlookup1_c = [x for x in vlookup1 if ~np.isnan(x)]
    ###################################################################
    # This following script contains the logic and calculations for
    # the proposed Lake Okeechobee Adaptive Protocol.
    # adaptive_protocol_df = pd.DataFrame(date_range_5, columns=["date"])
    adaptive_protocol_df = pd.DataFrame(daily_date_range, columns=["date"])
    _calculate_late_dry_season(adaptive_protocol_df, tc_lonino_df)

    # targ_stg_df = _calculate_targ_stg_df(config, data, date_range_5, model_variables)
    targ_stg_df = _calculate_targ_stg_df(
        config, data, daily_date_range, model_variables
    )

    # Outlet1_baseflow = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    outlet1_baseflow = 450  # cfs
    vlookup2 = basin_ro_daily["BS-C43RO"]
    vlookup2_c = [x for x in vlookup2 if ~np.isnan(x)]
    ###################################################################
    _initialize_model_variables_stage_levels_flags(
        model_variables, config, startdate, lo_model, begdateCS, forecast, workspace
    )
    if forecast:
        stage = pd.read_csv(os.path.join(workspace, "LO_Stage.csv"))
        start_stage = stage["Average_Stage"].iloc[-1]
        start_storage = stg_sto_ar.stg2sto(start_stage, 0)
    else:
        start_storage = stg_sto_ar.stg2sto(config["start_stage"], 0)
    _set_starting_storage(model_variables, start_storage)
    # Flood = np.zeros(n_rows, dtype = object)
    ##Here, I will insert the Storage Deviaiton Values as Input!
    # storage deviation is always 0 in forecast mode
    if forecast:
        storage_deviation = [0] * n_rows
    else:
        storage_deviation = data.Storage_dev_df["DS_dev"]
    # Create a Choose Function for AP Post Baseflow
    # if Pre_defined_Variables.Opt_AdapProt == 0:
    #     C = 450
    # elif Pre_defined_Variables.Opt_AdapProt == 1:
    #     C = Data.S77_RegRelRates['Zone_D0'].iloc[0]
    # Choose_1 = C
    choose_1 = 450  # cfs
    ###################################################################
    _calculate_initial_zone_code_and_lo_zone(model_variables, lo_functions, df_WSMs)
    for i in range(n_rows - 2):
        _calculate_wsm_zone(i, model_variables, lo_functions, df_WSMs)
        _calculate_max_supply(i, model_variables, lo_functions, water_demand, config)
        _calculate_losa_supply(i, model_variables, lo_functions, lo_model, config)
        _calculate_ni_supply(i, model_variables, lo_model)
        _calculate_cutback(i, model_variables, lo_model)
        _calculate_demand_not_supplied(i, model_variables, lo_model)
        _calculate_zone_code(i, model_variables, lo_functions, df_WSMs)
        _generate_zone_column(i, model_variables, lo_functions)
        _calculate_decision_tree_variables(
            i,
            model_variables,
            dec_tree_functions,
            tc_lonino_df,
            config,
            seasons,
            dec_tree_df,
            df_WSMs,
        )
        _calculate_release_levels(
            i, startdate, lo_model, model_variables, lo_functions, config, tc_lonino_df
        )
        _calculate_zone_codes(i, model_variables, lo_functions, df_WSMs)
        _calculate_release_level_codes(i, model_variables, lo_functions, config)
        _calculate_outlet2ds_mult(
            i, lo_model, model_variables, lo_functions, seasons, config
        )
        _calculate_outlet2dsrs(
            i,
            lo_model,
            model_variables,
            lo_functions,
            data,
            s80avg_l1,
            s80avg_l2,
            s80avg_l3,
            config,
        )
        _calculate_outlet2dsbs(
            i, model_variables, lo_functions, vlookup1_c, outlet2_baseflow, config
        )
        _calculate_outlet2usbk(
            i, model_variables, lo_functions, df_WSMs, lo_model, data, config
        )
        _calculate_outlet2usrg_code(
            i,
            model_variables,
            lo_functions,
            data,
            config,
            p1,
            p2,
            s308_dv,
            tp_lake_s,
            lo_model,
            date_range_one_day,
        )
        _calculate_outlet2ds(i, model_variables, lo_functions, data, config)
        _calculate_relevel_code_3_s77(i, model_variables, lo_functions, config)
        _calculate_outlet1us_mult(i, model_variables, lo_functions, seasons, config)
        _calculate_outlet1us_mult_2(i, model_variables, lo_functions, lo_model, config)
        _calculate_outlet1usrs(
            i,
            model_variables,
            lo_functions,
            lo_model,
            data,
            config,
            s77avg_l1,
            s77avg_l2,
            s77avg_l3,
        )
        _calculate_sum_outlet1usrs(i, model_variables, lo_functions, lo_model)
        _calculate_outlet1usbk(i, model_variables, lo_functions, lo_model, data, config)
        _calculate_rowest(i, model_variables, lo_model)
        _calculate_outlet1dsbs(
            i, model_variables, lo_functions, vlookup2_c, outlet1_baseflow, config
        )
        _calculate_outlet1usbs(i, model_variables, lo_functions, config)
        _define_thc_class_normal_or_above(
            i,
            n_rows,
            model_variables,
            config,
            adaptive_protocol_df,
            df_WSMs,
            targ_stg_df,
            data,
            choose_1,
        )
        _calculate_outlet1usbsap(i, model_variables, lo_functions, config)
        _calculate_outlet1usews(i, model_variables, lo_functions, data, config)
        _calculate_outlet1usreg(
            i,
            model_variables,
            lo_functions,
            data,
            config,
            p1,
            p2,
            s77_dv,
            tp_lake_s,
            date_range_one_day,
            lo_model,
            sensitivity_analysis_params,
        )
        _calculate_outlet1ds(i, model_variables, lo_functions, data, config)
        _calculate_tot_reg_ew(i, model_variables)
        _calculate_choose_wca(i, model_variables, lo_functions, data, config)
        _calculate_choose_l8c51(i, model_variables, lo_functions, data, config)
        _calculate_tot_reg_so(i, model_variables)
        _calculate_stage2ar_and_stage2marsh(i, model_variables, stg_sto_ar)
        _calculate_rf_and_et(i, model_variables, lo_functions, data, config)
        _calculate_wsa(
            i,
            model_variables,
            lo_functions,
            df_WSMs,
            wca_stages_df,
            tc_lonino_df,
            data,
            config,
        )
        _calculate_dsto_and_storage(
            i, model_variables, lo_functions, storage_deviation, start_storage
        )
        _calculate_lake_stage(
            i, model_variables, lo_functions, stg_sto_ar, data, config
        )
        _update_lo_model(lo_model, model_variables)

    # Write out the results to a file - Needed because
    # execute_loone.py calls this script as a subprocess and can't get
    # this data.
    if forecast:
        #these should be in cfs or acft depending on the column
        lo_model.to_csv(f"LOONE_Q_Outputs_{ensemble:02d}.csv")
    else:
        lo_model.to_csv("LOONE_Q_Outputs.csv")

    # Add scenario data to output df
    df_stage = pd.DataFrame(data=model_variables.Lake_Stage, columns=["Outputs"])
    df_caloosahatchee = pd.DataFrame(
        data=model_variables.Outlet1USREG, columns=["Outputs"]
    )
    df_saint_lucie = pd.DataFrame(data=model_variables.Outlet2USRG, columns=["Outputs"])
    df_south = pd.DataFrame(data=model_variables.TotRegSo / 1.9835, columns=["Outputs"])
    df_out = pd.concat([df_stage, df_caloosahatchee, df_saint_lucie, df_south])

    return [lo_model, df_out.T.loc["Outputs"]]


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "workspace",
        nargs=1,
        help="The path to the working directory.",
    )
    argparser.add_argument(
        "--sensitivity_analysis_params",
        nargs="*",
        type=float,
        help="Optional sensitivity analysis parameters: release rates for each month of the year",
    )

    argparser.add_argument(
        "--optimization_params",
        nargs="*",
        type=float,
        default=[0, 0, 0, 0, 0],
        help="Optional optimization parameters: p1, p2, s77_dv, s308_dv, tp_lake_s",
    )

    args = argparser.parse_args()
    workspace = args.workspace[0]
    all_args = [workspace]

    if args.sensitivity_analysis_params:
        sensitivity_analysis_params = args.sensitivity_analysis_params
        all_args.append(sensitivity_analysis_params)

    if args.optimization_params and len(args.optimization_params) == 5:
        optimization_params = {
            "p1": args.optimization_params[0],
            "p2": args.optimization_params[1],
            "s77_dv": args.optimization_params[2],
            "s308_dv": args.optimization_params[3],
            "tp_lake_s": args.optimization_params[4],
        }
        all_args.append(optimization_params)

    LOONE_Q(*all_args)
