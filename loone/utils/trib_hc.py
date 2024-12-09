import os
import pandas as pd
import numpy as np
from datetime import datetime
from loone.data.model_variables import M_var as MVarClass
from loone.utils import load_config, lonino_functions
from loone.data import Data as DClass


# I determine daily values for the Tributary conditions and Seasonal/Multi-Seasonal LONINO classes
# using a weekly Trib. Condition data and Monthly LONINO data.
def Trib_HC(workspace: str):
    """
    This function generates daily values for the Tributary conditions and Seasonal/Multi-Seasonal LONINO classes
    using a weekly Trib. Condition data and Monthly LONINO data.

    It generates the Tributary Condition Dataframe, Trib_Cond_df, and the LONINO Dataframe, LONINO_df.
    It uses the Trib_Cond_df to calculate the Tributary Condition Index (TCI) and the LONINO_df to calculate
    the LONINO seasonal and multi-seasonal classes.

    Finally, it generates a daily date range, TC_LONINO_df, and assigns the TCI, LONINO seasonal, and LONINO
    multi-seasonal classes to this dataframe.

    Returns:
        TC_LONINO_df (pd.DataFrame): A dataframe with the daily values of the Tributary conditions and Seasonal/Multi-Seasonal LONINO classes.

    """
    os.chdir(workspace)
    config = load_config(workspace)
    Data = DClass(workspace)
    M_var = MVarClass(config)
    # Generate weekly time step date column where frequency is 'W-Fri' to start on 01/01/2008.
    # FIXME: Always check here for start date, end date, and frequency to match with the Trib. Condition weekly data obtained.
    year, month, day = map(int, config["start_date_entry"])
    startdate = datetime(year, month, day).date()
    year, month, day = map(int, config["end_date_entry"])
    enddate = datetime(year, month, day).date()
    year, month, day = map(int, config["end_date_tc"])
    enddate_TC = datetime(year, month, day).date()
    date_rng_3 = pd.date_range(start=startdate, end=enddate_TC, freq="W-Fri")
    # Generate the Tributary Condition Dataframe.
    Trib_Cond_df = pd.DataFrame(date_rng_3, columns=["date"])
    TC_Count = len(Trib_Cond_df.index)
    for i in range(TC_Count):
        M_var.RF_Cls[i] = lonino_functions.RF_Cls(
            Data.Wkly_Trib_Cond["NetRF"].iloc[i]
        )
        M_var.MainTrib_Cls[i] = lonino_functions.MainTrib_Cls(
            Data.Wkly_Trib_Cond["S65E"].iloc[i]
        )
        M_var.Palmer_Cls[i] = lonino_functions.Palmer_Cls(
            Data.Wkly_Trib_Cond["Palmer"].iloc[i]
        )
        M_var.NetInflow_Cls[i] = lonino_functions.NetInflow_Cls(
            Data.Wkly_Trib_Cond["NetInf"].iloc[i]
        )
        M_var.Max_RF_MainTrib[i] = max(M_var.RF_Cls[i], M_var.MainTrib_Cls[i])
        M_var.Max_Palmer_NetInf[i] = max(
            M_var.Palmer_Cls[i], M_var.NetInflow_Cls[i]
        )
    if config["tci"] == 1:  # Tributary Condition Index
        Trib_Cond_df["TCI"] = M_var.Max_Palmer_NetInf
    else:
        Trib_Cond_df["TCI"] = M_var.Max_RF_MainTrib
    # Generate a monthly time step date column
    date_rng_4 = pd.date_range(start=startdate, end=enddate, freq="MS")
    # Create a LONINO Dataframe
    LONINO_df = pd.DataFrame(date_rng_4, columns=["date"])
    LONINO_Count = len(LONINO_df.index)
    for i in range(LONINO_Count):
        M_var.Seas[i] = Data.LONINO_Seas_data[
            str(LONINO_df["date"].iloc[i].month)
        ].iloc[LONINO_df["date"].iloc[i].year - config["start_year"]]
        M_var.M_Seas[i] = Data.LONINO_Mult_Seas_data[
            str(LONINO_df["date"].iloc[i].month)
        ].iloc[LONINO_df["date"].iloc[i].year - config["start_year"]]
    LONINO_df["LONINO_Seas"] = M_var.Seas
    LONINO_df["LONINO_Mult_Seas"] = M_var.M_Seas
    for i in range(config["month_n"]):
        M_var.LONINO_Seas_cls[i] = lonino_functions.LONINO_Seas_cls(
            LONINO_df["LONINO_Seas"].iloc[i]
        )
        M_var.LONINO_M_Seas_cls[i] = lonino_functions.LONINO_M_Seas_cls(
            LONINO_df["LONINO_Mult_Seas"].iloc[i]
        )
    LONINO_df["LONINO_Seasonal_Cls"] = M_var.LONINO_Seas_cls
    LONINO_df["LONINO_Mult_Seasonal_Cls"] = M_var.LONINO_M_Seas_cls
    # Generate a daily date range
    date_rng_5 = pd.date_range(start=startdate, end=enddate, freq="D")
    TC_LONINO_df = pd.DataFrame(date_rng_5, columns=["Date"])
    row_nm = len(TC_LONINO_df.index)
    Trib_Cond = np.zeros(row_nm)
    for i in range(row_nm):
        Trib_Cond[i] = Trib_Cond_df["TCI"].iloc[int(i / 7)]
    TC_LONINO_df["Tributary_Condition"] = Trib_Cond
    data_S_MS = [
        LONINO_df["date"],
        LONINO_df["LONINO_Seasonal_Cls"],
        LONINO_df["LONINO_Mult_Seasonal_Cls"],
    ]
    headers_S_MS = ["date", "LONINO_Seasonal_Class", "LONINO_MSeasonal_Class"]
    LONINO_Seas_MSeas_df = pd.concat(data_S_MS, axis=1, keys=headers_S_MS)
    LONINO_Seas_MSeas_df = LONINO_Seas_MSeas_df.set_index("date")
    LONINO_Seas_MSeas_df = LONINO_Seas_MSeas_df.reindex(
        date_rng_5, method="ffill"
    )
    LONINO_Seas_MSeas_df = LONINO_Seas_MSeas_df.reset_index()
    TC_LONINO_df["LONINO_Seasonal_Classes"] = LONINO_Seas_MSeas_df[
        "LONINO_Seasonal_Class"
    ]
    TC_LONINO_df["LONINO_MultiSeasonal_Classes"] = LONINO_Seas_MSeas_df[
        "LONINO_MSeasonal_Class"
    ]
    return TC_LONINO_df
