# -*- coding: utf-8 -*-
"""
Created on Tue May 17 00:37:38 2022

@author: osama
"""
# This Script Interpolates each Water Shortage Management (WSMs) and each Regulation Schedule Breakpoint Zone (D, C, B, and A).
import pandas as pd
from datetime import datetime
import numpy as np
from scipy import interpolate
import os
from Model_Config import Model_Config

Working_Path = Model_Config.Working_Path
os.chdir("%s" % Working_Path)
from Pre_defined_Variables import Pre_defined_Variables
from Data import Data

year, month, day = map(int, Pre_defined_Variables.startdate_entry)
startdate = datetime(year, month, day).date()
year, month, day = map(int, Pre_defined_Variables.startdate_entry)
begdateCS = datetime(year, month, day).date()
year, month, day = map(int, Pre_defined_Variables.enddate_entry)
enddate = datetime(year, month, day).date()


def WSMs():
    # Set time frame for model run such that it starts on the defined startdate but ends on 1/1/(endyear+1)
    # date_rng_1 = pd.date_range(start = startdate, end = '1/1/%d'%(Pre_defined_Variables.endyear+1), freq= 'D')
    date_rng_1 = pd.date_range(
        start=startdate, end="4/1/%d" % (Pre_defined_Variables.endyear), freq="D"
    )
    # Create a data frame with a date column
    df_WSMs = pd.DataFrame(date_rng_1, columns=["date"])
    # Generate an annual cumulative day count
    WSM_length = len(df_WSMs.index)
    Oper_Zones = list(Data.WSMs_RSBKs)
    Oper_Zones.remove("Date")
    Oper_Zones.remove("Day")
    for i in Oper_Zones:
        globals()[i] = np.zeros(WSM_length)
    WSM_Count = date_rng_1.strftime("%j")
    for i in range(WSM_length):
        for j in Oper_Zones:
            globals()[j][i] = interpolate.interp1d(
                Data.WSMs_RSBKs["Day"], Data.WSMs_RSBKs["%s" % j], kind="linear"
            )(WSM_Count[i])

    df_WSMs["count"] = WSM_Count
    for j in Oper_Zones:
        df_WSMs["%s" % j] = globals()[j]
    if Pre_defined_Variables.Opt_NewTree == 1:
        df_WSMs["C-b"] = Data.WSMs_RSBKs["C-b_NewTree"]
    else:
        df_WSMs["C-b"] = Data.WSMs_RSBKs["C-b_NoNewTree"]
    df_WSMs.drop(["C-b_NewTree", "C-b_NoNewTree"], axis=1, inplace=True)

    # Add one row in top of the dataframe (i.e. Dec/31/previous year) where the values = values of the original first row in the dataframe!
    First_row = pd.DataFrame(
        {
            "WSM4": df_WSMs["WSM4"].iloc[0],
            "WSM3": df_WSMs["WSM3"].iloc[0],
            "WSM2": df_WSMs["WSM2"].iloc[0],
            "WSM1": df_WSMs["WSM1"].iloc[0],
            "D0": df_WSMs["D0"].iloc[0],
            "D1": df_WSMs["D1"].iloc[0],
            "D2": df_WSMs["D2"].iloc[0],
            "D3": df_WSMs["D3"].iloc[0],
            "C": df_WSMs["C"].iloc[0],
            "B": df_WSMs["B"].iloc[0],
            "A": df_WSMs["A"].iloc[0],
            "C-b": df_WSMs["C-b"].iloc[0],
        },
        index=[0],
    )
    df_WSMs = pd.concat([First_row, df_WSMs]).reset_index(drop=True)

    df_WSMs.to_csv(os.path.join(Working_Path, "df_WSMs.csv"))
