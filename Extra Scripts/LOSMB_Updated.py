# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 11:14:10 2021

@author: osama
"""

# This Script represents a simple Mass Balance Model of Lake Okeechobee Inflows/Outflows!
# First Import Needed Packages
import pandas as pd
import numpy as np
import os

# Directory where the Script loactes
os.chdir("C:/Work/Research/Data Analysis/Tools/Python_Scripts")
from Data_Analyses_Fns import *

Working_dir = "C:/Work/Research/LOONE"
os.chdir("%s" % Working_dir)
import utils.stg_sto_ar

Working_dir = "C:/Work/Research/LOONE/LOONE_Data_Pre"
os.chdir("%s" % Working_dir)

# Read LO Flows Data File
Csv_Date_Range
Csv_Date_Range("Flow_df.csv", 2007, 10, 1, 2022, 12, 31)
Flow_df = pd.read_csv("./Flow_df_3MLag.csv")
# Read Netflows
NetFlows = Flow_df["Netflows"] / 1233  # m3 to acft (per day)
# Read Outflows
Outflows = Flow_df["Outflows"] / 1233
# Read S77 and S308 Outflows
S77Backflow = Flow_df["S77_In"] / 1233
S308Backflow = Flow_df["S308_In"] / 1233
# Read LO Stage Storage SA
Csv_Date_Range("LO_Stg_Sto_SA.csv", 2007, 10, 1, 2022, 12, 31)
LO_Stg_Sto_SA = pd.read_csv("./LO_Stg_Sto_SA_2008-2023.csv")
# Determine Rainfall Volume ac-ft
# Read Average Rainfall (In)
Csv_Date_Range("LAKE_RAINFALL_DATA.csv", 2007, 10, 1, 2022, 12, 31)  # inch
Rainfall_Z = pd.read_csv("./LAKE_RAINFALL_DATA_2008-2023.csv")
RFVolume = Rainfall_Z["average_rainfall"] * LO_Stg_Sto_SA["SA_acres"] / 12  # acft
# Read Average ETP (In)
Csv_Date_Range("LOONE_AVERAGE_ETPI_DATA.csv", 2007, 10, 1, 2022, 12, 31)  # inch
ETP_Z = pd.read_csv("./LOONE_AVERAGE_ETPI_DATA_2008-2023.csv")
ETVolume = ETP_Z["average_ETPI"] * LO_Stg_Sto_SA["SA_acres"] / 12  # acft

num_days = len(NetFlows)
date_range = Flow_df["date"]
NetFlows = NetFlows.fillna(0)
Outflows = Outflows.fillna(0)
S308BK = S308Backflow.fillna(0)
S77BK = S77Backflow.fillna(0)
RFVOL = RFVolume.fillna(0)
ETVOL = ETVolume.fillna(0)

# Exploring Stroage Deviations between Stroage calculated from Continuity Equation and Stroage observed (From Observed Stage)
# Read Observed Storage(acft)
LOSto_obs = LO_Stg_Sto_SA["Storage_acft"]
# Calculate deviation from Zero!
DS_cteq = np.zeros(num_days)
DS_obs = np.zeros(num_days)
DS_dev = np.zeros(num_days)
for i in range(1, num_days):
    DS_cteq[i] = (
        NetFlows[i] + RFVOL[i] - ETVOL[i] + (S77BK[i] + S308BK[i]) - Outflows[i]
    )
    DS_obs[i] = LOSto_obs[i] - LOSto_obs[i - 1]
    DS_dev[i] = DS_obs[i] - DS_cteq[i]
Storage_Exploring = pd.DataFrame(date_range, columns=["date"])
Storage_Exploring["DS_cteq"] = DS_cteq
Storage_Exploring["DS_obs"] = DS_obs
Storage_Exploring["DS_dev"] = DS_dev
Storage_Exploring.to_csv("./Storage_Exploring_3MLag.csv")

# #Read Delta Storage Deviations!
# Storage_Dev_df = pd.read_csv('./Storage_Exploring.csv')
# Storage_dev = Storage_Dev_df['DS_dev']

# DS = np.zeros(num_days)
# Storage = np.zeros(num_days)
# Storage[0] = LOSto_obs[0] #acft
# Stage = np.zeros(num_days)
# Stage[0] = LO_Stg_Sto_SA['Stage_ft'].iloc[0] #ft
# for i in range(1,num_days):
#     DS[i] = NetFlows[i] + RFVOL[i] - ETVOL[i] + (S77BK[i] + S308BK[i]) - Outflows[i] + Storage_dev[i]
#     Storage[i] = Storage[i-1] + DS[i]
#     Stage[i] = utils.stg_sto_ar.stg2sto(Storage[i],1)
# LOSMB = pd.DataFrame(date_range, columns=['date'])
# LOSMB['Storage_acft'] = Storage
# LOSMB['Stage_ft'] = Stage
# LOSMB.to_csv('./LOSMB_Outputs.csv')
