# This Script Interpolates each Water Shortage Management (WSMs) and each Regulation Schedule Breakpoint Zone (D, C, B, and A).
import os
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from scipy import interpolate
from loone.data import Data as DClass
from loone.utils import load_config


def WSMs(workspace: str, forecast: bool = False, ensemble: int = None) -> None:
    """Generate WSMs (Weather State Modifiers) based on the given workspace.

        Args:
            workspace (str): The path to the workspace.
    
        Returns:
            None

        Raises:
            None
    """

    os.chdir(workspace)
    config = load_config(workspace)
    data = DClass(workspace, forecast, ensemble)
    if forecast:
        today = datetime.today().date()
        start_date = today
        end_date = today + timedelta(days=16)
    else:
        year, month, day = map(int, config["start_date_entry"])
        start_date = datetime(year, month, day).date()
        year, month, day = map(int, config["end_date_entry"])
        # Config end date + 1 day
        end_date = datetime(year, month, day).date() + timedelta(days=1)

    date_rng_1 = pd.date_range(start=start_date, end=end_date, freq="D")
    # Create a data frame with a date column
    df_wsms = pd.DataFrame(date_rng_1, columns=["date"])
    # Generate an annual cumulative day count
    wsm_length = len(df_wsms.index)
    operational_zones = list(data.WSMs_RSBKs)
    operational_zones.remove("Date")
    operational_zones.remove("Day")
    for i in operational_zones:
        globals()[i] = np.zeros(wsm_length)
    WSM_Count = date_rng_1.strftime("%j")
    for i in range(wsm_length):
        for j in operational_zones:
            globals()[j][i] = interpolate.interp1d(
                data.WSMs_RSBKs["Day"],
                data.WSMs_RSBKs[j],
                kind="linear",
            )(WSM_Count[i])

    df_wsms["count"] = WSM_Count
    for j in operational_zones:
        df_wsms[j] = globals()[j]
    if config["opt_new_tree"] == 1:
        df_wsms["C-b"] = data.WSMs_RSBKs["C-b_NewTree"]
    else:
        df_wsms["C-b"] = data.WSMs_RSBKs["C-b_NoNewTree"]
    df_wsms.drop(["C-b_NewTree", "C-b_NoNewTree"], axis=1, inplace=True)

    # Add one row in top of the dataframe (i.e. Dec/31/previous year) where the values = values of the original first row in the dataframe!
    first_row = pd.DataFrame(
        {
            "WSM4": df_wsms["WSM4"].iloc[0],
            "WSM3": df_wsms["WSM3"].iloc[0],
            "WSM2": df_wsms["WSM2"].iloc[0],
            "WSM1": df_wsms["WSM1"].iloc[0],
            "D0": df_wsms["D0"].iloc[0],
            "D1": df_wsms["D1"].iloc[0],
            "D2": df_wsms["D2"].iloc[0],
            "D3": df_wsms["D3"].iloc[0],
            "C": df_wsms["C"].iloc[0],
            "B": df_wsms["B"].iloc[0],
            "A": df_wsms["A"].iloc[0],
            "C-b": df_wsms["C-b"].iloc[0],
        },
        index=[0],
    )
    df_wsms = pd.concat([first_row, df_wsms]).reset_index(drop=True)

    df_wsms.to_csv("df_WSMs.csv")
