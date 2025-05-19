import numpy as np
import pandas as pd
from datetime import datetime, timedelta


class M_var:
    """Class to represents model variables."""
    def __init__(self, config: dict, forecast: bool = False):
        """
        Initializes the M_var class with model variables.

        This constructor sets up various model variables and arrays based on the provided configuration dictionary.
        It calculates date ranges and initializes arrays for different time scales (daily, weekly, monthly)
        and model parameters, including tributary conditions, seasonal classes, lake stages, supply, and outflows.

        Args:
            config (dict): A dictionary containing configuration parameters, including:
                - "start_date_entry": A list of integers [year, month, day] for the start date.
                - "end_date_entry": A list of integers [year, month, day] for the end date.
                - "end_date_tc": A list of integers [year, month, day] for the end date of tributary conditions.
                - "month_n": An integer representing the number of months for the LONINO seasonal classes.
        """
        if forecast==True:
            today = datetime.today().date()
            startdate = today
            enddate = today + timedelta(days=15)
            enddate_TC = enddate
        else:
            year, month, day = map(int, config["start_date_entry"])
            startdate = datetime(year, month, day).date()
            year, month, day = map(int, config["end_date_entry"])
            enddate = datetime(year, month, day).date()
            year, month, day = map(int, config["end_date_tc"])
            enddate_TC = datetime(year, month, day).date()
        
        TC_Count = len(pd.date_range(
            start=startdate, end=enddate_TC, freq="W-Fri"
        ))
        self.RF_Cls = np.zeros(TC_Count)
        self.MainTrib_Cls = np.zeros(TC_Count)
        self.Palmer_Cls = np.zeros(TC_Count)
        self.NetInflow_Cls = np.zeros(TC_Count)
        self.Max_RF_MainTrib = np.zeros(TC_Count)
        self.Max_Palmer_NetInf = np.zeros(TC_Count)

        monthly_date_range = pd.date_range(start=startdate, end=enddate, freq="MS")
        if monthly_date_range.empty:
            monthly_date_range  = pd.DatetimeIndex([pd.to_datetime(startdate).replace(day=1)])
        LONINO_Count = len(monthly_date_range)
        self.Seas = np.zeros(LONINO_Count)
        self.M_Seas = np.zeros(LONINO_Count)
        if forecast:
            self.LONINO_Seas_cls = np.zeros(LONINO_Count)
            self.LONINO_M_Seas_cls = np.zeros(LONINO_Count)
        else:
            self.LONINO_Seas_cls = np.zeros(config["month_n"])
            self.LONINO_M_Seas_cls = np.zeros(config["month_n"])

        daily_date_range = pd.date_range(start=startdate, end=enddate, freq="D")
        Seas_Count = len(daily_date_range)
        self.Daily_Seasons = np.zeros(Seas_Count)
        self.Mon = np.zeros(Seas_Count)
        date_range_length = len(daily_date_range)

        n_rows = date_range_length + 1
        self.Lake_Stage = np.zeros(n_rows, dtype=object)
        self.WSM_Zone = np.zeros(n_rows, dtype=object)
        self.Max_Supply = np.zeros(n_rows, dtype=object)
        self.LOSA_Supply = np.zeros(n_rows, dtype=object)
        self.Cut_back = np.zeros(n_rows, dtype=object)
        self.Dem_N_Sup = np.zeros(n_rows, dtype=object)

        self.V10per = np.zeros(date_range_length)
        self.V20per = np.zeros(date_range_length)
        self.V25per = np.zeros(date_range_length)
        self.V30per = np.zeros(date_range_length)
        self.V40per = np.zeros(date_range_length)
        self.V45per = np.zeros(date_range_length)
        self.V50per = np.zeros(date_range_length)
        self.V60per = np.zeros(date_range_length)
        self.NI_Supply = np.zeros(n_rows, dtype=object)
        self.Zone_Code = np.zeros(n_rows, dtype=object)
        self.LO_Zone = np.zeros(n_rows, dtype=object)
        
        self.Zone_D_Trib = np.zeros(date_range_length)
        self.Zone_D_stage = np.zeros(date_range_length)
        self.Zone_D_Seas = np.zeros(date_range_length)
        self.Zone_D_MSeas = np.zeros(date_range_length)
        self.Zone_D_Branch_Code = np.zeros(date_range_length)
        self.Zone_D_Rel_Code = np.zeros(date_range_length)
        self.Zone_C_Trib = np.zeros(date_range_length)
        self.Zone_C_Seas = np.zeros(date_range_length)
        self.Zone_C_MSeas = np.zeros(date_range_length)
        self.Zone_C_MetFcast = np.zeros(date_range_length)
        self.Zone_C_Branch_Code = np.zeros(date_range_length)
        self.Zone_C_Rel_Code = np.zeros(date_range_length)
        self.Zone_B_Trib = np.zeros(date_range_length)
        self.Zone_B_Stage = np.zeros(date_range_length)
        self.Zone_B_Seas = np.zeros(date_range_length)
        self.Zone_B_Branch_Code = np.zeros(date_range_length)
        self.Zone_B_Rel_Code = np.zeros(date_range_length)
        self.DecTree_Relslevel = np.zeros(n_rows, dtype=object)
        self.DayFlags = np.zeros(n_rows, dtype=object)
        self.PlsDay = np.zeros(n_rows, dtype=object)
        self.Release_Level = np.zeros(n_rows, dtype=object)
        self.dh_7days = np.zeros(n_rows, dtype=object)
        self.ZoneCodeminus1Code = np.zeros(n_rows, dtype=object)
        self.ZoneCodeCode = np.zeros(n_rows, dtype=object)
        self.Fraction_of_Zone_height = np.zeros(n_rows, dtype=object)
        self.ReLevelCode_1 = np.zeros(n_rows, dtype=object)
        self.ReLevelCode_2 = np.zeros(n_rows, dtype=object)
        self.ReLevelCode_3_S80 = np.zeros(n_rows, dtype=object)
        self.Outlet2DS_Mult = np.zeros(n_rows, dtype=object)
        self.Outlet2DS_Mult_2 = np.zeros(n_rows, dtype=object)
        self.Outlet2DSRS = np.zeros(n_rows, dtype=object)
        self.Outlet2USRG1 = np.zeros(n_rows, dtype=object)
        self.Sum_Outlet2USRG1 = np.zeros(n_rows, dtype=object)
        self.Outlet2DSBS = np.zeros(n_rows, dtype=object)
        self.Outlet2USBK = np.zeros(n_rows, dtype=object)
        self.ROeast = np.zeros(n_rows, dtype=object)
        self.Outlet2USBS = np.zeros(n_rows, dtype=object)
        self.Sum_Outlet2USBK = np.zeros(n_rows, dtype=object)
        self.Outlet2USRG_Code = np.zeros(n_rows, dtype=object)
        self.Outlet2USRG = np.zeros(n_rows, dtype=object)
        self.Outlet2DS = np.zeros(n_rows, dtype=object)
        self.ReLevelCode_3_S77 = np.zeros(n_rows, dtype=object)
        self.Outlet1US_Mult = np.zeros(n_rows, dtype=object)
        self.Outlet1US_Mult_2 = np.zeros(n_rows, dtype=object)
        self.Outlet1USBSAP = np.zeros(n_rows, dtype=object)
        self.Outlet1USRS = np.zeros(n_rows, dtype=object)
        self.Sum_Outlet1USRS = np.zeros(n_rows, dtype=object)
        self.Outlet1USBK = np.zeros(n_rows, dtype=object)
        self.ROwest = np.zeros(n_rows, dtype=object)
        self.Outlet1DSBS = np.zeros(n_rows, dtype=object)
        self.Outlet1USBS = np.zeros(n_rows, dtype=object)
        self.Outlet1USEWS = np.zeros(n_rows, dtype=object)
        self.Outlet1USREG = np.zeros(n_rows, dtype=object)
        self.Outlet1DS = np.zeros(n_rows, dtype=object)
        self.TotRegEW = np.zeros(n_rows, dtype=object)
        self.Choose_WCA = np.zeros(n_rows, dtype=object)
        self.RegWCA = np.zeros(n_rows, dtype=object)
        self.Choose_L8C51 = np.zeros(n_rows, dtype=object)
        self.RegL8C51 = np.zeros(n_rows, dtype=object)
        self.TotRegSo = np.zeros(n_rows, dtype=object)
        self.Stage2ar = np.zeros(n_rows, dtype=object)
        self.Stage2marsh = np.zeros(n_rows, dtype=object)
        self.RF = np.zeros(n_rows, dtype=object)
        self.ET = np.zeros(n_rows, dtype=object)
        self.Choose_WSA_1 = np.zeros(n_rows, dtype=object)
        self.Choose_WSA_2 = np.zeros(n_rows, dtype=object)
        self.WSA_MIA = np.zeros(n_rows, dtype=object)
        self.WSA_NNR = np.zeros(n_rows, dtype=object)
        self.DSto = np.zeros(n_rows, dtype=object)
        self.Storage = np.zeros(n_rows, dtype=object)

        #Count AP is the same as date_range_length, so replaced it to that
        self.THC_Class_normal_or_above = np.zeros(date_range_length)
        self.Lake_O_Stage_AP = np.zeros(date_range_length)
        self.Lake_O_Schedule_Zone = np.zeros(date_range_length)
        self.LStgCorres = np.zeros(date_range_length)
        self.LowChance_Check = np.zeros(date_range_length)
        self.Outlet1USRS_AP = np.zeros(date_range_length)
        self.Outlet1USBS_AP = np.zeros(date_range_length)
        self.Outlet1USRS_Pre_AP_S77_Baseflow = np.zeros(date_range_length)
        self.Forecast_D_Sal = np.zeros(date_range_length)
        self.n30d_mavg = np.zeros(date_range_length)
        self.n30davgForecast = np.zeros(date_range_length)
        self.LORS08_bf_rel = np.zeros(date_range_length)
        self.LDS_LC6_1 = np.zeros(date_range_length)
        self.S_O = np.zeros(date_range_length)
        self.All_4 = np.zeros(date_range_length)
        self.Sabf = np.zeros(date_range_length)
        self.Swbf = np.zeros(date_range_length)
        self.Swbu = np.zeros(date_range_length)
        self.All_4andStage = np.zeros(date_range_length)
        self.All_4andStagein = np.zeros(date_range_length)
        self.P_AP_BF_Stg = np.zeros(date_range_length)
        self.Logic_test_1 = np.zeros(date_range_length)
        self.Post_Ap_Baseflow = np.zeros(date_range_length)
        self.Outlet1USRSplusPreAPS77bsf = np.zeros(date_range_length)
        self.AndEstNeedsLakeWater = np.zeros(date_range_length)
        self.Choose_PAPEWS_1 = np.zeros(date_range_length)
        self.Choose_PAPEWS_2 = np.zeros(date_range_length)
        self.AndLowChance61stagelessth11 = np.zeros(date_range_length)
        self.ATHCnora = np.zeros(date_range_length)
        self.Post_AP_EWS = np.zeros(date_range_length)
        self.Post_AP_Baseflow_EWS_cfs = np.zeros(date_range_length)
