import os
import pandas as pd
from loone.utils import load_config


class Data:
    """A class that represents the data used for running LOONE."""

    def __init__(self, working_path: str):
        """
        Initializes the Data object.

        Args:
            working_path (str): The working directory path.
        """
        config = load_config(working_path)
        self.data_dir = working_path

        # Read SFWMM Daily Output File
        self.SFWMM_Daily_Outputs = pd.read_csv(
            os.path.join(self.data_dir, config["sfwmm_daily_outputs"])
        )
        # read Water Shortage Management and Regulation Schedule Break points csv
        self.WSMs_RSBKs = pd.read_csv(
            os.path.join(self.data_dir, config["wsms_rsbps"])
        )
        # Read the weekly demand input file
        self.Weekly_dmd = pd.read_csv(
            os.path.join(self.data_dir, config["losa_wkly_dmd"])
        )
        # read Weekly tributary condition (NetRF, S65E runoff, Palmer index, and
        # Netinflow) data for the study period
        self.Wkly_Trib_Cond = pd.read_csv(
            os.path.join(self.data_dir, config["trib_cond_wkly_data"])
        )
        # Read Seasonal LONINO data (the latest values for use with LORS
        # simulations) defined as Option 4 in the LOOPS Model (TC-LONINO Tab).
        self.LONINO_Seas_data = pd.read_csv(
            os.path.join(self.data_dir, config["seasonal_lonino"])
        )
        self.LONINO_Mult_Seas_data = pd.read_csv(
            os.path.join(self.data_dir, config["multi_seasonal_lonino"])
        )
        # read Netinflow data(ac-ft/day)
        self.NetInf_Input = pd.read_csv(
            os.path.join(self.data_dir, config["netflows_acft"])
        )
        # read actual daily water demand (output of SFWMM)
        self.SFWMM_W_dmd = pd.read_csv(
            os.path.join(self.data_dir, config["water_dmd"])
        )
        # read Rainfall Volume data
        self.RF_Vol = pd.read_csv(
            os.path.join(self.data_dir, config["rf_vol"])
        )
        # read ET Vol data
        self.ET_Vol = pd.read_csv(
            os.path.join(self.data_dir, config["et_vol"])
        )
        # Read the C44 Runoff data which is output of SFWMM simulation.
        self.C44_Runoff = pd.read_csv(
            os.path.join(self.data_dir, config["c44ro"])
        )
        # Read the daily C43Runoff
        self.C43RO_Daily = pd.read_csv(
            os.path.join(self.data_dir, config["c43ro"])
        )
        # Read the mean monthly sum Basin runoffs file (i.e., sum of flow volume
        # for each month of the entire period)
        self.Sum_Basin_RO = pd.read_csv(
            os.path.join(self.data_dir, config["basin_ro_inputs"])
        )
        # Read Mean Monthly basin runoffs (cfs)
        self.C43RO = pd.read_csv(
            os.path.join(self.data_dir, config["c43ro_monthly"])
        )
        self.C44RO = pd.read_csv(
            os.path.join(self.data_dir, config["c44ro_nonthly"])
        )
        self.SLTRIB = pd.read_csv(
            os.path.join(self.data_dir, config["sltrib_monthly"])
        )
        # Read S80 and S77 Regulatory release rates for LORS2008 (Inputs from
        # ActiveSchedule Tab in LOOPS Model!)
        self.S77_RegRelRates = pd.read_csv(
            os.path.join(self.data_dir, config["s77_regulatory_release_rates"])
        )
        self.S80_RegRelRates = pd.read_csv(
            os.path.join(self.data_dir, config["s80_regulatory_release_rates"])
        )
        # FIXME: I will use the CE and SLE turns out of the LOOPS model as inputs
        # here (to figure out how to calculate it later in this model!)
        self.CE_SLE_turns = pd.read_csv(
            os.path.join(self.data_dir, config["ce_sle_turns_inputs"])
        )
        # Read the pulses input data (Tab Pulses in the LOOPS Spreadsheet Model!)
        self.Pulses = pd.read_csv(
            os.path.join(self.data_dir, config["pulses_inputs"])
        )
        # Note that I entered values of (-9999) in days of the year that do not
        # include data (i.e. June 2nd to September 30th)
        self.Targ_Stg_June_1st = pd.read_csv(
            os.path.join(
                self.data_dir,
                config["june_1st_lake_stage_below_11ft"],
            ),
            parse_dates=["Date"],
        )
        self.Targ_Stg_May_1st = pd.read_csv(
            os.path.join(
                self.data_dir,
                config["may_1st_lake_stage_below_11ft"],
            ),
            parse_dates=["Date"],
        )
        # FIXME
        # The Estuary needs water needs to be updated if needed!
        # Read the "Estuary needs water" for now (to be calculated later!)
        self.Estuary_needs_water = pd.read_csv(
            os.path.join(self.data_dir, config["estuary_needs_water_input"])
        )
        # Read EAA_MIA_RUNOFF in cfs
        self.EAA_MIA_RUNOFF = pd.read_csv(
            os.path.join(self.data_dir, config["eaa_mia_ro_inputs"])
        )
        # Read calculated Storage deviation daily values
        self.Stroage_dev_df = pd.read_csv(
            os.path.join(self.data_dir, config["storage_deviation"])
        )
        # Read calibration parameters
        self.Cal_Par = pd.read_csv(
            os.path.join(self.data_dir, config["calibration_parameters"])
        )
