import os
import pandas as pd
from loone.utils import load_config


class Data:
    def __init__(self, config_path: str):
        self.config = load_config(config_path)
        self.data_dir = self.config["working_path"]

        # Read SFWMM Daily Output File
        self.SFWMM_Daily_Outputs = pd.read_csv(
            os.path.join(self.data_dir, f"SFWMM_Daily_Outputs.csv")
        )
        # read Water Shortage Management and Regulation Schedule Break points csv
        self.WSMs_RSBKs = pd.read_csv(
            os.path.join(self.data_dir, f"WSMs_RSBPs.csv")
        )
        # Read the weekly demand input file
        self.Weekly_dmd = pd.read_csv(
            os.path.join(self.data_dir, f"LOSA_wkly_dmd.csv")
        )
        # read Weekly tributary condition (NetRF, S65E runoff, Palmer index, and
        # Netinflow) data for the study period
        self.Wkly_Trib_Cond = pd.read_csv(
            os.path.join(self.data_dir, f"Trib_cond_wkly_data.csv")
        )
        # Read Seasonal LONINO data (the latest values for use with LORS
        # simulations) defined as Option 4 in the LOOPS Model (TC-LONINO Tab).
        self.LONINO_Seas_data = pd.read_csv(
            os.path.join(self.data_dir, f"Seasonal_LONINO.csv")
        )
        self.LONINO_Mult_Seas_data = pd.read_csv(
            os.path.join(self.data_dir, f"Multi_Seasonal_LONINO.csv")
        )
        # read Netinflow data(ac-ft/day)
        self.NetInf_Input = pd.read_csv(
            os.path.join(self.data_dir, f"Netflows_acft.csv")
        )
        # read actual daily water demand (output of SFWMM)
        self.SFWMM_W_dmd = pd.read_csv(
            os.path.join(self.data_dir, f"Water_dmd.csv")
        )
        # read Rainfall Volume data
        self.RF_Vol = pd.read_csv(os.path.join(self.data_dir, f"RFVol.csv"))
        # read ET Vol data
        self.ET_Vol = pd.read_csv(os.path.join(self.data_dir, f"ETVol.csv"))
        # Read the C44 Runoff data which is output of SFWMM simulation.
        self.C44_Runoff = pd.read_csv(
            os.path.join(self.data_dir, f"C44RO.csv")
        )

        # Read the mean monthly sum Basin runoffs file (i.e., sum of flow volume
        # for each month of the entire period)
        self.Sum_Basin_RO = pd.read_csv(
            os.path.join(self.data_dir, f"Basin_RO_inputs.csv")
        )
        # Read Mean Monthly basin runoffs (cfs)
        self.C43RO = pd.read_csv(
            os.path.join(self.data_dir, f"C43RO_Monthly.csv")
        )
        self.C44RO = pd.read_csv(
            os.path.join(self.data_dir, f"C44RO_Monthly.csv")
        )
        self.SLTRIB = pd.read_csv(
            os.path.join(self.data_dir, f"SLTRIB_Monthly.csv")
        )
        # Read S80 and S77 Regulatory release rates for LORS2008 (Inputs from
        # ActiveSchedule Tab in LOOPS Model!)
        self.S77_RegRelRates = pd.read_csv(
            os.path.join(self.data_dir, f"S77_RegRelRates.csv")
        )
        self.S80_RegRelRates = pd.read_csv(
            os.path.join(self.data_dir, f"S80_RegRelRates.csv")
        )
        # Read the daily C43Runoff
        self.C43RO_Daily = pd.read_csv(
            os.path.join(self.data_dir, f"C43RO.csv")
        )
        # FIXME: I will use the CE and SLE turns out of the LOOPS model as inputs
        # here (to figure out how to calculate it later in this model!)
        self.CE_SLE_turns = pd.read_csv(
            os.path.join(self.data_dir, f"CE_SLE_turns_inputs.csv")
        )
        # Read the pulses input data (Tab Pulses in the LOOPS Spreadsheet Model!)
        self.Pulses = pd.read_csv(
            os.path.join(self.data_dir, f"Pulses_Inputs.csv")
        )
        # Note that I entered values of (-9999) in days of the year that do not
        # include data (i.e. June 2nd to September 30th)
        self.Targ_Stg_June_1st = pd.read_csv(
            os.path.join(
                self.data_dir,
                f"Chance of June 1st Lake stage falling below 11.0ft" ".csv",
            ),
            parse_dates=["Date"],
        )
        self.Targ_Stg_May_1st = pd.read_csv(
            os.path.join(
                self.data_dir,
                f"Chance of May 1st Lake stage falling below 11.0ft" ".csv",
            ),
            parse_dates=["Date"],
        )
        # FIXME
        # The Estuary needs water needs to be updated if needed!
        # Read the "Estuary needs water" for now (to be calculated later!)
        self.Estuary_needs_water = pd.read_csv(
            os.path.join(self.data_dir, f"Estuary_needs_water_Input.csv")
        )
        # Read EAA_MIA_RUNOFF in cfs
        self.EAA_MIA_RUNOFF = pd.read_csv(
            os.path.join(self.data_dir, f"EAA_MIA_RUNOFF_Inputs.csv")
        )
        # Read calculated Storage deviation daily values
        self.Stroage_dev_df = pd.read_csv(
            os.path.join(self.data_dir, f"Storage_Dev.csv")
        )
