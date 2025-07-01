import os
import pandas as pd
from loone.utils import load_config


class Data:
    """A class that represents the data used for running LOONE."""

    def __init__(self, working_path: str, forecast: bool = False, ensemble: int = None) -> None:
        """
        Initializes the Data object.

        Args:
            working_path (str): The working directory path.
        """
        # Load the configuration
        config = load_config(working_path)
        self.data_dir = working_path

        # Read the data from the configuration
        self.read_data(config, forecast, ensemble)

    def read_data(self, config: dict, forecast, ensemble) -> None:
        """
        Reads all the data from the configuration files.

        Args:
            config (dict): The configuration dictionary.
            forecast (bool): Whether to use forecast data or not.
            ensemble (int): The ensemble number for the forecast.
        """
        self.WSMs_RSBKs = self._read_csv(config, "wsms_rsbps")
        self.Weekly_dmd = self._read_csv(config, "losa_wkly_dmd")
        if forecast:
            self.SFWMM_Daily_Outputs = pd.read_csv(
                os.path.join(self.data_dir, f"SFWMM_Daily_Outputs_forecast.csv")
            )
            self.Wkly_Trib_Cond = pd.read_csv(os.path.join(self.data_dir, f"Trib_cond_predicted_{ensemble:02d}.csv"))
            self.LONINO_Seas_data = pd.read_csv(
                os.path.join(self.data_dir, f"Seasonal_LONINO_forecast.csv")
            )
            self.LONINO_Mult_Seas_data = pd.read_csv(
                os.path.join(self.data_dir, f"Multi_Seasonal_LONINO_forecast.csv")
            )
            self.NetInf_Input = pd.read_csv(os.path.join(self.data_dir, f"Netflows_acft_geoglows_{ensemble:02d}.csv"))
            self.Sum_Basin_RO = pd.read_csv(os.path.join(self.data_dir, f"Basin_RO_inputs_{ensemble:02d}.csv"))
            self.C44_Runoff = pd.read_csv(os.path.join(self.data_dir, f"C44RO_{ensemble:02d}.csv"))
            self.C43RO_Daily = pd.read_csv(os.path.join(self.data_dir, f"C43RO_{ensemble:02d}.csv"))
            self.C43RO = pd.read_csv(os.path.join(self.data_dir, f"C43RO_Monthly_{ensemble:02d}.csv"))
            self.C44RO = pd.read_csv(os.path.join(self.data_dir, f"C44RO_Monthly_{ensemble:02d}.csv"))
            self.Estuary_needs_water = pd.read_csv(
                os.path.join(self.data_dir, "Estuary_needs_water_Input_forecast.csv")
            )
            self.EAA_MIA_RUNOFF = pd.read_csv(
                os.path.join(self.data_dir, "EAA_MIA_RUNOFF_Inputs_forecast.csv")
            )
            self.SFWMM_W_dmd = pd.read_csv(
                os.path.join(self.data_dir, "Water_dmd_forecast.csv")
            )
        else:
            self.SFWMM_Daily_Outputs = self._read_csv(
                config, "sfwmm_daily_outputs"
            )
            self.Wkly_Trib_Cond = self._read_csv(config, "trib_cond_wkly_data")
            self.LONINO_Seas_data = self._read_csv(config, "seasonal_lonino")
            self.LONINO_Mult_Seas_data = self._read_csv(
                config, "multi_seasonal_lonino"
            )
            self.NetInf_Input = self._read_csv(config, "netflows_acft")
            self.Sum_Basin_RO = self._read_csv(config, "basin_ro_inputs")
            self.C44_Runoff = self._read_csv(config, "c44ro")
            self.C43RO_Daily = self._read_csv(config, "c43ro")
            self.C43RO = self._read_csv(config, "c43ro_monthly")
            self.C44RO = self._read_csv(config, "c44ro_nonthly")
            self.SFWMM_W_dmd = self._read_csv(config, "water_dmd")
            self.Estuary_needs_water = self._read_csv(
                config, "estuary_needs_water_input"
            )
            self.EAA_MIA_RUNOFF = self._read_csv(config, "eaa_mia_ro_inputs")
        self.RF_Vol = self._read_csv(config, "rf_vol")
        self.ET_Vol = self._read_csv(config, "et_vol")
        self.SLTRIB = self._read_csv(config, "sltrib_monthly")
        self.S80_RegRelRates = self._read_csv(
            config, "s80_regulatory_release_rates"
        )
        self.S77_RegRelRates = self._read_csv(
            config, "s77_regulatory_release_rates"
        )
        self.CE_SLE_turns = self._read_csv(config, "ce_sle_turns_inputs")
        self.Pulses = self._read_csv(config, "pulses_inputs")
        self.Targ_Stg_May_1st = self._read_csv(
            config, "may_1st_lake_stage_below_11ft", parse_dates=["Date"]
        )
        self.Targ_Stg_June_1st = self._read_csv(
            config, "june_1st_lake_stage_below_11ft", parse_dates=["Date"]
        )
        self.Storage_dev_df = self._read_csv(config, "storage_deviation")
        self.Cal_Par = self._read_csv(config, "calibration_parameters")

    def _read_csv(self, config: dict, key: str, **kwargs) -> pd.DataFrame:
        """
        Helper function to construct the file path and read the CSV file.

        Args:
            config (dict): The configuration dictionary.
            key (str): The key to look up the file path in the configuration.
            **kwargs: Additional arguments to pass to pd.read_csv.

        Returns:
            pd.DataFrame: The loaded DataFrame.
        """
        file_path = os.path.join(self.data_dir, config[key])
        return pd.read_csv(file_path, **kwargs)
