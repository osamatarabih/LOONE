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
        # Load the configuration
        config = load_config(working_path)
        self.data_dir = working_path

        # Read the data from the configuration
        self.read_data(config)

    def read_data(self, config: dict) -> None:
        """
        Reads all the data from the configuration files.

        Args:
            config (dict): The configuration dictionary.
        """
        self._read_sfwmm_daily_outputs(config)
        self._read_water_shortage_management_break_points(config)
        self._read_weekly_demand_input_file(config)
        self._read_weekly_tributary_condition_data(config)
        self._read_seasonal_lonino_data(config)
        self._read_netinflow_data(config)
        self._read_actual_water_demand(config)
        self._read_rainfall_volume_data(config)
        self._read_et_volume_data(config)
        self._read_c44_runoff_data(config)
        self._read_c43_runoff_data(config)
        self._read_sum_basin_runoffs(config)
        self._read_mean_monthly_basin_runoffs(config)
        self._read_s80_and_s77_regulatory_release_rates(config)
        self._read_ce_sle_turns(config)
        self._read_pulses_data(config)
        self._read_target_stage_may_1st(config)
        self._read_target_stage_june_1st(config)
        self._read_estuary_needs_water(config)
        self._read_eaa_mia_runoff(config)
        self._read_storage_deviation(config)
        self._read_calibration_parameters(config)

    def _read_sfwmm_daily_outputs(self, config: dict) -> None:
        """
        Reads the SFWMM daily output file.

        This file contains the daily output of SFWMM, including the daily
        water demand, rainfall, evapotranspiration, water supply, and the
        water levels in the lake and estuary.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'sfwmm_daily_outputs'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the SFWMM daily outputs
        file_path = os.path.join(self.data_dir, config["sfwmm_daily_outputs"])

        # Read the CSV file into a DataFrame
        self.SFWMM_Daily_Outputs = pd.read_csv(file_path)

    def _read_water_shortage_management_break_points(self, config: dict) -> None:
        """
        Reads the Water Shortage Management and Regulation Schedule Break points csv file.

        This method loads the break points data from a CSV file specified in the
        configuration and assigns it to the WSMs_RSBKs attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'wsms_rsbps' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the Water Shortage Management and Regulation Schedule Break points
        file_path = os.path.join(self.data_dir, config["wsms_rsbps"])

        # Read the CSV file into a DataFrame
        self.WSMs_RSBKs = pd.read_csv(file_path)

    def _read_weekly_demand_input_file(self, config: dict) -> None:
        """
        Reads the weekly demand input file.

        This method loads the weekly demand data from a CSV file specified in the
        configuration and assigns it to the Weekly_dmd attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'losa_wkly_dmd' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the weekly demand
        file_path = os.path.join(self.data_dir, config["losa_wkly_dmd"])

        # Read the CSV file into a DataFrame
        self.Weekly_dmd = pd.read_csv(file_path)

    def _read_weekly_tributary_condition_data(self, config: dict) -> None:
        """Reads weekly tributary condition (NetRF, S65E runoff, Palmer index, and
        Netinflow) data for the study period.

        This method loads the weekly tributary condition data from a CSV file
        specified in the configuration and assigns it to the Wkly_Trib_Cond
        attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'trib_cond_wkly_data' CSV
                           file.

        Returns:
            None
        """
        # Construct the full file path for the tributary condition data
        file_path = os.path.join(self.data_dir, config["trib_cond_wkly_data"])

        # Read the CSV file into a DataFrame
        self.Wkly_Trib_Cond = pd.read_csv(file_path)

    def _read_seasonal_lonino_data(self, config: dict) -> None:
        """
        Reads Seasonal LONINO data for LORS simulations.

        This method loads the latest values of Seasonal and Multi Seasonal
        LONINO data from CSV files specified in the configuration and assigns
        it to the LONINO_Seas_data and LONINO_Mult_Seas_data attributes.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'seasonal_lonino', and
                           'multi_seasonal_lonino' CSV files.

        Returns:
            None
        """
        # Construct the full file path for the seasonal LONINO data CSV
        seasonal_file_path = os.path.join(self.data_dir, config["seasonal_lonino"])
        multi_seasonal_file_path = os.path.join(self.data_dir, config["multi_seasonal_lonino"])

        # Read the CSV file into a DataFrame
        self.LONINO_Seas_data = pd.read_csv(seasonal_file_path)
        self.LONINO_Mult_Seas_data = pd.read_csv(multi_seasonal_file_path)

    def _read_netinflow_data(self, config: dict) -> None:
        """
        Reads Netinflow data (ac-ft/day) from a CSV file.

        This method loads the Netinflow data from a CSV file specified in the
        configuration and assigns it to the NetInf_Input attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'netflows_acft' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the Netinflow data CSV
        file_path = os.path.join(self.data_dir, config["netflows_acft"])

        # Read the CSV file into a DataFrame
        self.NetInf_Input = pd.read_csv(file_path)

    def _read_actual_water_demand(self, config: dict) -> None:
        """
        Reads actual daily water demand (output of SFWMM) from a CSV file.

        This method loads the actual daily water demand from a CSV file specified
        in the configuration and assigns it to the SFWMM_W_dmd attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'water_dmd' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the actual daily water demand data CSV
        file_path = os.path.join(self.data_dir, config["water_dmd"])

        # Read the CSV file into a DataFrame
        self.SFWMM_W_dmd = pd.read_csv(file_path)

    def _read_rainfall_volume_data(self, config: dict) -> None:
        """
        Reads Rainfall Volume data from a CSV file.

        This method loads the Rainfall Volume data from a CSV file specified in
        the configuration and assigns it to the RF_Vol attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'rf_vol' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the Rainfall Volume data CSV
        file_path = os.path.join(self.data_dir, config["rf_vol"])

        # Read the CSV file into a DataFrame
        self.RF_Vol = pd.read_csv(file_path)

    def _read_et_volume_data(self, config: dict) -> None:
        """
        Reads ET Volume data from a CSV file.

        This method loads the ET Volume data from a CSV file specified in the
        configuration and assigns it to the ET_Vol attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'et_vol' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the ET Volume data CSV
        file_path = os.path.join(self.data_dir, config["et_vol"])

        # Read the CSV file into a DataFrame
        self.ET_Vol = pd.read_csv(file_path)

    def _read_c44_runoff_data(self, config: dict) -> None:
        """
        Reads the C44 Runoff data from a CSV file.

        This method loads the C44 Runoff data from a CSV file specified in the
        configuration and assigns it to the C44RO attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'c44ro' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the C44 Runoff data CSV
        file_path = os.path.join(self.data_dir, config["c44ro"])

        # Read the CSV file into a DataFrame
        self.C44_Runoff = pd.read_csv(file_path)

    def _read_c43_runoff_data(self, config: dict) -> None:
        """
        Reads the daily C43 Runoff data from a CSV file.

        This method loads the daily C43 Runoff data from a CSV file specified in
        the configuration and assigns it to the C43RO_Daily attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'c43ro' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the daily C43 Runoff data CSV
        file_path = os.path.join(self.data_dir, config["c43ro"])

        # Read the CSV file into a DataFrame
        self.C43RO_Daily = pd.read_csv(file_path)

    def _read_sum_basin_runoffs(self, config: dict) -> None:
        """
        Reads the mean monthly sum Basin runoffs file (i.e., sum of flow volume
        for each month of the entire period).

        This file contains the sum of monthly flow volumes for each of the
        simulated basins in the SFWMM model.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'basin_ro_inputs'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the sum Basin runoffs CSV
        file_path = os.path.join(self.data_dir, config["basin_ro_inputs"])

        # Read the CSV file into a DataFrame
        self.Sum_Basin_RO = pd.read_csv(file_path)

    def _read_mean_monthly_basin_runoffs(self, config: dict) -> None:
        """Reads mean monthly basin runoffs (cfs).

        This method loads the mean monthly basin runoffs from three CSV files
        specified in the configuration and assigns them to the C43RO,
        C44RO, and SLTRIB attributes.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the paths to the 'c43ro_monthly',
                           'c44ro_nonthly', and 'sltrib_monthly' CSV files.

        Returns:
            None
        """
        # Construct the full file paths for the mean monthly basin runoffs CSVs
        c43ro_file_path = os.path.join(self.data_dir, config["c43ro_monthly"])
        c44ro_file_path = os.path.join(self.data_dir, config["c44ro_nonthly"])
        sltrib_file_path = os.path.join(self.data_dir, config["sltrib_monthly"])

        # Read the CSV files into DataFrames
        self.C43RO = pd.read_csv(c43ro_file_path)
        self.C44RO = pd.read_csv(c44ro_file_path)
        self.SLTRIB = pd.read_csv(sltrib_file_path)

    def _read_s80_and_s77_regulatory_release_rates(self, config: dict) -> None:
        """Reads S80 and S77 Regulatory release rates for LORS2008 (Inputs from
        ActiveSchedule Tab in LOOPS Model!).

        This method loads the S80 and S77 Regulatory release rates from CSV
        files specified in the configuration and assigns them to the
        S80_RegRelRates and S77_RegRelRates attributes.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 's80_regulatory_release_rates'
                           and 's77_regulatory_release_rates' CSV files.

        Returns:
            None
        """
        # Construct the full file path for the S80 and S77 Regulatory release rates
        s80_file_path = os.path.join(self.data_dir, config["s80_regulatory_release_rates"])
        s77_file_path = os.path.join(self.data_dir, config["s77_regulatory_release_rates"])

        # Read the CSV files into a DataFrame
        self.S80_RegRelRates = pd.read_csv(s80_file_path)
        self.S77_RegRelRates = pd.read_csv(s77_file_path)

    def _read_ce_sle_turns(self, config: dict) -> None:
        """FIXME: Reads CE and SLE turns out of the LOOPS model as inputs
        (to figure out how to calculate it later in this model!).

        This method loads the CE and SLE turns from a CSV file
        specified in the configuration and assigns them to the
        CE_SLE_turns attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'ce_sle_turns_inputs'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the CE and SLE turns
        file_path = os.path.join(self.data_dir, config["ce_sle_turns_inputs"])

        # Read the CSV file into a DataFrame
        self.CE_SLE_turns = pd.read_csv(file_path)

    def _read_pulses_data(self, config: dict) -> None:
        """
        Reads the pulses input data from a CSV file.

        This method loads the pulses input data from a CSV file specified in
        the configuration and assigns it to the Pulses attribute.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'pulses_inputs' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the pulses data CSV
        file_path = os.path.join(self.data_dir, config["pulses_inputs"])

        # Read the CSV file into a DataFrame
        self.Pulses = pd.read_csv(file_path)

    def _read_target_stage_may_1st(self, config: dict) -> None:
        """
        Reads the target stage for May 1st from a CSV file specified in the
        configuration. The data contains lake stages below 11ft, and the CSV
        file includes a 'Date' column which is parsed as datetime objects.

        This data is used in the LOOPS Spreadsheet Model in the 'Tab Pulses'
        sheet.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'may_1st_lake_stage_below_11ft'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the May 1st target stage CSV
        file_path = os.path.join(self.data_dir, config["may_1st_lake_stage_below_11ft"])

        # Read the CSV file into a DataFrame, parsing the 'Date' column
        self.Targ_Stg_May_1st = pd.read_csv(file_path, parse_dates=["Date"])

    def _read_target_stage_june_1st(self, config: dict) -> None:
        """
        Reads the target stage for June 1st from a CSV file specified in the
        configuration. The data contains lake stages below 11ft, and the CSV
        file includes a 'Date' column which is parsed as datetime objects.

        This data is used in the LOOPS Spreadsheet Model in the 'Tab Pulses'
        sheet.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'june_1st_lake_stage_below_11ft'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the June 1st target stage CSV
        file_path = os.path.join(self.data_dir, config["june_1st_lake_stage_below_11ft"])

        # Read the CSV file into a DataFrame, parsing the 'Date' column
        self.Targ_Stg_June_1st = pd.read_csv(file_path, parse_dates=["Date"])

    def _read_estuary_needs_water(self, config: dict) -> None:
        """
        Reads the Estuary needs water (to be calculated later in this model!).

        The Estuary needs water is a calculated value that determines whether the
        Estuary needs water from Lake Okeechobee. The calculation is done later in
        this model, and the result is stored in a DataFrame.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'estuary_needs_water_input'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the Estuary needs water CSV
        file_path = os.path.join(
            self.data_dir,
            config["estuary_needs_water_input"],
        )

        # Read the CSV file into a DataFrame, parsing the 'Date' column
        self.Estuary_needs_water = pd.read_csv(file_path, parse_dates=["Date"])

    def _read_eaa_mia_runoff(self, config: dict) -> None:
        """
        Reads EAA_MIA_RUNOFF in cfs.

        This data is used in the LOOPS Spreadsheet Model in the 'Tab Pulses'
        sheet.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'eaa_mia_ro_inputs'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the EAA MIA RUNOFF CSV
        file_path = os.path.join(
            self.data_dir,
            config["eaa_mia_ro_inputs"],
        )

        # Read the CSV file into a DataFrame, parsing the 'Date' column
        self.EAA_MIA_RUNOFF = pd.read_csv(file_path, parse_dates=["Date"])

    def _read_storage_deviation(self, config: dict) -> None:
        """
        Reads the calculated storage deviation daily values from a CSV file.

        This method loads the storage deviation data, which is essential for
        understanding the variations in storage capacity over time.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'storage_deviation' CSV file.

        Returns:
            None
        """
        # Construct the full file path for the storage deviation CSV
        file_path = os.path.join(self.data_dir, config["storage_deviation"])

        # Read the CSV file into a DataFrame
        self.Storage_dev_df = pd.read_csv(file_path)

    def _read_calibration_parameters(self, config: dict) -> None:
        """
        Reads calibration parameters.

        This method loads the calibration parameters, which are used to
        adjust the model parameters to fit the observed data.

        Args:
            config (dict): A dictionary containing the configuration paths,
                           specifically the path to the 'calibration_parameters'
                           CSV file.

        Returns:
            None
        """
        # Construct the full file path for the calibration parameters CSV
        file_path = os.path.join(self.data_dir, config["calibration_parameters"])

        # Read the CSV file into a DataFrame
        self.Cal_Par = pd.read_csv(file_path)
