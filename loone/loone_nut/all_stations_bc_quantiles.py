import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import argparse


def bias_correct_file(file_path_observed: str,
                      file_path_simulated: str,
                      station: str,
                      forecast_mode: bool = False) -> pd.DataFrame:
    """Bias correct the simulated total phosphorus values against a station in Lake Okeechobee.

    Args:
        file_path_observed: str: The path to the file containing the observed total phosphorus data
        file_path_simulated: str: The path to the file containing the simulated total phosphorus data
        station: str: The name of the station where the total phosphorus values were observed. One of: 'L001', 'L005', 'L008', 'L004', 'L006', 'L007', 'LZ40'
        forecast_mode (bool, optional): Whether to run the model in forecast mode. Defaults to False.

    Returns:
        pd.DataFrame: A dataframe containing the bias corrected total phosphorus values for the station.
    """
    # Station Regions
    # NOTE: Station L008 is considered to be in both the north and south regions of the lake
    north_stations = ['L001', 'L005', 'L008']
    south_stations = ['L004', 'L006', 'L007', 'L008', 'LZ40']

    # Validate the station
    station = station.upper()

    if station not in north_stations and station not in south_stations:
        raise ValueError(f"Station {station} is not a valid station. Please choose from the following stations: {north_stations + south_stations}")

    # Check that input files exist
    if not os.path.exists(file_path_observed):
        raise FileNotFoundError(f"File {file_path_observed} does not exist.")

    if not os.path.exists(file_path_simulated):
        raise FileNotFoundError(f"File {file_path_simulated} does not exist.")

    # Get the region of the station
    station_region = 'TP_Lake_N' if station in north_stations else 'TP_Lake_S'

    # Read simulated TP concentrations in Lake Okeechobee
    df_simulated_total_phosphorus = pd.read_csv(file_path_simulated)
    df_simulated_total_phosphorus['date'] = pd.to_datetime(df_simulated_total_phosphorus['Date'])

    # Read observations at given station
    df_station_observations = pd.read_csv(file_path_observed)
    df_station_observations['date'] = pd.to_datetime(df_station_observations['date'])

    # Station is located in both the north and south regions of the lake
    if station in north_stations and station in south_stations:
        # Run bias correction for both regions
        simulated_bias_corrected_quantile_entry_percentile_north = bias_correct_df(df_simulated_total_phosphorus,
                                                                                   df_station_observations,
                                                                                   station,
                                                                                   'TP_Lake_N',
                                                                                   forecast_mode)
        simulated_bias_corrected_quantile_entry_percentile_south = bias_correct_df(df_simulated_total_phosphorus,
                                                                                   df_station_observations,
                                                                                   station,
                                                                                   'TP_Lake_S',
                                                                                   forecast_mode)

        df_results = pd.merge(simulated_bias_corrected_quantile_entry_percentile_north,
                              simulated_bias_corrected_quantile_entry_percentile_south,
                              how='outer',
                              on='date')

        # Return the bias corrected values for both regions
        return df_results
    else:
        # Run bias correction for the region where the station is located
        simulated_bias_corrected_quantile_entry_percentile = bias_correct_df(df_simulated_total_phosphorus,
                                                                             df_station_observations,
                                                                             station,
                                                                             station_region,
                                                                             forecast_mode)

        # Return the bias corrected values
        return simulated_bias_corrected_quantile_entry_percentile


def bias_correct_df(df_simulated_total_phosphorus: pd.DataFrame,
                    df_station_observations: pd.DataFrame,
                    station: str,
                    station_region: str,
                    forecast_mode: bool = False) -> pd.DataFrame:
    """Bias correct the simulated total phosphorus values against a station in Lake Okeechobee.

    Args:
        df_simulated_total_phosphorus: pd.DataFrame: A dataframe containing the simulated total phosphorus values
        df_station_observations: pd.DataFrame: A dataframe containing the observed total phosphorus values
        station: str: The name of the station where the total phosphorus values were observed. One of: 'L001', 'L005', 'L008', 'L004', 'L006', 'L007', 'LZ40'
        station_region: str: The region of the lake where the station is located. One of: 'TP_Lake_N', 'TP_Lake_S'
        forecast_mode (bool, optional): Whether to run the model in forecast mode. Defaults to False.

    Returns:
        np.array: A numpy array containing the bias corrected total phosphorus values for the station.
    """
    # Create a dataframe for the TP simulations
    df_simulated_total_phosphorus_region = pd.DataFrame()
    df_simulated_total_phosphorus_region['date'] = df_simulated_total_phosphorus['date']
    df_simulated_total_phosphorus_region[station_region] = df_simulated_total_phosphorus[station_region]
    df_simulated_total_phosphorus_region = df_simulated_total_phosphorus_region.set_index(['date'])
    df_simulated_total_phosphorus_region.index = pd.to_datetime(df_simulated_total_phosphorus_region.index, unit='ns')
    df_simulated_total_phosphorus_monthly = df_simulated_total_phosphorus_region.resample('ME').mean()

    # Create a dataframe for the observed TP values
    df_station_observations_processed = pd.DataFrame()
    df_station_observations_processed['date'] = df_station_observations['date']
    df_station_observations_processed[f'{station}_Obs_mg/m3'] = df_station_observations[f'{station}_PHOSPHATE, TOTAL AS P_mg/L'] * 1000  # mg/m3
    df_station_observations_processed = df_station_observations_processed.set_index(['date'])
    df_station_observations_processed.index = pd.to_datetime(df_station_observations_processed.index, unit='ns')
    df_station_observations_monthly = df_station_observations_processed.resample('ME').mean()

    df_bias_correct = pd.merge(df_simulated_total_phosphorus_monthly, df_station_observations_monthly, how='left', on='date')
    df_bias_correct_no_nan = df_bias_correct.dropna()
    observed_values = df_bias_correct_no_nan[f'{station}_Obs_mg/m3'].values
    simulated_values = df_bias_correct_no_nan[station_region].values

    # Determine R2 for simulations related to observations at the station
    if simulated_values.size == 0 or observed_values.size == 0:
        print(f'No data available for bias correction calculation for station {station}_{station_region[-1]}.')
        return df_bias_correct_no_nan
    slope, intercept, r_value, p_value, std_err = stats.linregress(simulated_values, observed_values)

    # Calculate R-squared (coefficient of determination)
    r_squared = r_value**2
    print(f'R2_{station}_{station_region[-1]}:', r_squared)

    # Calculate the quantiles of the observed and simulated time series.
    obs_quantiles = np.percentile(observed_values, [25, 50, 75])

    # Calculate the average of values less than the 25th quantile
    average_below_25th = np.mean([value for value in observed_values if value <= obs_quantiles[0]])
    # Calculate the average of values less than the 50th quantile (median)
    average_below_50th = np.mean([value for value in observed_values if obs_quantiles[0] < value <= obs_quantiles[1]])
    # Calculate the average of values less than the 75th quantile
    average_below_75th = np.mean([value for value in observed_values if obs_quantiles[1] < value <= obs_quantiles[2]])
    # Calculate the average of values more than the 75th quantile
    average_above_75th = np.mean([value for value in observed_values if value > obs_quantiles[2]])

    # Average of Sim
    average_tp = simulated_values.mean()

    # Calculate the correction factor
    correction_factor_1 = average_below_25th / average_tp
    correction_factor_2 = average_below_50th / average_tp
    correction_factor_3 = average_below_75th / average_tp
    correction_factor_4 = average_above_75th / average_tp

    # Apply the correction factor
    sim_BC_Qnt_EntPer = np.zeros(len(observed_values))

    for i in range(len(observed_values)):
        if observed_values[i] <= obs_quantiles[0]:
            sim_BC_Qnt_EntPer[i] = simulated_values[i] * correction_factor_1
        elif observed_values[i] > obs_quantiles[0] and observed_values[i] <= obs_quantiles[1]:
            sim_BC_Qnt_EntPer[i] = simulated_values[i] * correction_factor_2
        elif observed_values[i] > obs_quantiles[1] and observed_values[i] <= obs_quantiles[2]:
            sim_BC_Qnt_EntPer[i] = simulated_values[i] * correction_factor_3
        elif observed_values[i] > obs_quantiles[2]:
            sim_BC_Qnt_EntPer[i] = simulated_values[i] * correction_factor_4
        else:
            sim_BC_Qnt_EntPer[i] = simulated_values[i]

    # Calculate R-squared (coefficient of determination)
    if sim_BC_Qnt_EntPer.size == 0 or observed_values.size == 0:
        print(f'No data available for bias correction calculation for station {station}_{station_region[-1]}.')
        return df_bias_correct_no_nan
    slope, intercept, r_value, p_value, std_err = stats.linregress(sim_BC_Qnt_EntPer, observed_values)
    r_squared = r_value**2
    print(f'R2_{station}_{station_region[-1]}_BC_Qnt:', r_squared)

    # Check that length of bias corrected array matches the length of the date column in the dataframe
    assert len(sim_BC_Qnt_EntPer) == len(df_bias_correct_no_nan.index), 'The length of the bias corrected array does not match the length of the date column in the dataframe.'

    # Return the bias corrected values as a dataframe
    return pd.DataFrame(sim_BC_Qnt_EntPer, index=df_bias_correct_no_nan.index, columns=[station_region])


def main(input_dir: str,
         output_dir: str,
         simulated_file: str = 'LOONE_Nut_Output.csv',
         forecast_mode: bool = False) -> None:
    """Bias correct the simulated total phosphorus values against each station in Lake Okeechobee.

    Expects the input files to be in the format 'water_quality_{station}_PHOSPHATE, TOTAL AS P.csv' and 'LOONE_Nut_Output.csv'.
    This function expects the input files to be named 'water_quality_{station}_PHOSPHATE, TOTAL AS P.csv' for the observed data
    and 'LOONE_Nut_Output.csv' for the simulated data. It applies bias correction to the simulated total phosphorus values. The corrected
    values are saved to the specified output directory under the file name 'LOONE_Nut_Output_{station}_corrected.csv'.

    Args:
        input_dir: str: The path to the directory containing the input files.
        output_dir: str: The path to the directory where the output files will be saved.
        simulated_file: str: The name of the file containing the simulated Nut Ouput data. Defaults to 'LOONE_Nut_Output.csv'.
        forecast_mode (bool, optional): Whether to run the model in forecast mode. Defaults to False.

    Returns:
        None
    """

    # Bias correct the simulated total phosphorus values against each station
    for station in ['L001', 'L005', 'L008', 'L004', 'L006', 'L007', 'LZ40']:
        # Get the file paths to the observed and simulated total phosphorus data
        file_path_observed = os.path.join(input_dir, f'water_quality_{station}_PHOSPHATE, TOTAL AS P.csv')
        file_path_simulated = os.path.join(input_dir, simulated_file)

        # Bias correct the simulated total phosphorus values against the station
        df_results = bias_correct_file(file_path_observed, file_path_simulated, station, forecast_mode)

        # Save the bias corrected values to a file
        df_results.to_csv(os.path.join(output_dir, f'LOONE_Nut_Output_{station}_corrected.csv'))


if __name__ == '__main__':
    # Parse the command line arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        'input_dir',
        help='The path to the directory containing the input files.',
        nargs=1
    )
    argparser.add_argument(
        'output_dir',
        help='The path to the directory where the output files will be saved.',
        nargs=1
    )
    argparser.add_argument(
        '--simulated_file',
        help='The name of the file containing the simulated Nut Ouput data.',
        default='LOONE_Nut_Output.csv',
        nargs=1
    )
    argparser.add_argument(
        "--forecast_mode",
        action="store_true",
        help="Flag to indicate that the bias correction is running in forecast mode.",
    )

    args = argparser.parse_args()

    input_dir = args.input_dir[0]
    output_dir = args.output_dir[0]
    simulated_file = args.simulated_file[0]
    forecast_mode = args.forecast_mode

    # Run the main function
    main(input_dir, output_dir, simulated_file, forecast_mode)
