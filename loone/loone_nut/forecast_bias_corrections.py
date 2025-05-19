import argparse
import os
import numpy as np
import pandas as pd


def bias_correct(station: str, forecasted_tp_lake: list, most_recent_tp_obs: float) -> np.ndarray:
    """Bias corrects the forecasted TP concentrations in Lake Okeechobee based on observed data.
    
    Args:
        station (str): The station name of the station where the values were forecasted and observed. One of: 'L001', 'L005', 'L008', 'L004', 'L006', 'L007', 'LZ40'.
        forecasted_tp_lake (list): A list of forecasted TP concentrations in Lake Okeechobee.
        most_recent_tp_obs (float): The most recent observed TP concentration at the given station.
        
    Returns:
        sim_bc_qnt_ent_per (np.ndarray): A list of bias corrected TP concentrations in Lake Okeechobee.
    """
    STATION_PARAMETERS = {
        'L001': {'obs_quantiles': [98, 122, 152], 'correction_factors': [0.5913678934015472, 0.8353793753288427, 1.0524001582697413, 1.4605983806263183]},
        'L005': {'obs_quantiles': [56.25, 85, 115.25], 'correction_factors': [0.29828881937048657, 0.5283271275855149, 0.7632318707941915, 1.188700366068547]},
        'L008': {'obs_quantiles': [98, 140.5, 174], 'correction_factors': [0.6062543848348515, 0.9218176856485235, 1.2001724608466868, 1.7926071099770216]},
        'L004': {'obs_quantiles': [122, 145, 177.5], 'correction_factors': [0.8955760993304536, 1.1510189122660914, 1.3689041532929282, 1.9750371175582473]},
        'L006': {'obs_quantiles': [106.5, 135.5, 174], 'correction_factors': [0.7347463391896812, 1.0422794056109943, 1.3076781699380404, 1.9930823718916406]},
        'L007': {'obs_quantiles': [86, 116, 151.125], 'correction_factors': [0.5927850757579236, 0.9145550431423951, 1.1594227707702824, 1.6784629147274226]},
        'LZ40': {'obs_quantiles': [132.25, 158.75, 183.75], 'correction_factors': [0.9184276576796921, 1.2443912322787118, 1.4516529344806937, 2.171551720078875]}
    }
    
    # Validate the station
    station = station.upper()
    
    if station not in STATION_PARAMETERS:
        raise ValueError(f"Invalid station: {station}  Expected one of: {list(STATION_PARAMETERS.keys())}")
    
    # Get the station parameters
    obs_quantiles = STATION_PARAMETERS[station]['obs_quantiles']
    correction_factors = STATION_PARAMETERS[station]['correction_factors']
    
    sim_bc_qnt_ent_per = np.zeros(len(forecasted_tp_lake))
    
    # Apply the correction factor to Forecasted TP_Lake_(N/S)
    for i in range(len(forecasted_tp_lake)):
        if most_recent_tp_obs <= obs_quantiles[0]:
            sim_bc_qnt_ent_per[i] = forecasted_tp_lake[i] * correction_factors[0] 
        elif most_recent_tp_obs > obs_quantiles[0] and most_recent_tp_obs <= obs_quantiles[1]:
            sim_bc_qnt_ent_per[i] = forecasted_tp_lake[i] * correction_factors[1] 
        elif most_recent_tp_obs > obs_quantiles[1] and most_recent_tp_obs <= obs_quantiles[2]:
            sim_bc_qnt_ent_per[i] = forecasted_tp_lake[i] * correction_factors[2]
        elif most_recent_tp_obs > obs_quantiles[2]:
            sim_bc_qnt_ent_per[i] = forecasted_tp_lake[i] * correction_factors[3]
        else:
            sim_bc_qnt_ent_per[i] = forecasted_tp_lake[i]
    
    return sim_bc_qnt_ent_per


def bias_correct_file(station: str, file_tp_forecast: str, file_tp_observed: str) -> pd.DataFrame:
    """Bias corrects the forecasted TP concentrations in Lake Okeechobee based on observed data.
    
    Args:
        station (str): The station name of the station where the values were forecasted and observed. One of: 'L001', 'L005', 'L008', 'L004', 'L006', 'L007', 'LZ40'.
        file_tp_forecast (str): The file path to the forecasted TP concentrations in Lake Okeechobee.
        file_tp_observed (str): The file path to the observed TP concentrations in Lake Okeechobee.
        
    Returns:
        sim_bc_qnt_ent_per (np.ndarray): A list of bias corrected TP concentrations in Lake Okeechobee.
    """
    # Get the station in expected format
    station = station.upper()
    
    # Read the forecasted and observed data
    df_forecast = pd.read_csv(file_tp_forecast)
    df_observed = pd.read_csv(file_tp_observed)

    # Get whether the station is in the north or south of the lake
    forecast_column_name = ''
    
    if station in ['L001', 'L005', 'L008']:
        forecast_column_name = 'TP_Lake_N'
    elif station in ['L004', 'L006', 'L007', 'LZ40']:
        forecast_column_name = 'TP_Lake_S'
    
    # Get the most recent observed TP concentration
    tp_observed = df_observed.iloc[-1, 2]
    
    # Get the dates of the forecasted TP concentrations
    tp_forecast_dates = df_forecast['Date'].values
    
    # Get the forecasted TP concentrations
    tp_forecast = df_forecast[forecast_column_name].values
    
    # Bias correct the forecasted TP concentrations
    tp_corrected = bias_correct(station, tp_forecast, tp_observed)

    df_corrected = pd.DataFrame({
        'Date': tp_forecast_dates,
        forecast_column_name: tp_corrected
    })

    return df_corrected


def main(input_dir: str, output_dir: str):
    """Bias corrects the forecasted TP data using the files LOONE_Nut_Output_ens<ensemble>.csv and 'water_quality_<station>_PHOSPHATE, TOTAL AS P.csv'.
    
    Args:
        input_dir (str): The directory path where the input files are located.
        output_dir (str): The directory path where the output files will be saved.
    """
    tp_forecast_file = ''
    tp_observed_file = ''
    
    # Bias correct each ensemble using each station's observed data
    for station in ['L001', 'L005', 'L008', 'L004', 'L006', 'L007', 'LZ40']:
        for i in range(1, 52):
            # Read the forecasted and observed data
            tp_forecast_file = os.path.join(input_dir, f'LOONE_Nut_Output_ens{i:02d}.csv')
            tp_observed_file = os.path.join(input_dir, f'water_quality_{station}_PHOSPHATE, TOTAL AS P.csv')
            
            # Bias correct the forecasted TP concentrations
            df = bias_correct_file(station, tp_forecast_file, tp_observed_file)
            
            # Save the bias corrected data
            df.to_csv(os.path.join(output_dir, f'LOONE_Nut_Output_ens{i:02d}_{station}_corrected.csv'), index=False)


if __name__ == '__main__':
    # Parse the command line arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument('input_dir', nargs=1, help='The path to the directory containing the input files.')
    argparser.add_argument('output_dir', nargs=1, help='The path to the directory where the output files will be saved.')
    
    args = argparser.parse_args()
    
    input_dir = args.input_dir[0]
    output_dir = args.output_dir[0]
    
    # Run the main function
    main(input_dir, output_dir)