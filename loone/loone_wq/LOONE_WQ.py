# Daily Nitrate-Nitrite NOx and Chlorophyll-a Modeling

# import required packages
import os
import argparse
import pandas as pd
import numpy as np
from scipy.integrate import odeint
from loone.utils import loone_nchla_fns, load_config
from loone.data import Data as DClass


def _calculate_summer_month_loads(constit_loads_m: pd.DataFrame) -> tuple:
    """
    Calculate the summer month loads for NOx and Chla.

    Args:
        constit_loads_m (pd.DataFrame): The constituent loads DataFrame.

    Returns:
        tuple: A tuple containing arrays of summer month loads for NOx and Chla (StL and Cal).
    """
    Smr_Mnth_NOx_StL = []
    Smr_Mnth_NOx_Cal = []
    Smr_Mnth_Chla_StL = []
    Smr_Mnth_Chla_Cal = []
    for i in range(len(constit_loads_m.index)):
        if constit_loads_m['date'].iloc[i].month in [5, 6, 7, 8, 9, 10]:
            Smr_Mnth_NOx_StL.append(constit_loads_m['NO_Load_StL'].iloc[i])
            Smr_Mnth_Chla_StL.append(constit_loads_m['Chla_Load_StL'].iloc[i])
            Smr_Mnth_NOx_Cal.append(constit_loads_m['NO_Load_Cal'].iloc[i])
            Smr_Mnth_Chla_Cal.append(constit_loads_m['Chla_Load_Cal'].iloc[i])

    Smr_Mnth_NOx_StL_arr = np.asarray(Smr_Mnth_NOx_StL)
    Smr_Mnth_Chla_StL_arr = np.asarray(Smr_Mnth_Chla_StL)
    Smr_Mnth_NOx_Cal_arr = np.asarray(Smr_Mnth_NOx_Cal)
    Smr_Mnth_Chla_Cal_arr = np.asarray(Smr_Mnth_Chla_Cal)

    return Smr_Mnth_NOx_StL_arr, Smr_Mnth_Chla_StL_arr, Smr_Mnth_NOx_Cal_arr, Smr_Mnth_Chla_Cal_arr


def _prepare_constit_loads_df(inflows: pd.DataFrame, NO_Load_Cal: pd.Series, NO_Load_StL: pd.Series,
                              NO_Load_South: pd.Series, Chla_Load_Cal: pd.Series, Chla_Load_StL: pd.Series,
                              Chla_Load_South: pd.Series) -> tuple:
    """
    Prepare the constituent loads DataFrame and its monthly resampled version.

    Args:
        inflows (pd.DataFrame): The inflows DataFrame.
        NO_Load_Cal (pd.Series): The NO Load Cal series.
        NO_Load_StL (pd.Series): The NO Load StL series.
        NO_Load_South (pd.Series): The NO Load South series.
        Chla_Load_Cal (pd.Series): The Chla Load Cal series.
        Chla_Load_StL (pd.Series): The Chla Load StL series.
        Chla_Load_South (pd.Series): The Chla Load South series.

    Returns:
        tuple: A tuple containing the constituent loads DataFrame and its monthly resampled version.
    """
    constit_loads_df = pd.DataFrame(inflows['date'], columns=['date'])
    constit_loads_df['date'] = pd.to_datetime(constit_loads_df['date'])

    constit_loads_df['NO_Load_Cal'] = pd.to_numeric(NO_Load_Cal) / 1E9  # tons
    constit_loads_df['NO_Load_StL'] = pd.to_numeric(NO_Load_StL) / 1E9  # tons
    constit_loads_df['NO_Load_South'] = pd.to_numeric(NO_Load_South) / 1E9  # tons
    constit_loads_df['Chla_Load_Cal'] = pd.to_numeric(Chla_Load_Cal) / 1E6  # Kgs
    constit_loads_df['Chla_Load_StL'] = pd.to_numeric(Chla_Load_StL) / 1E6  # Kgs
    constit_loads_df['Chla_Load_South'] = pd.to_numeric(Chla_Load_South) / 1E6  # Kgs

    constit_loads_df = constit_loads_df.set_index('date')
    constit_loads_df.index = pd.to_datetime(constit_loads_df.index, unit='ns')

    constit_loads_m = constit_loads_df.resample('ME').sum()
    constit_loads_df = constit_loads_df.reset_index()
    constit_loads_m = constit_loads_m.reset_index()

    return constit_loads_df, constit_loads_m


def _update_nitro_model_output(nitro_model_output: pd.DataFrame, NO_N: pd.Series, NO_S: pd.Series, NO_MEAN: pd.Series, Sim_Chla: pd.Series, Sim_Chla_N: pd.Series, Sim_Chla_S: pd.Series, temperature: pd.Series) -> tuple:
    """
    Update the nitro model output DataFrame with the given series and resample it monthly.

    Args:
        nitro_model_output (pd.DataFrame): The nitro model output DataFrame.
        NO_N (pd.Series): The NO_N series.
        NO_S (pd.Series): The NO_S series.
        NO_MEAN (pd.Series): The NO_MEAN series.
        Sim_Chla (pd.Series): The Sim_Chla series.
        Sim_Chla_N (pd.Series): The Sim_Chla_N series.
        Sim_Chla_S (pd.Series): The Sim_Chla_S series.
        temperature (pd.Series): The temperature series.

    Returns:
        tuple: A tuple containing the updated nitro model output DataFrame and its monthly resampled version.
    """
    nitro_model_output['NO_N'] = pd.to_numeric(NO_N)
    nitro_model_output['NO_S'] = pd.to_numeric(NO_S)
    nitro_model_output['NO_M'] = pd.to_numeric(NO_MEAN)
    nitro_model_output['Sim_Chla'] = pd.to_numeric(Sim_Chla)
    nitro_model_output['Sim_Chla_N'] = pd.to_numeric(Sim_Chla_N)
    nitro_model_output['Sim_Chla_S'] = pd.to_numeric(Sim_Chla_S)
    nitro_model_output['Water Temp'] = pd.to_numeric(temperature)  # C

    nitro_model_output = nitro_model_output.set_index('date')
    nitro_model_output.index = pd.to_datetime(nitro_model_output.index, unit='ns')

    nitro_mod_out_m = nitro_model_output.resample('ME').mean()
    nitro_model_output = nitro_model_output.reset_index()
    nitro_mod_out_m = nitro_mod_out_m.reset_index()

    return nitro_model_output, nitro_mod_out_m


def _load_data(workspace: str, flow_path: str, forecast_mode: bool, photo_period_filename: str, config: dict, ensemble_member: int = None) -> dict:
    """
    Load the necessary data files from the workspace.

    Args:
        workspace (str): The path to the workspace directory.
        flow_path (str): The path to the flow data file.
        forecast_mode (bool): Whether to use forecast mode.
        photo_period_filename (str): The filename for the photo period data.
        config (dict): The configuration dictionary that holds the contents of config.yaml.
        ensemble_member (int): The ensemble member number.

    Returns:
        dict: A dictionary containing the loaded data.
    """
    data = {}
    data['inflows'] = pd.read_csv(os.path.join(workspace, flow_path))

    if forecast_mode:
        data['temperature_data'] = pd.read_csv(os.path.join(workspace, 'Filled_WaterT_predicted.csv'))
        #TODO: Predict this
        data['dissolved_oxygen'] = pd.read_csv(os.path.join(workspace, 'LO_DO_Clean_daily.csv'))
        data['radiation_data'] = pd.read_csv(os.path.join(workspace, 'LO_RADT_data.csv'))
        data['storage_data'] = pd.read_csv(os.path.join(workspace, f'Average_LO_Storage_3MLag_{ensemble_member:02d}.csv'))
        data['chlorophyll_a_north_data'] = pd.read_csv(os.path.join(workspace, 'N_Merged_Chla.csv'))  # microgram/L
        data['chlorophyll_a_south_data'] = pd.read_csv(os.path.join(workspace, 'S_Merged_Chla.csv'))  # microgram/L
        data['external_nitrate_loadings'] = pd.read_csv(os.path.join(workspace, 'LO_External_Loadings_NO.csv'))
    else:
        data['temperature_data'] = pd.read_csv(os.path.join(workspace, 'Filled_WaterT.csv'))
        data['dissolved_oxygen'] = pd.read_csv(os.path.join(workspace, 'LO_DO_Clean_daily.csv'))
        data['radiation_data'] = pd.read_csv(os.path.join(workspace, 'LO_RADT_data.csv'))
        data['storage_data'] = pd.read_csv(os.path.join(workspace, config['sto_stage']))
        data['chlorophyll_a_north_data'] = pd.read_csv(os.path.join(workspace, 'N_Merged_Chla.csv'))  # microgram/L
        data['chlorophyll_a_south_data'] = pd.read_csv(os.path.join(workspace, 'S_Merged_Chla.csv'))  # microgram/L
        data['external_nitrate_loadings'] = pd.read_csv(os.path.join(workspace, 'LO_External_Loadings_NO.csv'))  # mg

    s65e_basename = 'water_quality_S65E_NITRATE+NITRITE-N_Interpolated_forecast.csv' if forecast_mode else 'water_quality_S65E_NITRATE+NITRITE-N_Interpolated.csv'
    data['s65e_nitrate_data'] = pd.read_csv(os.path.join(workspace, s65e_basename))  # mg/m3

    data['chlorophyll_a_loads_in'] = pd.read_csv(os.path.join(workspace, 'Chla_Loads_In.csv'))  # mg

    s65e_chlorophyll_a_basename = 'S65E_Chla_Merged_forecast.csv' if forecast_mode else 'S65E_Chla_Merged.csv'
    data['s65e_chlorophyll_a_data'] = pd.read_csv(os.path.join(workspace, s65e_chlorophyll_a_basename))  # mg/m3

    data['photo_period'] = pd.read_csv(os.path.join(workspace, f'{photo_period_filename}.csv'))
    data['lo_orthophosphate_north_data'] = pd.read_csv(os.path.join(workspace, 'N_OP.csv'))  # mg/m3
    data['lo_orthophosphate_south_data'] = pd.read_csv(os.path.join(workspace, 'S_OP.csv'))  # mg/m3
    data['lo_dissolved_inorganic_nitrogen_north_data'] = pd.read_csv(os.path.join(workspace, 'N_DIN.csv'))  # mg/m3
    data['lo_dissolved_inorganic_nitrogen_south_data'] = pd.read_csv(os.path.join(workspace, 'S_DIN.csv'))  # mg/m3

    return data


def _calculate_observed_values(lo_dissolved_inorganic_nitrogen_north_data: pd.DataFrame,
                               lo_dissolved_inorganic_nitrogen_south_data: pd.DataFrame,
                               chlorophyll_a_north_data: pd.DataFrame,
                               chlorophyll_a_south_data: pd.DataFrame,
                               lo_orthophosphate_north_data: pd.DataFrame,
                               lo_orthophosphate_south_data: pd.DataFrame) -> dict:
    """
    Calculate the observed values for various parameters.

    Args:
        lo_dissolved_inorganic_nitrogen_north_data (pd.DataFrame): The north dissolved inorganic nitrogen data.
        lo_dissolved_inorganic_nitrogen_south_data (pd.DataFrame): The south dissolved inorganic nitrogen data.
        chlorophyll_a_north_data (pd.DataFrame): The north chlorophyll-a data.
        chlorophyll_a_south_data (pd.DataFrame): The south chlorophyll-a data.
        lo_orthophosphate_north_data (pd.DataFrame): The north orthophosphate data.
        lo_orthophosphate_south_data (pd.DataFrame): The south orthophosphate data.

    Returns:
        dict: A dictionary containing the calculated observed values.
    """
    observed_values = {}

    observed_values['ammonium_north'] = lo_dissolved_inorganic_nitrogen_north_data['NH4'].astype(float)
    observed_values['ammonium_south'] = lo_dissolved_inorganic_nitrogen_south_data['NH4'].astype(float)
    observed_values['ammonium'] = (observed_values['ammonium_north'] + observed_values['ammonium_south']) / 2

    observed_values['nitrogen_oxide_north'] = lo_dissolved_inorganic_nitrogen_north_data['NO'].astype(float)
    observed_values['nitrogen_oxide_south'] = lo_dissolved_inorganic_nitrogen_south_data['NO'].astype(float)
    observed_values['nitrogen_oxide'] = (observed_values['nitrogen_oxide_north'] + observed_values['nitrogen_oxide_south']) / 2

    observed_values['chlorophyll_a_north'] = chlorophyll_a_north_data['Chla'].astype(float)
    observed_values['chlorophyll_a_south'] = chlorophyll_a_south_data['Chla'].astype(float)
    observed_values['chlorophyll_a'] = (observed_values['chlorophyll_a_north'] + observed_values['chlorophyll_a_south']) / 2

    observed_values['dissolved_inorganic_nitrogen_north'] = lo_dissolved_inorganic_nitrogen_north_data['DIN'].astype(float)
    observed_values['dissolved_inorganic_nitrogen_south'] = lo_dissolved_inorganic_nitrogen_south_data['DIN'].astype(float)
    observed_values['dissolved_inorganic_nitrogen'] = (observed_values['dissolved_inorganic_nitrogen_north'] + observed_values['dissolved_inorganic_nitrogen_south']) / 2

    observed_values['dissolved_inorganic_phosphorus_north'] = lo_orthophosphate_north_data['OP'].astype(float)
    observed_values['dissolved_inorganic_phosphorus_south'] = lo_orthophosphate_south_data['OP'].astype(float)
    observed_values['dissolved_inorganic_phosphorus'] = (observed_values['dissolved_inorganic_phosphorus_north'] + observed_values['dissolved_inorganic_phosphorus_south']) / 2

    return observed_values


def _get_light_attenuation(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Get light attenuation parameters.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the light attenuation parameters.
    """
    light_attenuation = {
        "Kw_NO": calibration_parameter_no[4],
        "Kw_Chla": calibration_parameter_chla[4],
        "Kc_NO": calibration_parameter_no[5],
        "Kc_Chla": calibration_parameter_chla[5],
    }
    return light_attenuation


def _get_light_limitation(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Get light limitation and inhibition coefficients.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the light limitation and inhibition coefficients.
    """
    light_limitation = {
        "K1_NO": calibration_parameter_no[6],
        "K1_Chla": calibration_parameter_chla[6],
        "K2_NO": calibration_parameter_no[7],
        "K2_Chla": calibration_parameter_chla[7],
    }
    return light_limitation


def _get_phytoplankton_growth(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Get maximum phytoplankton growth rate.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the maximum phytoplankton growth rate.
    """
    phytoplankton_growth = {
        "G_max_NO": calibration_parameter_no[8],
        "G_max_Chla": calibration_parameter_chla[8],
    }
    return phytoplankton_growth


def _get_half_saturation_coefficients(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Get half saturation coefficients.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the half saturation coefficients.
    """
    half_saturation = {
        "K_NH_NO": calibration_parameter_no[9],
        "K_NH_Chla": calibration_parameter_chla[9],
        "K_Nitr_NO": calibration_parameter_no[10],
        "K_Nitr_Chla": calibration_parameter_chla[10],
        "K_DIN_NO": calibration_parameter_no[9] + calibration_parameter_no[10],
        "K_DIN_Chla": calibration_parameter_chla[9] + calibration_parameter_chla[10],
        "KP_NO": calibration_parameter_no[11],
        "KP_Chla": calibration_parameter_chla[11],
        "K_TN_NO": calibration_parameter_no[12],
        "K_TN_Chla": calibration_parameter_chla[12],
        "YNOChla_NO": calibration_parameter_no[13],
        "YNOChla_Chla": calibration_parameter_chla[13],
    }
    return half_saturation


def _get_temperatures(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Get temperature parameters.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the temperature parameters.
    """
    temperatures = {
        "T_opt_NO": calibration_parameter_no[14],
        "T_min_NO": calibration_parameter_no[15],
        "T_max_NO": calibration_parameter_no[16],
        "T_opt_Chla": calibration_parameter_chla[14],
        "T_min_Chla": calibration_parameter_chla[15],
        "T_max_Chla": calibration_parameter_chla[16],
    }
    return temperatures


def _calculate_dissolved_oxygen(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Calculate dissolved oxygen parameters.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the dissolved oxygen parameters.
    """
    dissolved_oxygen = {
        "KDO_NO": calibration_parameter_no[17],
        "KDO_Chla": calibration_parameter_chla[17],
    }
    return dissolved_oxygen


def _calculate_sediment_release(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Calculate sediment release parameters.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the sediment release parameters.
    """
    sediment_release = {
        "S_NO_NO": calibration_parameter_no[18],
        "S_NO_Chla": calibration_parameter_chla[18],
    }
    return sediment_release


def _get_nitrification_parameters(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Get nitrification parameters for the detailed method.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the nitrification parameters.
    """
    nitrification_parameters = {
        "kni0_no": calibration_parameter_no[0],
        "kni0_chla": calibration_parameter_chla[0],
        "kni_no": calibration_parameter_no[1],
        "kni_chla": calibration_parameter_chla[1],
    }
    return nitrification_parameters


def _get_denitrification_parameters(calibration_parameter_no: list, calibration_parameter_chla: list) -> dict:
    """
    Get denitrification parameters for the detailed method.

    Args:
        calibration_parameter_no (list): The calibration parameters for NO.
        calibration_parameter_chla (list): The calibration parameters for Chla.

    Returns:
        dict: A dictionary containing the denitrification parameters.
    """
    denitrification_parameters = {
        "kden0_no": calibration_parameter_no[2],
        "kden0_chla": calibration_parameter_chla[2],
        "kden_no": calibration_parameter_no[3],
        "kden_chla": calibration_parameter_chla[3],
    }
    return denitrification_parameters


def _calculate_nutrient_limitation_factors(dissolved_inorganic_phosphorus_observed: float,
                                           dissolved_inorganic_nitrogen_observed: float,
                                           dissolved_inorganic_phosphorus_north_observed: float,
                                           dissolved_inorganic_phosphorus_south_observed: float,
                                           dissolved_inorganic_nitrogen_north_observed: float,
                                           dissolved_inorganic_nitrogen_south_observed: float,
                                           KP_NO: float, KP_Chla: float, K_DIN_NO: float, K_DIN_Chla: float) -> dict:
    """
    Calculate the nutrient limitation factors.

    Args:
        dissolved_inorganic_phosphorus_observed (float): The observed dissolved inorganic phosphorus.
        dissolved_inorganic_nitrogen_observed (float): The observed dissolved inorganic nitrogen.
        dissolved_inorganic_phosphorus_north_observed (float): The observed dissolved inorganic phosphorus in the north.
        dissolved_inorganic_phosphorus_south_observed (float): The observed dissolved inorganic phosphorus in the south.
        dissolved_inorganic_nitrogen_north_observed (float): The observed dissolved inorganic nitrogen in the north.
        dissolved_inorganic_nitrogen_south_observed (float): The observed dissolved inorganic nitrogen in the south.
        KP_NO (float): The KP_NO value.
        KP_Chla (float): The KP_Chla value.
        K_DIN_NO (float): The K_DIN_NO value.
        K_DIN_Chla (float): The K_DIN_Chla value.

    Returns:
        dict: A dictionary containing the nutrient limitation factors.
    """
    nutrient_limitation_factors = {
        "f_P_NO": dissolved_inorganic_phosphorus_observed / (KP_NO + dissolved_inorganic_phosphorus_observed),
        "f_P_Chla": dissolved_inorganic_phosphorus_observed / (KP_Chla + dissolved_inorganic_phosphorus_observed),
        "f_N_NO": dissolved_inorganic_nitrogen_observed / (K_DIN_NO + dissolved_inorganic_nitrogen_observed),
        "f_N_Chla": dissolved_inorganic_nitrogen_observed / (K_DIN_Chla + dissolved_inorganic_nitrogen_observed),
        "f_P_N_NO": dissolved_inorganic_phosphorus_north_observed / (KP_NO + dissolved_inorganic_phosphorus_north_observed),
        "f_P_N_Chla": dissolved_inorganic_phosphorus_north_observed / (KP_Chla + dissolved_inorganic_phosphorus_north_observed),
        "f_P_S_NO": dissolved_inorganic_phosphorus_south_observed / (KP_NO + dissolved_inorganic_phosphorus_south_observed),
        "f_P_S_Chla": dissolved_inorganic_phosphorus_south_observed / (KP_Chla + dissolved_inorganic_phosphorus_south_observed),
        "f_N_N_NO": dissolved_inorganic_nitrogen_north_observed / (K_DIN_NO + dissolved_inorganic_nitrogen_north_observed),
        "f_N_N_Chla": dissolved_inorganic_nitrogen_north_observed / (K_DIN_Chla + dissolved_inorganic_nitrogen_north_observed),
        "f_N_S_NO": dissolved_inorganic_nitrogen_south_observed / (K_DIN_NO + dissolved_inorganic_nitrogen_south_observed),
        "f_N_S_Chla": dissolved_inorganic_nitrogen_south_observed / (K_DIN_Chla + dissolved_inorganic_nitrogen_south_observed),
    }

    return nutrient_limitation_factors


def _calculate_inflows_and_outflows(i: int, storage_dev: list, q_i: list, q_o: list, external_nitrite_nitrate: list,
                                    s65e_nitrite_nitrate: list, external_chlorophyll_a: list,
                                    s65e_chlorophyll_a: list) -> tuple:
    """
    Calculate the inflows and outflows for the given index.

    Args:
        i (int): The current index.
        storage_dev (list): The storage deviation list.
        q_i (list): The inflow list.
        q_o (list): The outflow list.
        external_nitrite_nitrate (list): The external nitrite nitrate list.
        s65e_nitrite_nitrate (list): The S65E nitrite nitrate list.
        external_chlorophyll_a (list): The external chlorophyll-a list.
        s65e_chlorophyll_a (list): The S65E chlorophyll-a list.

    Returns:
        tuple: A tuple containing the calculated inflows and outflows.
    """
    if storage_dev[i] >= 0:
        Q_I_M = q_i[i] + storage_dev[i] * 1233.48  # m3/d
        Q_O_M = q_o[i]
        External_NO_M = external_nitrite_nitrate[i] + Q_I_M * s65e_nitrite_nitrate[i]
        External_Chla_M = external_chlorophyll_a[i] + Q_I_M * s65e_chlorophyll_a[i]
    else:
        Q_O_M = q_o[i] - storage_dev[i] * 1233.48  # m3/d
        Q_I_M = q_i[i]
        External_NO_M = external_nitrite_nitrate[i]
        External_Chla_M = external_chlorophyll_a[i]

    return Q_I_M, Q_O_M, External_NO_M, External_Chla_M


def _calculate_nitrification_rates(i: int, loone_nchla_fns: object, nit_denit_opt: str, total_nitrification_rate: float,
                                   kni0_no: float, kni_no: float, ammonium_north_observed: list,
                                   ammonium_south_observed: list, K_NH_NO: float, dissolved_oxygen: list, KDO_NO: float,
                                   fox_min: float, DO_Cr_n: float, DO_Opt_n: float, a: float, theta_ni: float,
                                   temperature: list) -> tuple:
    """
    Calculate the nitrification rates for the given index.

    Args:
        i (int): The current index.
        loone_nchla_fns (object): The LOONE NChla functions object.
        nit_denit_opt (str): The nitrification-denitrification option.
        total_nitrification_rate (float): The total nitrification rate.
        kni0_no (float): The kni0_no value.
        kni_no (float): The kni_no value.
        ammonium_north_observed (list): The observed ammonium in the north.
        ammonium_south_observed (list): The observed ammonium in the south.
        K_NH_NO (float): The K_NH_NO value.
        dissolved_oxygen (list): The dissolved oxygen list.
        KDO_NO (float): The KDO_NO value.
        fox_min (float): The fox_min value.
        DO_Cr_n (float): The DO_Cr_n value.
        DO_Opt_n (float): The DO_Opt_n value.
        a (float): The a value.
        theta_ni (float): The theta_ni value.
        temperature (list): The temperature list.

    Returns:
        tuple: A tuple containing the nitrification rates for the north and south regions.
    """
    Nit_R_N_NO = loone_nchla_fns.Nit_Rate('%s' % nit_denit_opt, total_nitrification_rate, kni0_no, kni_no,
                                          ammonium_north_observed[i], K_NH_NO, dissolved_oxygen[i], KDO_NO,
                                          fox_min, DO_Cr_n, DO_Opt_n, a)
    Nit_R_N_T_NO = Nit_R_N_NO * theta_ni ** (temperature[i] - 20)
    Nit_R_S_NO = loone_nchla_fns.Nit_Rate('%s' % nit_denit_opt, total_nitrification_rate, kni0_no, kni_no,
                                          ammonium_south_observed[i], K_NH_NO, dissolved_oxygen[i], KDO_NO,
                                          fox_min, DO_Cr_n, DO_Opt_n, a)
    Nit_R_S_T_NO = Nit_R_S_NO * theta_ni ** (temperature[i] - 20)

    return Nit_R_N_T_NO, Nit_R_S_T_NO


def _calculate_denitrification_rates(i: int, loone_nchla_fns: object, nit_denit_opt: str,
                                     total_denitrification_rate: float, kden0_no: float, kden_no: float, NO_N: list,
                                     NO_S: list, K_Nitr_NO: float, dissolved_oxygen: list, KDO_NO: float,
                                     DO_Cr_d: float, DO_Opt_d: float, theta_deni: float, temperature: list) -> tuple:
    """
    Calculate the denitrification rates for the given index.

    Args:
        i (int): The current index.
        loone_nchla_fns (object): The LOONE NChla functions object.
        nit_denit_opt (str): The nitrification-denitrification option.
        total_denitrification_rate (float): The total denitrification rate.
        kden0_no (float): The kden0_no value.
        kden_no (float): The kden_no value.
        NO_N (list): The NO in the north.
        NO_S (list): The NO in the south.
        K_Nitr_NO (float): The K_Nitr_NO value.
        dissolved_oxygen (list): The dissolved oxygen list.
        KDO_NO (float): The KDO_NO value.
        DO_Cr_d (float): The DO_Cr_d value.
        DO_Opt_d (float): The DO_Opt_d value.
        theta_deni (float): The theta_deni value.
        temperature (list): The temperature list.

    Returns:
        tuple: A tuple containing the denitrification rates for the north and south regions.
    """
    Denit_R_N_NO = loone_nchla_fns.Denit_Rate('%s' % nit_denit_opt, total_denitrification_rate, kden0_no,
                                              kden_no, NO_N[i], K_Nitr_NO, dissolved_oxygen[i], KDO_NO, DO_Cr_d,
                                              DO_Opt_d)
    Denit_R_N_T_NO = Denit_R_N_NO * theta_deni ** (temperature[i] - 20)
    Denit_R_S_NO = loone_nchla_fns.Denit_Rate('%s' % nit_denit_opt, total_denitrification_rate, kden0_no,
                                              kden_no, NO_S[i], K_Nitr_NO, dissolved_oxygen[i], KDO_NO, DO_Cr_d,
                                              DO_Opt_d)
    Denit_R_S_T_NO = Denit_R_S_NO * theta_deni ** (temperature[i] - 20)

    return Denit_R_N_T_NO, Denit_R_S_T_NO


def _calculate_nitrogen_oxide(i: int, loone_nchla_fns: object, temperature: list, T_opt_NO: float, T_min_NO: float,
                              T_max_NO: float, photo_period: pd.DataFrame, rad: list, Kw_NO: float, Kc_NO: float,
                              Sim_Chla: list, z: list, K1_NO: float, K2_NO: float, NO_N: list, t: np.ndarray,
                              external_nitrite_nitrate: list, atmospheric_deposition_north: float, Q_N2S: list,
                              Nit_R_N_T_NO: list, Denit_R_N_T_NO: list, ammonium_north_observed: list,
                              volume_north: list, G_max_NO: float, f_P_N_NO: list, f_N_N_NO: list, K_NH_NO: float,
                              K_TN_NO: float, YNOChla_NO: float, Sim_Chla_N: list, S_NO_NO: float, Area_N: float,
                              NO_Temp_Adj: list, atmospheric_deposition_south: float, Q_O_M: list, NO_S: list,
                              Nit_R_S_T_NO: list, Denit_R_S_T_NO: list, ammonium_south_observed: list,
                              volume_south: list, Sim_Chla_S: list, Area_S: float, NO_MEAN: list,
                              f_P_S_NO: list, f_N_S_NO: list) -> None:
    """
    Calculate the nitrogen oxide values for the given index.

    Args:
        i (int): The current index.
        loone_nchla_fns (object): The LOONE NChla functions object.
        temperature (list): The temperature list.
        T_opt_NO (float): The T_opt_NO value.
        T_min_NO (float): The T_min_NO value.
        T_max_NO (float): The T_max_NO value.
        photo_period (pd.DataFrame): The photo period DataFrame.
        rad (list): The radiation list.
        Kw_NO (float): The Kw_NO value.
        Kc_NO (float): The Kc_NO value.
        Sim_Chla (list): The simulated chlorophyll-a list.
        z (list): The depth list.
        K1_NO (float): The K1_NO value.
        K2_NO (float): The K2_NO value.
        NO_N (list): The NO in the north.
        t (np.ndarray): The time array.
        external_nitrite_nitrate (list): The external nitrite nitrate list.
        atmospheric_deposition_north (float): The atmospheric deposition in the north.
        Q_N2S (list): The Q_N2S list.
        Nit_R_N_T_NO (list): The nitrification rate in the north.
        Denit_R_N_T_NO (list): The denitrification rate in the north.
        ammonium_north_observed (list): The observed ammonium in the north.
        volume_north (list): The volume in the north.
        G_max_NO (float): The maximum growth rate for NO.
        f_P_N_NO (list): The phosphorus limitation factor for NO in the north.
        f_N_N_NO (list): The nitrogen limitation factor for NO in the north.
        K_NH_NO (float): The K_NH_NO value.
        K_TN_NO (float): The K_TN_NO value.
        YNOChla_NO (float): The YNOChla_NO value.
        Sim_Chla_N (list): The simulated chlorophyll-a in the north.
        S_NO_NO (float): The S_NO_NO value.
        Area_N (float): The area in the north.
        NO_Temp_Adj (list): The temperature adjustment for NO.
        atmospheric_deposition_south (float): The atmospheric deposition in the south.
        Q_O_M (list): The Q_O_M list.
        NO_S (list): The NO in the south.
        Nit_R_S_T_NO (list): The nitrification rate in the south.
        Denit_R_S_T_NO (list): The denitrification rate in the south.
        ammonium_south_observed (list): The observed ammonium in the south.
        volume_south (list): The volume in the south.
        Sim_Chla_S (list): The simulated chlorophyll-a in the south.
        Area_S (float): The area in the south.
        NO_MEAN (list): The mean NO.
        f_P_S_NO (list): The phosphorus limitation factor for NO in the south.
        f_N_S_NO (list): The nitrogen limitation factor for NO in the south.

    Returns:
        None
    """
    fT_NO = loone_nchla_fns.f_T_alt1(temperature[i], T_opt_NO, T_min_NO, T_max_NO)
    fL_NO = loone_nchla_fns.f_L_alt1(photo_period['Data'].iloc[i], rad[i], Kw_NO, Kc_NO, Sim_Chla[i], z[i],
                                     K1_NO, K2_NO)

    DfEq_Res_N = odeint(loone_nchla_fns.NOx_N_DiffEq, NO_N[i], t,
                        args=(external_nitrite_nitrate[i], atmospheric_deposition_north, Q_N2S[i], Nit_R_N_T_NO[i],
                              Denit_R_N_T_NO[i], ammonium_north_observed[i], volume_north[i], G_max_NO, fT_NO,
                              fL_NO, f_P_N_NO[i], f_N_N_NO[i], K_NH_NO, K_TN_NO, YNOChla_NO, Sim_Chla_N[i],
                              S_NO_NO, Area_N, NO_Temp_Adj[i], ))
    NO_N[i + 1] = DfEq_Res_N[:, 0][1]
    DfEq_Res_S = odeint(loone_nchla_fns.NOx_S_DiffEq, NO_S[i], t,
                        args=(atmospheric_deposition_south, Q_N2S[i], Q_O_M[i], NO_N[i], Nit_R_S_T_NO[i],
                              Denit_R_S_T_NO[i], ammonium_south_observed[i], volume_south[i], G_max_NO, fT_NO,
                              fL_NO, f_P_S_NO[i], f_N_S_NO[i], K_NH_NO, K_TN_NO, YNOChla_NO, Sim_Chla_S[i],
                              S_NO_NO, Area_S, NO_Temp_Adj[i], ))
    NO_S[i + 1] = DfEq_Res_S[:, 0][1]
    NO_MEAN[i + 1] = (NO_N[i + 1] + NO_S[i + 1]) * 0.5


def _calculate_no_loads(i: int, NO_S: list, s77_outflow: list, s308_outflow: list,
                        total_regional_outflow_south: list) -> tuple:
    """
    Calculate the NO loads for the given index.

    Args:
        i (int): The current index.
        NO_S (list): The NO in the south.
        s77_outflow (list): The S77 outflow list.
        s308_outflow (list): The S308 outflow list.
        total_regional_outflow_south (list): The total regional outflow in the south.

    Returns:
        tuple: A tuple containing the NO loads for Cal, StL, and South.
    """
    NO_Load_Cal = s77_outflow[i] * NO_S[i]  # mg/d P
    NO_Load_StL = s308_outflow[i] * NO_S[i]  # mg/d P
    NO_Load_South = total_regional_outflow_south[i] * 1233.48 * NO_S[i]  # mg/d P

    return NO_Load_Cal, NO_Load_StL, NO_Load_South


def _calculate_chlorophyll_a(i: int, loone_nchla_fns: object, temperature: list, T_opt_Chla: float, T_min_Chla: float,
                             T_max_Chla: float, nitro_model_output: pd.DataFrame, photo_period: pd.DataFrame,
                             rad: list, Kw_Chla: float, Kc_Chla: float, Sim_Chla: list, z: list, K1_Chla: float,
                             K2_Chla: float, External_Chla_M: list, Q_N2S: list, Sim_Chla_N: list, vc: float,
                             K_m_T: list, K_r_T: list, G_max_Chla: float, f_P_N_Chla: list, f_N_N_Chla: list,
                             volume_north: list, Q_O_M: list, Sim_Chla_S: list, f_P_S_Chla: list, f_N_S_Chla: list,
                             volume_south: list) -> None:
    """
    Calculate the chlorophyll-a values for the given index.

    Args:
        i (int): The current index.
        loone_nchla_fns (object): The LOONE NChla functions object.
        temperature (list): The temperature list.
        T_opt_Chla (float): The T_opt_Chla value.
        T_min_Chla (float): The T_min_Chla value.
        T_max_Chla (float): The T_max_Chla value.
        nitro_model_output (pd.DataFrame): The nitro model output DataFrame.
        photo_period (pd.DataFrame): The photo period DataFrame.
        rad (list): The radiation list.
        Kw_Chla (float): The Kw_Chla value.
        Kc_Chla (float): The Kc_Chla value.
        Sim_Chla (list): The simulated chlorophyll-a list.
        z (list): The depth list.
        K1_Chla (float): The K1_Chla value.
        K2_Chla (float): The K2_Chla value.
        External_Chla_M (list): The external chlorophyll-a list.
        Q_N2S (list): The Q_N2S list.
        Sim_Chla_N (list): The simulated chlorophyll-a in the north.
        vc (float): The vc value.
        K_m_T (list): The K_m_T value.
        K_r_T (list): The K_r_T value.
        G_max_Chla (float): The maximum growth rate for chlorophyll-a.
        f_P_N_Chla (list): The phosphorus limitation factor for chlorophyll-a in the north.
        f_N_N_Chla (list): The nitrogen limitation factor for chlorophyll-a in the north.
        volume_north (list): The volume in the north.
        Q_O_M (list): The Q_O_M list.
        Sim_Chla_S (list): The simulated chlorophyll-a in the south.
        f_P_S_Chla (list): The phosphorus limitation factor for chlorophyll-a in the south.
        f_N_S_Chla (list): The nitrogen limitation factor for chlorophyll-a in the south.
        volume_south (list): The volume in the south.

    Returns:
        None
    """
    fT_Chla = loone_nchla_fns.f_T__Chla_alt1(temperature[i], T_opt_Chla, T_min_Chla, T_max_Chla, nitro_model_output['date'].iloc[i].month)
    fL_Chla = loone_nchla_fns.f_L_alt1(photo_period['Data'].iloc[i], rad[i], Kw_Chla, Kc_Chla, Sim_Chla[i], z[i], K1_Chla, K2_Chla) * 1.2 if nitro_model_output['date'].iloc[i].month in (6, 7, 8, 9, 10) else loone_nchla_fns.f_L_alt1(photo_period['Data'].iloc[i], rad[i], Kw_Chla, Kc_Chla, Sim_Chla[i], z[i], K1_Chla, K2_Chla) * 1

    Sim_Chla_N[i + 1] = loone_nchla_fns.Chla_N_alt1(External_Chla_M[i], Q_N2S[i], Sim_Chla_N[i], vc, z[i], K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla, fL_Chla, f_P_N_Chla[i], f_N_N_Chla[i], volume_north[i])
    Sim_Chla_S[i + 1] = loone_nchla_fns.Chla_S_alt1(Q_N2S[i], Q_O_M[i], Sim_Chla_N[i], Sim_Chla_S[i], vc, z[i], K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla, fL_Chla, f_P_S_Chla[i], f_N_S_Chla[i], volume_south[i])
    Sim_Chla[i + 1] = (Sim_Chla_N[i + 1] + Sim_Chla_S[i + 1]) / 2


def _calculate_chla_loads(i: int, Sim_Chla_S: list, s77_outflow: list, s308_outflow: list,
                          total_regional_outflow_south: list) -> tuple:
    """
    Calculate the chlorophyll-a loads for the given index.

    Args:
        i (int): The current index.
        Sim_Chla_S (list): The simulated chlorophyll-a in the south.
        s77_outflow (list): The S77 outflow list.
        s308_outflow (list): The S308 outflow list.
        total_regional_outflow_south (list): The total regional outflow in the south.

    Returns:
        tuple: A tuple containing the chlorophyll-a loads for Cal, StL, and South.
    """
    Chla_Load_Cal = s77_outflow[i] * Sim_Chla_S[i]  # mg/d P
    Chla_Load_StL = s308_outflow[i] * Sim_Chla_S[i]  # mg/d P
    Chla_Load_South = total_regional_outflow_south[i] * 1233.48 * Sim_Chla_S[i]  # mg/d P

    return Chla_Load_Cal, Chla_Load_StL, Chla_Load_South


def LOONE_WQ(workspace: str, photo_period_filename: str = 'PhotoPeriod', forecast_mode: bool = False, ensemble_number: int = None) -> list:
    """Daily Nitrate-Nitrite NOx and Chlorophyll-a Modeling

    Args:
        workspace (str): The working directory path.
        photo_period_filename (str): The filename of the photo period data. Default is 'PhotoPeriod'.
        forecast_mode (bool): The forecast mode flag.
        ensemble (int): The ensemble number for forecast mode.

    Returns:
        list: The results of the simulation.
    """
    # Read in the config file
    config = load_config(workspace)

    # Initialize the Data object
    data = DClass(workspace, forecast_mode, ensemble_number)

    # Read Required Data
    flow_path = os.path.join(workspace, f'LO_Inflows_BK_forecast_{ensemble_number:02}.csv' if forecast_mode else 'LO_Inflows_BK.csv')

    data_dict = _load_data(workspace, flow_path, forecast_mode, photo_period_filename, config, ensemble_number)

    inflows = data_dict['inflows']
    temperature_data = data_dict['temperature_data']
    dissolved_oxygen = data_dict['dissolved_oxygen']
    radiation_data = data_dict['radiation_data']
    storage_data = data_dict['storage_data']
    chlorophyll_a_north_data = data_dict['chlorophyll_a_north_data']
    chlorophyll_a_south_data = data_dict['chlorophyll_a_south_data']
    external_nitrate_loadings = data_dict['external_nitrate_loadings']
    s65e_nitrate_data = data_dict['s65e_nitrate_data']
    chlorophyll_a_loads_in = data_dict['chlorophyll_a_loads_in']
    s65e_chlorophyll_a_data = data_dict['s65e_chlorophyll_a_data']
    photo_period = data_dict['photo_period']
    lo_orthophosphate_north_data = data_dict['lo_orthophosphate_north_data']
    lo_orthophosphate_south_data = data_dict['lo_orthophosphate_south_data']
    lo_dissolved_inorganic_nitrogen_north_data = data_dict['lo_dissolved_inorganic_nitrogen_north_data']
    lo_dissolved_inorganic_nitrogen_south_data = data_dict['lo_dissolved_inorganic_nitrogen_south_data']

    date_start = inflows['date'].iloc[0]

    rad = radiation_data['Mean_RADT'].astype(float) * 4.6 * 1000

    temperature = temperature_data['Water_T'].astype(float)
    dissolved_oxygen = dissolved_oxygen['Mean_DO'].astype(float)
    # external_nitrite_nitrate = external_nitrate_loadings['External_NO_Ld_mg'].astype(float)  # mg
    # s65e_nitrite_nitrate = (s65e_nitrate_data[s65e_nitrate_data['date'] >= date_start]['Data'] * 1000).astype(float).tolist()  # mg/m3

    # external_chlorophyll_a = chlorophyll_a_loads_in['Chla_Loads'].astype(float) * 3  # mg
    # s65e_chlorophyll_a = s65e_chlorophyll_a_data[s65e_chlorophyll_a_data['date'] >= date_start]['Data'].astype(float).tolist()

    # N-S Procedure
    N_Per = 0.43
    S_Per = 0.57

    storage_dev = data.Storage_dev_df['DS_dev'].astype(float)  # acft
    # q_i = inflows['Inflows_cmd'].astype(float)  # m3

    # Simulated Q
    # S77_Q = LOONE_Q_Outputs['S77_Q'] * 0.0283168 * 86400  # cfs to cubic meters per day
    # S308_Q = LOONE_Q_Outputs['S308_Q'] * 0.0283168 * 86400  # cfs to cubic meters per day
    # TotRegSo = LOONE_Q_Outputs['TotRegSo']  # acft/day
    # TotRegEW = LOONE_Q_Outputs['TotRegEW']  # acft/day

    # Simulated Stage and Storage
    # Stage = LOONE_Q_Outputs['Stage'] * 0.3048  # ft to m
    # volume = LOONE_Q_Outputs['Storage'] * 1233.48  # acft to m3

    # Observed S77 S308 South
    # TODO: This should have ensembles
    outflows_observed = pd.read_csv(os.path.join(workspace, config['outflows_observed']))
    s77_outflow = outflows_observed['S77_Out']
    s308_outflow = outflows_observed['S308_Out']
    total_regional_outflow_south = outflows_observed[['S351_Out', 'S354_Out', 'S352_Out', 'L8_Out']].sum(axis=1) / 1233.48    # m3/day to acft

    # Observed Stage and Storage
    stage = storage_data['Stage_ft'].astype(float) * 0.3048  # m
    volume = storage_data['Storage_cmd'].astype(float)  # m3

    volume_north = volume * N_Per
    volume_south = volume * S_Per

    # q_o = s77_outflow + s308_outflow + total_regional_outflow_south * 1233.48  # cmd

    #TODO: This is reading in mostly just historical values
    observed_values = _calculate_observed_values(lo_dissolved_inorganic_nitrogen_north_data, lo_dissolved_inorganic_nitrogen_south_data, chlorophyll_a_north_data, chlorophyll_a_south_data, lo_orthophosphate_north_data, lo_orthophosphate_south_data)

    ammonium_north_observed = observed_values['ammonium_north']
    ammonium_south_observed = observed_values['ammonium_south']
    ammonium_observed = observed_values['ammonium']
    nitrogen_oxide_north_observed = observed_values['nitrogen_oxide_north']
    nitrogen_oxide_south_observed = observed_values['nitrogen_oxide_south']
    nitrogen_oxide_observed = observed_values['nitrogen_oxide']
    chlorophyll_a_north = observed_values['chlorophyll_a_north']
    chlorophyll_a_south = observed_values['chlorophyll_a_south']
    chlorophyll_a = observed_values['chlorophyll_a']
    dissolved_inorganic_nitrogen_north_observed = observed_values['dissolved_inorganic_nitrogen_north']
    dissolved_inorganic_nitrogen_south_observed = observed_values['dissolved_inorganic_nitrogen_south']
    dissolved_inorganic_nitrogen_observed = observed_values['dissolved_inorganic_nitrogen']
    dissolved_inorganic_phosphorus_north_observed = observed_values['dissolved_inorganic_phosphorus_north']
    dissolved_inorganic_phosphorus_south_observed = observed_values['dissolved_inorganic_phosphorus_south']
    dissolved_inorganic_phosphorus_observed = observed_values['dissolved_inorganic_phosphorus']

    nit_denit_opt = 'Opt1'

    # NO Atmospheric Deposition
    atmospheric_deposition = (714.6 * 1000 * 1000)  # mg/day
    atmospheric_deposition_north = N_Per * atmospheric_deposition
    atmospheric_deposition_south = S_Per * atmospheric_deposition

    # Calibration parameters
    calibration_parameter_no = data.Cal_Par['Par_NO']
    calibration_parameter_chla = data.Cal_Par['Par_Chla']

    # nitrification
    # Either a simple method where Nit_R_T = nitrification rate in water column (Tot_nit_R) * Temp Coeff for nitrification * (T-20)
    # Or a more detailed method R = k0 + k * fam * fox
    # two options for fam and fox
    # Option1
    # fam = (Cam/Ks+Cam) and fox = (Cox/Ksox+Cox)
    # Option 2
    # Fam = Cam and fox = (1-foxmin) * (Cox-Coxc/Coxo-Coxc) ^10a + foxmin
    # The nitrification rate for the simple method
    total_nitrification_rate = 0.3
    # nitrification parameters for the detailed method
    nitrification_parameters = _get_nitrification_parameters(calibration_parameter_no, calibration_parameter_chla)
    kni0_no = nitrification_parameters["kni0_no"]
    kni0_chla = nitrification_parameters["kni0_chla"]
    kni_no = nitrification_parameters["kni_no"]
    kni_chla = nitrification_parameters["kni_chla"]

    theta_ni = 1.06

    # Denitrification
    # Either a simple method where Denit_R_T = denitrification rate in water column (Tot_denit_R) * Temp Coeff for denitrification * (T-20)
    # Or a more detailed method R = k0 + k * fni * fox
    # two options for fam and fox
    # Option1
    # fni = (Cni/Ks+Cni) and fox = 1- (Cox/Ksox+Cox)
    # Option 2
    # Fni = Cni and fox = (Coxc-Cox/Coxc-Coxo)
    # The denitrification rate for the simple method
    total_denitrification_rate = 0.1
    # The denitrification rate for the simple method
    denitrification_parameters = _get_denitrification_parameters(calibration_parameter_no, calibration_parameter_chla)
    kden0_no = denitrification_parameters["kden0_no"]
    kden0_chla = denitrification_parameters["kden0_chla"]
    kden_no = denitrification_parameters["kden_no"]
    kden_chla = denitrification_parameters["kden_chla"]

    theta_deni = 1.06

    # Light Attenuation due to water Kw and algae Kc
    light_attenuation = _get_light_attenuation(calibration_parameter_no, calibration_parameter_chla)
    Kw_NO = light_attenuation["Kw_NO"]
    Kw_Chla = light_attenuation["Kw_Chla"]
    Kc_NO = light_attenuation["Kc_NO"]
    Kc_Chla = light_attenuation["Kc_Chla"]

    # Light limitation and inhibition coefficients K1 and K2
    light_limitation = _get_light_limitation(calibration_parameter_no, calibration_parameter_chla)
    K1_NO = light_limitation["K1_NO"]
    K1_Chla = light_limitation["K1_Chla"]
    K2_NO = light_limitation["K2_NO"]
    K2_Chla = light_limitation["K2_Chla"]

    # Maximum Phytoplankton growth rate
    phytoplankton_growth = _get_phytoplankton_growth(calibration_parameter_no, calibration_parameter_chla)
    G_max_NO = phytoplankton_growth["G_max_NO"]
    G_max_Chla = phytoplankton_growth["G_max_Chla"]

    # Half Saturation Coefficients
    half_saturation = _get_half_saturation_coefficients(calibration_parameter_no, calibration_parameter_chla)
    K_NH_NO = half_saturation["K_NH_NO"]
    K_NH_Chla = half_saturation["K_NH_Chla"]
    K_Nitr_NO = half_saturation["K_Nitr_NO"]
    K_Nitr_Chla = half_saturation["K_Nitr_Chla"]
    K_DIN_NO = half_saturation["K_DIN_NO"]
    K_DIN_Chla = half_saturation["K_DIN_Chla"]
    KP_NO = half_saturation["KP_NO"]
    KP_Chla = half_saturation["KP_Chla"]
    K_TN_NO = half_saturation["K_TN_NO"]
    K_TN_Chla = half_saturation["K_TN_Chla"]
    YNOChla_NO = half_saturation["YNOChla_NO"]
    YNOChla_Chla = half_saturation["YNOChla_Chla"]

    # Temperatures
    temperatures = _get_temperatures(calibration_parameter_no, calibration_parameter_chla)
    T_opt_NO = temperatures["T_opt_NO"]
    T_min_NO = temperatures["T_min_NO"]
    T_max_NO = temperatures["T_max_NO"]
    T_opt_Chla = temperatures["T_opt_Chla"]
    T_min_Chla = temperatures["T_min_Chla"]
    T_max_Chla = temperatures["T_max_Chla"]

    # Dissolved Oxygen
    dissolved_oxygen_params = _calculate_dissolved_oxygen(calibration_parameter_no, calibration_parameter_chla)
    KDO_NO = dissolved_oxygen_params["KDO_NO"]
    KDO_Chla = dissolved_oxygen_params["KDO_Chla"]

    # Sediment Release
    sediment_release = _calculate_sediment_release(calibration_parameter_no, calibration_parameter_chla)
    S_NO_NO = sediment_release["S_NO_NO"]
    S_NO_Chla = sediment_release["S_NO_Chla"]

    Theta_NO = 1.06
    NO_Temp_Adj = Theta_NO ** (temperature - 20)
    Area = 1730 * 1E6  # m2
    Area_N = Area * N_Per
    Area_S = Area * S_Per

    X = len(inflows.index)

    # Nitrogen and Phosphorus Limiting
    nutrient_limitation_factors = _calculate_nutrient_limitation_factors(dissolved_inorganic_phosphorus_observed, dissolved_inorganic_nitrogen_observed, dissolved_inorganic_phosphorus_north_observed, dissolved_inorganic_phosphorus_south_observed, dissolved_inorganic_nitrogen_north_observed, dissolved_inorganic_nitrogen_south_observed, KP_NO, KP_Chla, K_DIN_NO, K_DIN_Chla)

    f_P_NO = nutrient_limitation_factors["f_P_NO"]
    f_P_Chla = nutrient_limitation_factors["f_P_Chla"]
    f_N_NO = nutrient_limitation_factors["f_N_NO"]
    f_N_Chla = nutrient_limitation_factors["f_N_Chla"]
    f_P_N_NO = nutrient_limitation_factors["f_P_N_NO"]
    f_P_N_Chla = nutrient_limitation_factors["f_P_N_Chla"]
    f_P_S_NO = nutrient_limitation_factors["f_P_S_NO"]
    f_P_S_Chla = nutrient_limitation_factors["f_P_S_Chla"]
    f_N_N_NO = nutrient_limitation_factors["f_N_N_NO"]
    f_N_N_Chla = nutrient_limitation_factors["f_N_N_Chla"]
    f_N_S_NO = nutrient_limitation_factors["f_N_S_NO"]
    f_N_S_Chla = nutrient_limitation_factors["f_N_S_Chla"]

    Q_I_M = np.zeros(X, dtype=object)
    Q_O_M = np.zeros(X, dtype=object)
    Q_N2S = np.zeros(X, dtype=object)
    External_NO_M = np.zeros(X, dtype=object)
    External_Chla_M = np.zeros(X, dtype=object)

    Nit_R_N_NO = np.zeros(X, dtype=object)
    Nit_R_N_T_NO = np.zeros(X, dtype=object)
    Nit_R_S_NO = np.zeros(X, dtype=object)
    Nit_R_S_T_NO = np.zeros(X, dtype=object)
    Denit_R_N_NO = np.zeros(X, dtype=object)
    Denit_R_N_T_NO = np.zeros(X, dtype=object)
    Denit_R_S_NO = np.zeros(X, dtype=object)
    Denit_R_S_T_NO = np.zeros(X, dtype=object)

    fT_NO = np.zeros(X, dtype=object)
    fL_NO = np.zeros(X, dtype=object)

    fT_Chla = np.zeros(X, dtype=object)
    fL_Chla = np.zeros(X, dtype=object)

    NO_N = np.zeros(X, dtype=object)
    NO_S = np.zeros(X, dtype=object)
    NO_MEAN = np.zeros(X, dtype=object)

    NO_Load_Cal = np.zeros(X, dtype=object)
    NO_Load_StL = np.zeros(X, dtype=object)
    NO_Load_South = np.zeros(X, dtype=object)

    Sim_Chla_N = np.zeros(X, dtype=object)
    Sim_Chla_S = np.zeros(X, dtype=object)
    Sim_Chla = np.zeros(X, dtype=object)

    Chla_Load_Cal = np.zeros(X, dtype=object)
    Chla_Load_StL = np.zeros(X, dtype=object)
    Chla_Load_South = np.zeros(X, dtype=object)

    # Water depth could be stage - B.L. (variable) or a constant average depth (e.g., 2.5 m)

    z = np.zeros(X, dtype=object)
    # B_L = 1  # m
    # z = Stage - B_L  # m

    for i in range(X):
        z[i] = 2.5

    fox_min = calibration_parameter_no[19]
    DO_Cr_n = calibration_parameter_no[20]
    DO_Opt_n = calibration_parameter_no[21]
    a = calibration_parameter_no[22]
    DO_Cr_d = calibration_parameter_no[23]
    DO_Opt_d = calibration_parameter_no[24]

    vc = calibration_parameter_chla[28]  # phytoplankton settling velocity (m/d).

    K_r = calibration_parameter_chla[25]
    Theta_r = 1.06
    K_r_T = K_r * Theta_r ** (temperature - 20)

    YNHChla = calibration_parameter_chla[26]  # 0.05
    # Chla mortality
    k_m = calibration_parameter_chla[27]
    Theta_m = 1.06
    K_m_T = k_m * Theta_m ** (temperature - 20)

    Im = calibration_parameter_chla[31]

    Grazing = 0
    Theta_G = 1.06
    Grazing_T = Grazing * Theta_G ** (temperature - 20)
    # Initial values
    NO_N[0] = nitrogen_oxide_north_observed[0]
    NO_S[0] = nitrogen_oxide_south_observed[0]

    NO_MEAN[0] = (NO_N[0] + NO_S[0]) * 0.5

    Sim_Chla_N[0] = chlorophyll_a_north[0]
    Sim_Chla_S[0] = chlorophyll_a_south[0]
    Sim_Chla[0] = chlorophyll_a[0]

    nitro_model_output = pd.DataFrame(inflows['date'], columns=['date'])
    nitro_model_output['date'] = pd.to_datetime(nitro_model_output['date'])
    print("LOONE Nitrogen Module is Running!")
    # Filter all datasets by date_start
    datasets = [s65e_nitrate_data, s65e_chlorophyll_a_data, chlorophyll_a_loads_in, inflows, outflows_observed, external_nitrate_loadings]
    for df in datasets:
        df.drop(df[df['date'] < date_start].index, inplace=True)

    # Merge all data on date - this ensures that the dates will line up
    merged = inflows[['date', 'Inflows_cmd']].merge(
        chlorophyll_a_loads_in[['date', 'Chla_Loads']], on='date'
    ).merge(
        s65e_chlorophyll_a_data[['date', 'Data']], on='date', suffixes=('', '_s65e_chla')
    ).merge(
        s65e_nitrate_data[['date', 'Data']], on='date', suffixes=('', '_s65e_nitrate')
    ).merge(
        outflows_observed[['date', 'S77_Out', 'S308_Out', 'S351_Out', 'S354_Out', 'S352_Out', 'L8_Out']], on='date'
    ).merge(
        external_nitrate_loadings[['date', 'External_NO_Ld_mg']], on='date'
    )

    # Rename for clarity
    merged.rename(columns={
        'Data': 'S65E_Chla',
        'Data_s65e_nitrate': 'S65E_NO'
    }, inplace=True)

    # Compute q_o (outflows in m³/day)
    merged['total_regional_outflow_south'] = merged[['S351_Out', 'S354_Out', 'S352_Out', 'L8_Out']].sum(axis=1) / 1233.48
    merged['q_o'] = merged['S77_Out'] + merged['S308_Out'] + merged['total_regional_outflow_south'] * 1233.48

    # Prepare input lists
    q_i = merged['Inflows_cmd'].astype(float).tolist()
    q_o = merged['q_o'].astype(float).tolist()
    s65e_nitrite_nitrate = (merged['S65E_NO'] * 1000).astype(float).tolist()
    s65e_chlorophyll_a = merged['S65E_Chla'].astype(float).tolist()
    external_chlorophyll_a = (merged['Chla_Loads'].astype(float) * 3).tolist()
    external_nitrite_nitrate = merged['External_NO_Ld_mg'].astype(float).tolist()
    for i in range(len(merged.index) - 1):
        # print(Nitro_Model_Output['date'].iloc[i])

        Q_I_M[i], Q_O_M[i], External_NO_M[i], External_Chla_M[i] = _calculate_inflows_and_outflows(i, storage_dev,
                                                                                                   q_i, q_o,
                                                                                                   external_nitrite_nitrate,
                                                                                                   s65e_nitrite_nitrate,
                                                                                                   external_chlorophyll_a,
                                                                                                   s65e_chlorophyll_a)

        Q_N2S[i] = (Q_I_M[i] * 1 + Q_O_M[i] * 0)

        t = np.linspace(1, 2, num=2)

        # Nitrification calcs
        Nit_R_N_T_NO[i], Nit_R_S_T_NO[i] = _calculate_nitrification_rates(i, loone_nchla_fns, nit_denit_opt,
                                                                          total_nitrification_rate, kni0_no, kni_no,
                                                                          ammonium_north_observed,
                                                                          ammonium_south_observed, K_NH_NO,
                                                                          dissolved_oxygen, KDO_NO, fox_min, DO_Cr_n,
                                                                          DO_Opt_n, a, theta_ni, temperature)

        # Denitrification calcs
        Denit_R_N_T_NO[i], Denit_R_S_T_NO[i] = _calculate_denitrification_rates(i, loone_nchla_fns, nit_denit_opt,
                                                                                total_denitrification_rate, kden0_no,
                                                                                kden_no, NO_N, NO_S, K_Nitr_NO,
                                                                                dissolved_oxygen, KDO_NO, DO_Cr_d,
                                                                                DO_Opt_d, theta_deni, temperature)

        _calculate_nitrogen_oxide(i, loone_nchla_fns, temperature, T_opt_NO, T_min_NO, T_max_NO, photo_period, rad,
                                  Kw_NO, Kc_NO, Sim_Chla, z, K1_NO, K2_NO, NO_N, t, external_nitrite_nitrate,
                                  atmospheric_deposition_north, Q_N2S, Nit_R_N_T_NO, Denit_R_N_T_NO,
                                  ammonium_north_observed, volume_north, G_max_NO, f_P_N_NO, f_N_N_NO, K_NH_NO, K_TN_NO,
                                  YNOChla_NO, Sim_Chla_N, S_NO_NO, Area_N, NO_Temp_Adj, atmospheric_deposition_south,
                                  Q_O_M, NO_S, Nit_R_S_T_NO, Denit_R_S_T_NO, ammonium_south_observed, volume_south,
                                  Sim_Chla_S, Area_S, NO_MEAN, f_P_S_NO, f_N_S_NO)

        NO_Load_Cal[i], NO_Load_StL[i], NO_Load_South[i] = _calculate_no_loads(i, NO_S, s77_outflow, s308_outflow,
                                                                               total_regional_outflow_south)

        ##### Chla
        _calculate_chlorophyll_a(i, loone_nchla_fns, temperature, T_opt_Chla, T_min_Chla, T_max_Chla,
                                 nitro_model_output, photo_period, rad, Kw_Chla, Kc_Chla, Sim_Chla, z,
                                 K1_Chla, K2_Chla, External_Chla_M, Q_N2S, Sim_Chla_N, vc, K_m_T, K_r_T, G_max_Chla,
                                 f_P_N_Chla, f_N_N_Chla, volume_north, Q_O_M, Sim_Chla_S, f_P_S_Chla, f_N_S_Chla,
                                 volume_south)

        Chla_Load_Cal[i], Chla_Load_StL[i], Chla_Load_South[i] = _calculate_chla_loads(i, Sim_Chla_S, s77_outflow,
                                                                                       s308_outflow,
                                                                                       total_regional_outflow_south)

    print("Exporting Module Outputs!")

    nitro_model_output, nitro_mod_out_m = _update_nitro_model_output(nitro_model_output, NO_N, NO_S, NO_MEAN, Sim_Chla,
                                                                     Sim_Chla_N, Sim_Chla_S, temperature)

    # Constituent Loads
    constit_loads_df, constit_loads_m = _prepare_constit_loads_df(inflows, NO_Load_Cal, NO_Load_StL, NO_Load_South,
                                                                  Chla_Load_Cal, Chla_Load_StL, Chla_Load_South)

    # Summer month loads
    Smr_Mnth_NOx_StL_arr, Smr_Mnth_Chla_StL_arr, Smr_Mnth_NOx_Cal_arr, Smr_Mnth_Chla_Cal_arr = _calculate_summer_month_loads(constit_loads_m)

    if config['sim_type'] in [0, 1, 3]:
        variables_dict = {
            'Constit_Loads': constit_loads_df,
            'Nitro_Model_Output': nitro_model_output,
            'Constit_Loads_M': constit_loads_m,
            'Nitro_Mod_Out_M': nitro_mod_out_m,
        }
        return_list = list(variables_dict.values())

        #TODO: Will this end up having ensembles after we have the forecast data reading in correctly?
        for k, v in variables_dict.items():
            file_name = f'{k}_forecast' if forecast_mode else k
            v.to_csv(os.path.join(workspace, f'{file_name}.csv'), index=False)
        return return_list
    else:
        return_list = [Smr_Mnth_NOx_StL_arr, Smr_Mnth_Chla_StL_arr, Smr_Mnth_NOx_Cal_arr, Smr_Mnth_Chla_Cal_arr,
                       nitro_model_output]
    
    return return_list

    # Algae_Opt_Mnth_NOx_StL = []
    # Algae_Opt_Mnth_NOx_Cal = []
    # Algae_Opt_Mnth_Chla_StL = []
    # Algae_Opt_Mnth_Chla_Cal = []

    # for i in range(len(Constit_Loads_M.index)):
    #     if Nitro_Model_Output['Water Temp'].iloc[i] >= 25:
    #         Algae_Opt_Mnth_NOx_StL.append(Constit_Loads_M['NO_Load_StL'].iloc[i])
    #         Algae_Opt_Mnth_Chla_StL.append(Constit_Loads_M['Chla_Load_StL'].iloc[i])
    #         Algae_Opt_Mnth_NOx_Cal.append(Constit_Loads_M['NO_Load_Cal'].iloc[i])
    #         Algae_Opt_Mnth_Chla_Cal.append(Constit_Loads_M['Chla_Load_Cal'].iloc[i])

    # Algae_Opt_Mnth_NOx_StL_arr = np.asarray(Algae_Opt_Mnth_NOx_StL)
    # Algae_Opt_Mnth_Chla_StL_arr = np.asarray(Algae_Opt_Mnth_Chla_StL)

    # Algae_Opt_Mnth_NOx_Cal_arr = np.asarray(Algae_Opt_Mnth_NOx_Cal)
    # Algae_Opt_Mnth_Chla_Cal_arr = np.asarray(Algae_Opt_Mnth_Chla_Cal)
    # if config['sim_type'] in [0, 1]:
    #     return[Constit_Loads_M,Nitro_Mod_Out_M,Algae_Opt_Mnth_NOx_StL_arr,Algae_Opt_Mnth_Chla_StL_arr,Algae_Opt_Mnth_NOx_Cal_arr,Algae_Opt_Mnth_Chla_Cal_arr]
    # else:
    #     return[Algae_Opt_Mnth_NOx_StL_arr,Algae_Opt_Mnth_Chla_StL_arr,Algae_Opt_Mnth_NOx_Cal_arr,Algae_Opt_Mnth_Chla_Cal_arr,Nitro_Model_Output]


if __name__ == "__main__":
    # Parse Arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "workspace",
        nargs=1,
        help="The path to the working directory.",
    )

    args = argparser.parse_args()
    workspace = args.workspace[0]

    # Run LOONE_WQ
    LOONE_WQ(workspace)
