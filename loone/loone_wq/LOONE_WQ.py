# Daily Nitrate-Nitrite NOx and Chlorophyll-a Modeling

# import required packages
import os
import argparse
import pandas as pd
import numpy as np
from scipy.integrate import odeint
from loone.utils import loone_nchla_fns, load_config
from loone.data import Data as DClass


def LOONE_WQ(workspace: str, photo_period_filename: str = 'PhotoPeriod', forecast_mode: bool = False) -> list:
    """Daily Nitrate-Nitrite NOx and Chlorophyll-a Modeling

    Args:
        workspace (str): The working directory path.
        photo_period_filename (str): The filename of the photo period data. Default is 'PhotoPeriod'.
        forecast_mode (bool): The forecast mode flag.

    Returns:
        list: The results of the simulation.
    """
    # Read in the config file
    config = load_config(workspace)

    # Initialize the Data object
    data = DClass(workspace)

    # Read Required Data
    flow_path = os.path.join(workspace, 'LO_Inflows_BK_forecast.csv' if forecast_mode else 'LO_Inflows_BK.csv')

    inflows = pd.read_csv(os.path.join(workspace, flow_path))
    temperature_data = pd.read_csv(os.path.join(workspace, 'Filled_WaterT.csv'))
    dissolved_oxygen = pd.read_csv(os.path.join(workspace, 'LO_DO_Clean_daily.csv'))
    radiation_data = pd.read_csv(os.path.join(workspace, 'LO_RADT_data.csv'))
    storage_data = pd.read_csv(os.path.join(workspace, 'Average_LO_Storage_3MLag.csv'))
    chlorophyll_a_north_data = pd.read_csv(os.path.join(workspace, 'N_Merged_Chla.csv'))  # microgram/L
    chlorophyll_a_south_data = pd.read_csv(os.path.join(workspace, 'S_Merged_Chla.csv'))  # microgram/L
    external_nitrate_loadings = pd.read_csv(os.path.join(workspace, 'LO_External_Loadings_NO.csv'))  # mg
    s65e_basename = 'water_quality_S65E_NITRATE+NITRITE-N_Interpolated_forecast.csv' if forecast_mode else 'water_quality_S65E_NITRATE+NITRITE-N_Interpolated.csv'
    s65e_nitrate_data = pd.read_csv(os.path.join(workspace, s65e_basename))  # mg/m3
    chlorophyll_a_loads_in = pd.read_csv(os.path.join(workspace, 'Chla_Loads_In.csv'))  # mg
    s65e_chlorophyll_a_basename = 'S65E_Chla_Merged_forecast.csv' if forecast_mode else 'S65E_Chla_Merged.csv'
    s65e_chlorophyll_a_data = pd.read_csv(os.path.join(workspace, s65e_chlorophyll_a_basename))  # mg/m3

    photoperiod = pd.read_csv(os.path.join(workspace, f'{photo_period_filename}.csv'))
    lo_orthophosphate_north_data = pd.read_csv(os.path.join(workspace, 'N_OP.csv'))  # mg/m3
    lo_orthophosphate_south_data = pd.read_csv(os.path.join(workspace, 'S_OP.csv'))  # mg/m3
    lo_dissolved_inorganic_nitrogen_north_data = pd.read_csv(os.path.join(workspace, 'N_DIN.csv'))  # mg/m3
    lo_dissolved_inorganic_nitrogen_south_data = pd.read_csv(os.path.join(workspace, 'S_DIN.csv'))  # mg/m3

    date_start = inflows['date'].iloc[0]

    rad = radiation_data['Mean_RADT'].astype(float) * 4.6 * 1000

    temperature = temperature_data['Water_T'].astype(float)
    dissolved_oxygen = dissolved_oxygen['Mean_DO'].astype(float)
    external_nitrite_nitrate = external_nitrate_loadings['External_NO_Ld_mg'].astype(float)  # mg
    s65e_nitrite_nitrate = (s65e_nitrate_data[s65e_nitrate_data['date'] >= date_start]['Data'] * 1000).astype(float).tolist()  # mg/m3

    external_chlorophyll_a = chlorophyll_a_loads_in['Chla_Loads'].astype(float) * 3  # mg
    s65e_chlorophyll_a = s65e_chlorophyll_a_data[s65e_chlorophyll_a_data['date'] >= date_start]['Data'].astype(float).tolist()

    # N-S Procedure
    N_Per = 0.43
    S_Per = 0.57

    storage_dev = data.Storage_dev_df['DS_dev'].astype(float)  # acft
    q_i = inflows['Inflows_cmd'].astype(float)  # m3

    # Simulated Q
    # S77_Q = LOONE_Q_Outputs['S77_Q'] * 0.0283168 * 86400  # cfs to cubic meters per day
    # S308_Q = LOONE_Q_Outputs['S308_Q'] * 0.0283168 * 86400  # cfs to cubic meters per day
    # TotRegSo = LOONE_Q_Outputs['TotRegSo']  # acft/day
    # TotRegEW = LOONE_Q_Outputs['TotRegEW']  # acft/day

    # Simulated Stage and Storage
    # Stage = LOONE_Q_Outputs['Stage'] * 0.3048  # ft to m
    # volume = LOONE_Q_Outputs['Storage'] * 1233.48  # acft to m3

    # Observed S77 S308 South
    outflows_observed = pd.read_csv(os.path.join(workspace, 'Flow_df_3MLag.csv'))
    s77_outflow = outflows_observed['S77_Out']
    s308_outflow = outflows_observed['S308_Out']
    total_regional_outflow_south = outflows_observed[['S351_Out', 'S354_Out', 'S352_Out', 'L8_Out']].sum(axis=1) / 1233.48    # m3/day to acft

    # Observed Stage and Storage
    stage = storage_data['Stage_ft'].astype(float) * 0.3048  # m
    volume = storage_data['Storage_cmd'].astype(float)  # m3

    volume_north = volume * N_Per
    volume_south = volume * S_Per

    q_o = s77_outflow + s308_outflow + total_regional_outflow_south * 1233.48  # cmd

    ammonium_north_observed = lo_dissolved_inorganic_nitrogen_north_data['NH4'].astype(float)  # mg/m3
    ammonium_south_observed = lo_dissolved_inorganic_nitrogen_south_data['NH4'].astype(float)  # mg/m3
    ammonium_observed = (ammonium_north_observed + ammonium_south_observed) / 2
    nitrogen_oxide_north_observed = lo_dissolved_inorganic_nitrogen_north_data['NO'].astype(float)  # mg/m3
    nitrogen_oxide_south_observed = lo_dissolved_inorganic_nitrogen_south_data['NO'].astype(float)  # mg/m3
    nitrogen_oxide_observed = (nitrogen_oxide_north_observed + nitrogen_oxide_south_observed) / 2
    chlorophyll_a_north = chlorophyll_a_north_data['Chla'].astype(float)  # microgram/L = mg/m3
    chlorophyll_a_south = chlorophyll_a_south_data['Chla'].astype(float)  # microgram/L = mg/m3
    chlorophyll_a = (chlorophyll_a_north + chlorophyll_a_south) / 2
    dissolved_inorganic_nitrogen_north_observed = lo_dissolved_inorganic_nitrogen_north_data['DIN'].astype(float)
    dissolved_inorganic_nitrogen_south_observed = lo_dissolved_inorganic_nitrogen_south_data['DIN'].astype(float)
    dissolved_inorganic_nitrogen_observed = (dissolved_inorganic_nitrogen_north_observed +
                                             dissolved_inorganic_nitrogen_south_observed) / 2
    dissolved_inorganic_phosphorus_north_observed = lo_orthophosphate_north_data['OP'].astype(float)
    dissolved_inorganic_phosphorus_south_observed = lo_orthophosphate_south_data['OP'].astype(float)
    dissolved_inorganic_phosphorus_observed = (dissolved_inorganic_phosphorus_north_observed +
                                               dissolved_inorganic_phosphorus_south_observed) / 2

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
    kni0_no = calibration_parameter_no[0]
    kni0_chla = calibration_parameter_chla[0]
    kni_no = calibration_parameter_no[1]
    kni_chla = calibration_parameter_chla[1]

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
    kden0_no = calibration_parameter_no[2]
    kden0_chla = calibration_parameter_chla[2]

    kden_no = calibration_parameter_no[3]
    kden_chla = calibration_parameter_chla[3]

    theta_deni = 1.06

    # Light Attenuation due to water Kw and algae Kc
    Kw_NO = calibration_parameter_no[4]  # 1.7
    Kw_Chla = calibration_parameter_chla[4]  # 1.7
    Kc_NO = calibration_parameter_no[5]  # 0.029
    Kc_Chla = calibration_parameter_chla[5]  # 0.029
    # Light limitation and inhibition coefficients K1 and K2
    K1_NO = calibration_parameter_no[6]  # 100
    K1_Chla = calibration_parameter_chla[6]  # 100
    K2_NO = calibration_parameter_no[7]  # 700
    K2_Chla = calibration_parameter_chla[7]  # 700
    # Maximum Phytoplankton growth rate
    G_max_NO = calibration_parameter_no[8]  # 1.5
    G_max_Chla = calibration_parameter_chla[8]  # 1.5

    # Half Saturation Coefficients
    K_NH_NO = calibration_parameter_no[9]  # 0.1
    K_NH_Chla = calibration_parameter_chla[9]  # 0.1
    K_Nitr_NO = calibration_parameter_no[10]
    K_Nitr_Chla = calibration_parameter_chla[10]
    K_DIN_NO = K_NH_NO + K_Nitr_NO
    K_DIN_Chla = K_NH_Chla + K_Nitr_Chla
    KP_NO = calibration_parameter_no[11]
    KP_Chla = calibration_parameter_chla[11]
    K_TN_NO = calibration_parameter_no[12]  # 0.1
    K_TN_Chla = calibration_parameter_chla[12]  # 0.1
    YNOChla_NO = calibration_parameter_no[13]  # 0.1
    YNOChla_Chla = calibration_parameter_chla[13]  # 0.1

    # Temperatures
    T_opt_NO = calibration_parameter_no[14]  # 15 # C
    T_min_NO = calibration_parameter_no[15]  # 10 # C
    T_max_NO = calibration_parameter_no[16]  # 30 # C

    T_opt_Chla = calibration_parameter_chla[14]  # 15 # C
    T_min_Chla = calibration_parameter_chla[15]  # 10 # C
    T_max_Chla = calibration_parameter_chla[16]  # 30 # C

    # Dissolved Oxygen
    KDO_NO = calibration_parameter_no[17]
    KDO_Chla = calibration_parameter_chla[17]

    # Sediment Release
    S_NO_NO = calibration_parameter_no[18]  # mg/m2/d
    S_NO_Chla = calibration_parameter_chla[18]  # mg/m2/d

    Theta_NO = 1.06
    NO_Temp_Adj = Theta_NO ** (temperature - 20)
    Area = 1730 * 1E6  # m2
    Area_N = Area * N_Per
    Area_S = Area * S_Per

    X = len(inflows.index)

    # Nitrogen and Phosphorus Limiting
    f_P_NO = dissolved_inorganic_phosphorus_observed / (KP_NO + dissolved_inorganic_phosphorus_observed)
    f_P_Chla = dissolved_inorganic_phosphorus_observed / (KP_Chla + dissolved_inorganic_phosphorus_observed)

    f_N_NO = dissolved_inorganic_nitrogen_observed / (K_DIN_NO + dissolved_inorganic_nitrogen_observed)
    f_N_Chla = dissolved_inorganic_nitrogen_observed / (K_DIN_Chla + dissolved_inorganic_nitrogen_observed)

    f_P_N_NO = dissolved_inorganic_phosphorus_north_observed / (KP_NO + dissolved_inorganic_phosphorus_north_observed)
    f_P_N_Chla = dissolved_inorganic_phosphorus_north_observed / (KP_Chla + dissolved_inorganic_phosphorus_north_observed)

    f_P_S_NO = dissolved_inorganic_phosphorus_south_observed / (KP_NO + dissolved_inorganic_phosphorus_south_observed)
    f_P_S_Chla = dissolved_inorganic_phosphorus_south_observed / (KP_Chla + dissolved_inorganic_phosphorus_south_observed)

    f_N_N_NO = dissolved_inorganic_nitrogen_north_observed / (K_DIN_NO + dissolved_inorganic_nitrogen_north_observed)
    f_N_N_Chla = dissolved_inorganic_nitrogen_north_observed / (K_DIN_Chla + dissolved_inorganic_nitrogen_north_observed)

    f_N_S_NO = dissolved_inorganic_nitrogen_south_observed / (K_DIN_NO + dissolved_inorganic_nitrogen_south_observed)
    f_N_S_Chla = dissolved_inorganic_nitrogen_south_observed / (K_DIN_Chla + dissolved_inorganic_nitrogen_south_observed)

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
    for i in range(X - 1):
        # print(Nitro_Model_Output['date'].iloc[i])

        if storage_dev[i] >= 0:
            Q_I_M[i] = q_i[i] + storage_dev[i] * 1233.48  # m3/d
            Q_O_M[i] = q_o[i]

            External_NO_M[i] = external_nitrite_nitrate[i] + Q_I_M[i] * s65e_nitrite_nitrate[i]
            External_Chla_M[i] = external_chlorophyll_a[i] + Q_I_M[i] * s65e_chlorophyll_a[i]

        else:
            Q_O_M[i] = q_o[i] - storage_dev[i] * 1233.48  # m3/d
            Q_I_M[i] = q_i[i]
            External_NO_M[i] = external_nitrite_nitrate[i]
            External_Chla_M[i] = external_chlorophyll_a[i]

        Q_N2S[i] = (Q_I_M[i] * 1 + Q_O_M[i] * 0)

        t = np.linspace(1, 2, num=2)
        Nit_R_N_NO[i] = loone_nchla_fns.Nit_Rate('%s' % nit_denit_opt, total_nitrification_rate, kni0_no, kni_no,
                                                 ammonium_north_observed[i], K_NH_NO, dissolved_oxygen[i], KDO_NO,
                                                 fox_min, DO_Cr_n, DO_Opt_n, a)
        Nit_R_N_T_NO[i] = Nit_R_N_NO[i] * theta_ni ** (temperature[i] - 20)
        Nit_R_S_NO[i] = loone_nchla_fns.Nit_Rate('%s' % nit_denit_opt, total_nitrification_rate, kni0_no, kni_no,
                                                 ammonium_south_observed[i], K_NH_NO, dissolved_oxygen[i], KDO_NO,
                                                 fox_min, DO_Cr_n, DO_Opt_n, a)
        Nit_R_S_T_NO[i] = Nit_R_S_NO[i] * theta_ni ** (temperature[i] - 20)

        Denit_R_N_NO[i] = loone_nchla_fns.Denit_Rate('%s' % nit_denit_opt, total_denitrification_rate, kden0_no,
                                                     kden_no, NO_N[i], K_Nitr_NO, dissolved_oxygen[i], KDO_NO, DO_Cr_d,
                                                     DO_Opt_d)
        Denit_R_N_T_NO[i] = Denit_R_N_NO[i] * theta_deni ** (temperature[i] - 20)
        Denit_R_S_NO[i] = loone_nchla_fns.Denit_Rate('%s' % nit_denit_opt, total_denitrification_rate, kden0_no,
                                                     kden_no, NO_S[i], K_Nitr_NO, dissolved_oxygen[i], KDO_NO, DO_Cr_d,
                                                     DO_Opt_d)
        Denit_R_S_T_NO[i] = Denit_R_S_NO[i] * theta_deni ** (temperature[i] - 20)

        fT_NO[i] = loone_nchla_fns.f_T_alt1(temperature[i], T_opt_NO, T_min_NO, T_max_NO)
        fL_NO[i] = loone_nchla_fns.f_L_alt1(photoperiod['Data'].iloc[i], rad[i], Kw_NO, Kc_NO, Sim_Chla[i], z[i],
                                            K1_NO, K2_NO)

        DfEq_Res_N = odeint(loone_nchla_fns.NOx_N_DiffEq, NO_N[i], t,
                            args=(external_nitrite_nitrate[i], atmospheric_deposition_north, Q_N2S[i], Nit_R_N_T_NO[i],
                                  Denit_R_N_T_NO[i], ammonium_north_observed[i], volume_north[i], G_max_NO, fT_NO[i],
                                  fL_NO[i], f_P_N_NO[i], f_N_N_NO[i], K_NH_NO, K_TN_NO, YNOChla_NO, Sim_Chla_N[i],
                                  S_NO_NO, Area_N, NO_Temp_Adj[i], ))
        NO_N[i + 1] = DfEq_Res_N[:, 0][1]
        DfEq_Res_S = odeint(loone_nchla_fns.NOx_S_DiffEq, NO_S[i], t,
                            args=(atmospheric_deposition_south, Q_N2S[i], Q_O_M[i], NO_N[i], Nit_R_S_T_NO[i],
                                  Denit_R_S_T_NO[i], ammonium_south_observed[i], volume_south[i], G_max_NO, fT_NO[i],
                                  fL_NO[i], f_P_S_NO[i], f_N_S_NO[i], K_NH_NO, K_TN_NO, YNOChla_NO, Sim_Chla_S[i],
                                  S_NO_NO, Area_S, NO_Temp_Adj[i], ))
        NO_S[i + 1] = DfEq_Res_S[:, 0][1]
        NO_MEAN[i + 1] = (NO_N[i + 1] + NO_S[i + 1]) * 0.5

        NO_Load_Cal[i] = s77_outflow[i] * NO_S[i]  # mg/d P
        NO_Load_StL[i] = s308_outflow[i] * NO_S[i]  # mg/d P
        NO_Load_South[i] = total_regional_outflow_south[i] * 1233.48 * NO_S[i]  # mg/d P

        ##### Chla
        fT_Chla[i] = loone_nchla_fns.f_T__Chla_alt1(temperature[i], T_opt_Chla, T_min_Chla, T_max_Chla,
                                                    nitro_model_output['date'].iloc[i].month)

        fL_Chla[i] = loone_nchla_fns.f_L_alt1(photoperiod['Data'].iloc[i], rad[i], Kw_Chla, Kc_Chla, Sim_Chla[i], z[i], K1_Chla, K2_Chla) * 1.2 if nitro_model_output['date'].iloc[i].month in (6, 7, 8, 9, 10) else loone_nchla_fns.f_L_alt1(photoperiod['Data'].iloc[i], rad[i], Kw_Chla, Kc_Chla, Sim_Chla[i], z[i], K1_Chla, K2_Chla) * 1

        Sim_Chla_N[i + 1] = loone_nchla_fns.Chla_N_alt1(External_Chla_M[i], Q_N2S[i], Sim_Chla_N[i], vc, z[i],
                                                        K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla[i], fL_Chla[i],
                                                        f_P_N_Chla[i], f_N_N_Chla[i], volume_north[i])
        Sim_Chla_S[i + 1] = loone_nchla_fns.Chla_S_alt1(Q_N2S[i], Q_O_M[i], Sim_Chla_N[i], Sim_Chla_S[i], vc, z[i],
                                                        K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla[i], fL_Chla[i],
                                                        f_P_S_Chla[i], f_N_S_Chla[i], volume_south[i])
        Sim_Chla[i + 1] = (Sim_Chla_N[i + 1] + Sim_Chla_S[i + 1]) / 2

        Chla_Load_Cal[i] = s77_outflow[i] * Sim_Chla_S[i]  # mg/d P
        Chla_Load_StL[i] = s308_outflow[i] * Sim_Chla_S[i]  # mg/d P
        Chla_Load_South[i] = total_regional_outflow_south[i] * 1233.48 * Sim_Chla_S[i]  # mg/d P

    print("Exporting Module Outputs!")

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

    if config['sim_type'] in [0, 1]:
        variables_dict = {
            'Constit_Loads': constit_loads_df,
            'Nitro_Model_Output': nitro_model_output,
            'Constit_Loads_M': constit_loads_m,
            'Nitro_Mod_Out_M': nitro_mod_out_m,
        }
        return_list = list(variables_dict.values())

        for k, v in variables_dict.items():
            file_name = f'{k}_forecast' if forecast_mode else k
            v.to_csv(os.path.join(workspace, f'{file_name}.csv'), index=False)
        return return_list
    else:
        return_list = [Smr_Mnth_NOx_StL_arr, Smr_Mnth_Chla_StL_arr, Smr_Mnth_NOx_Cal_arr, Smr_Mnth_Chla_Cal_arr,
                       nitro_model_output]

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
