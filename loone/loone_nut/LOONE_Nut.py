# Standard library imports
from datetime import datetime, timedelta
import os

# Third-party library imports
import numpy as np
import pandas as pd

# Local application/library specific imports
from data.data import Data
from data.pre_defined_variables import Pre_defined_Variables
from data.tp_variables_regions import TP_Variables
from loone_config.Model_Config import Model_Config
import utils.stg_sto_ar
import utils.tp_mass_balance_functions_regions as TP_MBFR

SECONDS_IN_DAY = 86400
CUBIC_METERS_IN_ACRE_FOOT = 1233.48
SQUARE_METERS_IN_ACRE = 4046.85642
METERS_IN_FOOT = 0.3048
METERS_IN_FOOT_3 = 0.305
SECONDS_IN_HOUR = 3600
HOURS_IN_DAY = 24
CUBIC_METERS_IN_CUBIC_FOOT = 0.028316847
MILLIGRAMS_IN_TON = 1e9

def LOONE_Nut(
    loone_q_path: str,
    loads_external_filename: str,
    flow_df_filename: str,
    data_dir: str | None = None,
) -> pd.DataFrame:
    """
    Simulates nutrient (phosphorus) dynamics in the water column

    Args:
        loone_q_path (str): Path to LOONE Q CSV file.
        loads_external_filename (str): The name of the file that holds the external loads for the model.
        flow_df_filename (str): The name of the file that holds the flow data for the model.
        data_dir (str | None, optional): Path to the data directory. Defaults to Model_Config.Working_Path.

    Returns:
        pd.DataFrame: A Dataframe containing an estimate of the total phosphorus concentration in the lake for a certain time series.
    """
    print("LOONE Nut Module is Running!")
    
    # Read in config values
    
    
    data_dir = data_dir if data_dir else Model_Config.Working_Path
    loone_q = pd.read_csv(loone_q_path)
    # Based on the defined Start and End year, month, and day on the
    # Pre_defined_Variables File, Startdate and enddate are defined.
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    startdate = datetime.now().date()  # datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    year, month, day = map(int, Pre_defined_Variables.enddate_entry)
    enddate = startdate + timedelta(days=15)  # datetime(year, month, day).date()

    date_rng_0 = pd.date_range(start=startdate, end=enddate, freq="D")
    load_ext = pd.read_csv(os.path.join(data_dir, loads_external_filename))
    q_in = pd.read_csv(os.path.join(data_dir, f"LO_Inflows_BK.csv"))
    flow_df = pd.read_csv(os.path.join(data_dir, flow_df_filename))
    q_o = flow_df["Outflows"].values
    s77_q = loone_q["S77_Q"].values
    s308_q = loone_q["S308_Q"].values
    tot_reg_so = flow_df[["S351_Out", "S352_Out", "S354_Out"]].sum(axis=1) * (
        70.0456 / SECONDS_IN_DAY
    )
    sto_stage = pd.read_csv(os.path.join(data_dir, f"Average_LO_Storage_3MLag.csv"))
    stage_lo = sto_stage["Stage_ft"].values
    storage = sto_stage["Storage_acft"].values
    n_rows = len(q_in.index)
    storage_dev = Data.Stroage_dev_df["DS_dev"]
    l_ext = load_ext["TP_Loads_In_mg"]  # mg
    atm_dep_n = TP_Variables.N_Per * load_ext["Atm_Loading_mg"]
    atm_dep_s = TP_Variables.S_Per * load_ext["Atm_Loading_mg"]

    # Read Shear Stress driven by Wind Speed
    wind_shear_str = pd.read_csv(os.path.join(data_dir, f"WindShearStress.csv"))
    w_ss = wind_shear_str["ShearStress"]  # Dyne/cm2
    nu_ts = pd.read_csv(os.path.join(data_dir, "nu.csv"))
    LO_BL = 0.5  # m (Bed Elevation of LO)
    g = 9.8  # m/s2 gravitational acceleration
    cal_res = pd.read_csv(os.path.join(data_dir, f"nondominated_Sol_var.csv"))
    par = cal_res["Par"]
    d_c = par[20]  # m (particle diameter 10 microm /1E6 to convert to m) clay
    d_s = par[21]  # m sand
    nu_d = nu_ts["nu"]

    R = 1.65  # submerged specific gravity (1.65 for quartz in water)
    c_1_c = par[16]
    c_2_c = par[17]
    c_1_s = par[18]
    c_2_s = par[19]

    # Parameters associated with sediment resuspension
    E_0 = 1e-4
    E_1 = 2
    E_2 = 3
    crtcl_sh_str = par[22]  # 0.32 #Dyne/cm2
    td = par[23]  # days
    l_ext_m = np.zeros(n_rows, dtype=object)
    q_n2s = np.zeros(n_rows, dtype=object)
    stage_2_ar = np.zeros(n_rows, dtype=object)
    lo_wd = np.zeros(n_rows, dtype=object)
    lake_o_storage_n = np.zeros(n_rows, dtype=object)
    lake_o_storage_s = np.zeros(n_rows, dtype=object)
    lake_o_a_n = np.zeros(n_rows, dtype=object)
    lake_o_a_s = np.zeros(n_rows, dtype=object)
    lake_o_a_m_n = np.zeros(n_rows, dtype=object)
    lake_o_a_s_n = np.zeros(n_rows, dtype=object)
    lake_o_a_r_n = np.zeros(n_rows, dtype=object)
    lake_o_a_p_n = np.zeros(n_rows, dtype=object)
    lake_o_a_m_s = np.zeros(n_rows, dtype=object)
    lake_o_a_s_s = np.zeros(n_rows, dtype=object)
    lake_o_a_r_s = np.zeros(n_rows, dtype=object)
    lake_o_a_p_s = np.zeros(n_rows, dtype=object)
    dip_lake_n = np.zeros(n_rows, dtype=object)
    dip_lake_s = np.zeros(n_rows, dtype=object)
    tp_lake_mean = np.zeros(n_rows, dtype=object)
    j_des_m_n = np.zeros(n_rows, dtype=object)
    j_des_s_n = np.zeros(n_rows, dtype=object)
    j_des_r_n = np.zeros(n_rows, dtype=object)
    j_des_p_n = np.zeros(n_rows, dtype=object)
    j_des_m_s = np.zeros(n_rows, dtype=object)
    j_des_s_s = np.zeros(n_rows, dtype=object)
    j_des_r_s = np.zeros(n_rows, dtype=object)
    j_des_p_s = np.zeros(n_rows, dtype=object)
    j_ads_m_n = np.zeros(n_rows, dtype=object)
    j_ads_s_n = np.zeros(n_rows, dtype=object)
    j_ads_r_n = np.zeros(n_rows, dtype=object)
    j_ads_p_n = np.zeros(n_rows, dtype=object)
    j_ads_m_s = np.zeros(n_rows, dtype=object)
    j_ads_s_s = np.zeros(n_rows, dtype=object)
    j_ads_r_s = np.zeros(n_rows, dtype=object)
    j_ads_p_s = np.zeros(n_rows, dtype=object)
    p_sed_m_n = np.zeros(n_rows, dtype=object)
    p_sed_s_n = np.zeros(n_rows, dtype=object)
    p_sed_r_n = np.zeros(n_rows, dtype=object)
    p_sed_p_n = np.zeros(n_rows, dtype=object)
    p_sed_m_s = np.zeros(n_rows, dtype=object)
    p_sed_s_s = np.zeros(n_rows, dtype=object)
    p_sed_r_s = np.zeros(n_rows, dtype=object)
    p_sed_p_s = np.zeros(n_rows, dtype=object)
    j_sedburial_m_n = np.zeros(n_rows, dtype=object)
    j_sedburial_s_n = np.zeros(n_rows, dtype=object)
    j_sedburial_r_n = np.zeros(n_rows, dtype=object)
    j_sedburial_p_n = np.zeros(n_rows, dtype=object)
    j_sedburial_m_s = np.zeros(n_rows, dtype=object)
    j_sedburial_s_s = np.zeros(n_rows, dtype=object)
    j_sedburial_r_s = np.zeros(n_rows, dtype=object)
    j_sedburial_p_s = np.zeros(n_rows, dtype=object)
    j_Γburial_m_n = np.zeros(n_rows, dtype=object)
    j_Γburial_s_n = np.zeros(n_rows, dtype=object)
    j_Γburial_r_n = np.zeros(n_rows, dtype=object)
    j_Γburial_p_n = np.zeros(n_rows, dtype=object)
    j_Γburial_m_s = np.zeros(n_rows, dtype=object)
    j_Γburial_s_s = np.zeros(n_rows, dtype=object)
    j_Γburial_r_s = np.zeros(n_rows, dtype=object)
    j_Γburial_p_s = np.zeros(n_rows, dtype=object)
    Γ_m_n = np.zeros(n_rows, dtype=object)
    Γ_s_n = np.zeros(n_rows, dtype=object)
    Γ_r_n = np.zeros(n_rows, dtype=object)
    Γ_p_n = np.zeros(n_rows, dtype=object)
    Γ_m_s = np.zeros(n_rows, dtype=object)
    Γ_s_s = np.zeros(n_rows, dtype=object)
    Γ_r_s = np.zeros(n_rows, dtype=object)
    Γ_p_s = np.zeros(n_rows, dtype=object)
    dip_pore_m_n = np.zeros(n_rows, dtype=object)
    dip_pore_s_n = np.zeros(n_rows, dtype=object)
    dip_pore_r_n = np.zeros(n_rows, dtype=object)
    dip_pore_p_n = np.zeros(n_rows, dtype=object)
    dip_pore_m_s = np.zeros(n_rows, dtype=object)
    dip_pore_s_s = np.zeros(n_rows, dtype=object)
    dip_pore_r_s = np.zeros(n_rows, dtype=object)
    dip_pore_p_s = np.zeros(n_rows, dtype=object)
    tp_lake_n = np.zeros(n_rows, dtype=object)
    tp_lake_s = np.zeros(n_rows, dtype=object)
    sed_resusp_m_n = np.zeros(n_rows, dtype=object)
    sed_resusp_s_n = np.zeros(n_rows, dtype=object)
    sed_resusp_r_n = np.zeros(n_rows, dtype=object)
    sed_resusp_p_n = np.zeros(n_rows, dtype=object)
    sed_resusp_m_s = np.zeros(n_rows, dtype=object)
    sed_resusp_s_s = np.zeros(n_rows, dtype=object)
    sed_resusp_r_s = np.zeros(n_rows, dtype=object)
    sed_resusp_p_s = np.zeros(n_rows, dtype=object)

    j_decomp_m_n = np.zeros(n_rows, dtype=object)
    j_decomp_s_n = np.zeros(n_rows, dtype=object)
    j_decomp_r_n = np.zeros(n_rows, dtype=object)
    j_decomp_p_n = np.zeros(n_rows, dtype=object)
    j_decomp_m_s = np.zeros(n_rows, dtype=object)
    j_decomp_s_s = np.zeros(n_rows, dtype=object)
    j_decomp_r_s = np.zeros(n_rows, dtype=object)
    j_decomp_p_s = np.zeros(n_rows, dtype=object)

    settling_p_n = np.zeros(n_rows, dtype=object)
    settling_p_s = np.zeros(n_rows, dtype=object)

    p_diff_m_n = np.zeros(n_rows, dtype=object)
    p_diff_s_n = np.zeros(n_rows, dtype=object)
    p_diff_r_n = np.zeros(n_rows, dtype=object)
    p_diff_p_n = np.zeros(n_rows, dtype=object)
    p_diff_m_s = np.zeros(n_rows, dtype=object)
    p_diff_s_s = np.zeros(n_rows, dtype=object)
    p_diff_r_s = np.zeros(n_rows, dtype=object)
    p_diff_p_s = np.zeros(n_rows, dtype=object)

    q_i = q_in["Inflows_cmd"]
    q_i_m = np.zeros(n_rows, dtype=object)
    q_o = np.zeros(n_rows, dtype=object)
    q_o_m = np.zeros(n_rows, dtype=object)

    p_load_cal = np.zeros(n_rows, dtype=object)
    p_load_stl = np.zeros(n_rows, dtype=object)
    p_load_south = np.zeros(n_rows, dtype=object)

    v_settle_n = np.zeros(n_rows, dtype=object)
    v_settle_s = np.zeros(n_rows, dtype=object)

    ##Initial Values##
    # S.A. is calculated based on the Lake's previous time step Stage, but for
    # the S.A. at i=0 I used same time step Stage!
    start_storage = utils.stg_sto_ar.stg2sto(Pre_defined_Variables.startstage, 0)
    stage_2_ar[0] = utils.stg_sto_ar.stg2ar(stage_lo[0], 0)
    stage_2_ar[1] = utils.stg_sto_ar.stg2ar(stage_lo[1], 0)
    storage[0] = start_storage  # ac-ft
    storage[1] = utils.stg_sto_ar.stg2sto(stage_lo[1], 0)  # ac-ft
    # TP_MassBalanceModel Initial Values.
    tp_lake_n[0] = 225  # mg/m3
    tp_lake_s[0] = 275  # mg/m3
    tp_lake_mean[0] = (tp_lake_n[0] + tp_lake_s[0]) / 2
    Γ_m_n[0] = 25  # mg/kg
    Γ_s_n[0] = 25  # mg/kg
    Γ_r_n[0] = 25  # mg/kg
    Γ_p_n[0] = 25  # mg/kg
    Γ_m_s[0] = 25  # mg/kg
    Γ_s_s[0] = 25  # mg/kg
    Γ_r_s[0] = 25  # mg/kg
    Γ_p_s[0] = 25  # mg/kg
    dip_pore_m_n[0] = 700  # 760 #mg/m3
    dip_pore_s_n[0] = 240  # 205 #mg/m3
    dip_pore_r_n[0] = 240  # 205 #mg/m3
    dip_pore_p_n[0] = 160  # 160 #mg/m3
    dip_pore_m_s[0] = 700  # 760 #mg/m3
    dip_pore_s_s[0] = 240  # 205 #mg/m3
    dip_pore_r_s[0] = 240  # 205 #mg/m3
    dip_pore_p_s[0] = 160  # 160 #mg/m3
    p_sed_m_n[0] = 1100  # mg/kg
    p_sed_s_n[0] = 300  # mg/kg
    p_sed_r_n[0] = 300  # mg/kg
    p_sed_p_n[0] = 200  # mg/kg
    p_sed_m_s[0] = 1100  # mg/kg
    p_sed_s_s[0] = 300  # mg/kg
    p_sed_r_s[0] = 300  # mg/kg
    p_sed_p_s[0] = 200  # mg/kg
    Θ_m = _calculate_porosity(TP_Variables.Bulk_density_M, TP_Variables.Particle_density_M, TP_Variables.Per_H2O_M)
    Θ_s = _calculate_porosity(TP_Variables.Bulk_density_S, TP_Variables.Particle_density_S, TP_Variables.Per_H2O_S)
    Θ_r = _calculate_porosity(TP_Variables.Bulk_density_R, TP_Variables.Particle_density_R, TP_Variables.Per_H2O_R)
    Θ_p = _calculate_porosity(TP_Variables.Bulk_density_P, TP_Variables.Particle_density_P, TP_Variables.Per_H2O_P)
    
    # Mass of sediment in surfacial mix Mud layer in the North Region(kg)
    mass_sed_m_n = _calculate_mass_sediment(TP_Variables.A_Mud_N, TP_Variables.Z_sed, TP_Variables.Per_H2O_M, TP_Variables.Bulk_density_M)

    # Mass of sediment in surfacial mix Sand layer in the North Region(kg)
    mass_sed_s_n = _calculate_mass_sediment(TP_Variables.A_Sand_N, TP_Variables.Z_sed, TP_Variables.Per_H2O_S, TP_Variables.Bulk_density_S)

    # Mass of sediment in surfacial mix Rock layer in the North Region(kg)
    mass_sed_r_n = _calculate_mass_sediment(TP_Variables.A_Rock_N, TP_Variables.Z_sed, TP_Variables.Per_H2O_R, TP_Variables.Bulk_density_R)

    # Mass of sediment in surfacial mix Peat layer in the North Region(kg)
    mass_sed_p_n = _calculate_mass_sediment(TP_Variables.A_Peat_N, TP_Variables.Z_sed, TP_Variables.Per_H2O_P, TP_Variables.Bulk_density_P)

    # Mass of sediment in surfacial mix Mud layer in the South Region(kg)
    mass_sed_m_s = _calculate_mass_sediment(TP_Variables.A_Mud_S, TP_Variables.Z_sed, TP_Variables.Per_H2O_M, TP_Variables.Bulk_density_M)

    # Mass of sediment in surfacial mix Sand layer in the South Region(kg)
    mass_sed_s_s = _calculate_mass_sediment(TP_Variables.A_Sand_S, TP_Variables.Z_sed, TP_Variables.Per_H2O_S, TP_Variables.Bulk_density_S)

    # Mass of sediment in surfacial mix Rock layer in the South Region(kg)
    mass_sed_r_s = _calculate_mass_sediment(TP_Variables.A_Rock_S, TP_Variables.Z_sed, TP_Variables.Per_H2O_R, TP_Variables.Bulk_density_R)

    # Mass of sediment in surfacial mix Peat layer in the South Region(kg)
    mass_sed_p_s = _calculate_mass_sediment(TP_Variables.A_Peat_S, TP_Variables.Z_sed, TP_Variables.Per_H2O_P, TP_Variables.Bulk_density_P)

    for i in range(n_rows - 2):
        if storage_dev[i] >= 0:
            q_i_m[i] = q_i[i] + storage_dev[i] * CUBIC_METERS_IN_ACRE_FOOT  # m3/d
            q_o_m[i] = q_o[i]
            l_ext_m[i] = l_ext[i] + q_i_m[i] * tp_lake_n[i]
        else:
            q_o_m[i] = q_o[i] - storage_dev[i] * CUBIC_METERS_IN_ACRE_FOOT  # m3/d
            q_i_m[i] = q_i[i]
            l_ext_m[i] = l_ext[i]
        q_n2s[i] = (q_i_m[i] + q_o_m[i]) / 2
        stage_2_ar[i + 2] = utils.stg_sto_ar.stg2ar(stage_lo[i + 2], 0)
        lo_wd[i] = stage_lo[i] * METERS_IN_FOOT - LO_BL
        lake_o_storage_n[i] = storage[i] * TP_Variables.N_Per * SQUARE_METERS_IN_ACRE * METERS_IN_FOOT_3  # m3
        lake_o_storage_s[i] = storage[i] * TP_Variables.S_Per * SQUARE_METERS_IN_ACRE * METERS_IN_FOOT_3  # m3
        lake_o_a_n[i] = stage_2_ar[i] * TP_Variables.N_Per * SQUARE_METERS_IN_ACRE  # m2
        lake_o_a_s[i] = stage_2_ar[i] * TP_Variables.S_Per * SQUARE_METERS_IN_ACRE  # m2
        lake_o_a_m_n[i] = _calculate_sediment_area(lake_o_a_n[i], TP_Variables.A_Mud_N, TP_Variables.A_Mud_N, TP_Variables.A_Sand_N, TP_Variables.A_Rock_N, TP_Variables.A_Peat_N)
        lake_o_a_s_n[i] = _calculate_sediment_area(lake_o_a_n[i], TP_Variables.A_Sand_N, TP_Variables.A_Mud_N, TP_Variables.A_Sand_N, TP_Variables.A_Rock_N, TP_Variables.A_Peat_N)
        lake_o_a_r_n[i] = _calculate_sediment_area(lake_o_a_n[i], TP_Variables.A_Rock_N, TP_Variables.A_Mud_N, TP_Variables.A_Sand_N, TP_Variables.A_Rock_N, TP_Variables.A_Peat_N)
        lake_o_a_p_n[i] = _calculate_sediment_area(lake_o_a_n[i], TP_Variables.A_Peat_N, TP_Variables.A_Mud_N, TP_Variables.A_Sand_N, TP_Variables.A_Rock_N, TP_Variables.A_Peat_N)
        lake_o_a_m_s[i] = _calculate_sediment_area(lake_o_a_s[i], TP_Variables.A_Mud_S, TP_Variables.A_Mud_S, TP_Variables.A_Sand_S, TP_Variables.A_Rock_S, TP_Variables.A_Peat_S)
        lake_o_a_s_s[i] = _calculate_sediment_area(lake_o_a_s[i], TP_Variables.A_Sand_S, TP_Variables.A_Mud_S, TP_Variables.A_Sand_S, TP_Variables.A_Rock_S, TP_Variables.A_Peat_S)
        lake_o_a_r_s[i] = _calculate_sediment_area(lake_o_a_s[i], TP_Variables.A_Rock_S, TP_Variables.A_Mud_S, TP_Variables.A_Sand_S, TP_Variables.A_Rock_S, TP_Variables.A_Peat_S)
        lake_o_a_p_s[i] = _calculate_sediment_area(lake_o_a_s[i], TP_Variables.A_Peat_S, TP_Variables.A_Mud_S, TP_Variables.A_Sand_S, TP_Variables.A_Rock_S, TP_Variables.A_Peat_S)

        dip_lake_n[i] = TP_MBFR.DIP_Lake(tp_lake_n[i])
        dip_lake_s[i] = TP_MBFR.DIP_Lake(tp_lake_s[i])
        
        # Calculate the settling velocity for the North and South regions
        v_settle_n[i] = _calculate_settling_velocity(
            R, g, d_c, d_s, c_1_c, c_1_s, c_2_c, c_2_s, nu_d[i],
            TP_Variables.A_Mud_N, TP_Variables.A_Peat_N, TP_Variables.A_Sand_N, TP_Variables.A_Rock_N, TP_Variables.A_N
        )
        
        v_settle_s[i] = _calculate_settling_velocity(
            R, g, d_c, d_s, c_1_c, c_1_s, c_2_c, c_2_s, nu_d[i],
            TP_Variables.A_Mud_S, TP_Variables.A_Peat_S, TP_Variables.A_Sand_S, TP_Variables.A_Rock_S, TP_Variables.A_S
        )

        # Calculate desorption fluxes for the North and South regions (for mud, sand, rock, and peat)
        j_des_m_n[i] = TP_MBFR.Des_flux(Γ_m_n[i], mass_sed_m_n, TP_Variables.K_des_M)
        j_des_s_n[i] = TP_MBFR.Des_flux(Γ_s_n[i], mass_sed_s_n, TP_Variables.K_des_S)
        j_des_r_n[i] = TP_MBFR.Des_flux(Γ_r_n[i], mass_sed_r_n, TP_Variables.K_des_R)
        j_des_p_n[i] = TP_MBFR.Des_flux(Γ_p_n[i], mass_sed_p_n, TP_Variables.K_des_P)
        j_des_m_s[i] = TP_MBFR.Des_flux(Γ_m_s[i], mass_sed_m_s, TP_Variables.K_des_M)
        j_des_s_s[i] = TP_MBFR.Des_flux(Γ_s_s[i], mass_sed_s_s, TP_Variables.K_des_S)
        j_des_r_s[i] = TP_MBFR.Des_flux(Γ_r_s[i], mass_sed_r_s, TP_Variables.K_des_R)
        j_des_p_s[i] = TP_MBFR.Des_flux(Γ_p_s[i], mass_sed_p_s, TP_Variables.K_des_P)

        # Calculate adsorption fluxes for the North and South regions
        j_ads_m_n[i] = TP_MBFR.Ads_flux(
            dip_pore_m_n[i],
            Γ_m_n[i],
            mass_sed_m_n,
            TP_Variables.K_ads_M,
            TP_Variables.Γ_inf,
        )
        j_ads_s_n[i] = TP_MBFR.Ads_flux(
            dip_pore_s_n[i],
            Γ_s_n[i],
            mass_sed_s_n,
            TP_Variables.K_ads_S,
            TP_Variables.Γ_inf,
        )
        j_ads_r_n[i] = TP_MBFR.Ads_flux(
            dip_pore_r_n[i],
            Γ_r_n[i],
            mass_sed_r_n,
            TP_Variables.K_ads_R,
            TP_Variables.Γ_inf,
        )
        j_ads_p_n[i] = TP_MBFR.Ads_flux(
            dip_pore_p_n[i],
            Γ_p_n[i],
            mass_sed_p_n,
            TP_Variables.K_ads_P,
            TP_Variables.Γ_inf,
        )
        j_ads_m_s[i] = TP_MBFR.Ads_flux(
            dip_pore_m_s[i],
            Γ_m_s[i],
            mass_sed_m_s,
            TP_Variables.K_ads_M,
            TP_Variables.Γ_inf,
        )
        j_ads_s_s[i] = TP_MBFR.Ads_flux(
            dip_pore_s_s[i],
            Γ_s_s[i],
            mass_sed_s_s,
            TP_Variables.K_ads_S,
            TP_Variables.Γ_inf,
        )
        j_ads_r_s[i] = TP_MBFR.Ads_flux(
            dip_pore_r_s[i],
            Γ_r_s[i],
            mass_sed_r_s,
            TP_Variables.K_ads_R,
            TP_Variables.Γ_inf,
        )
        j_ads_p_s[i] = TP_MBFR.Ads_flux(
            dip_pore_p_s[i],
            Γ_p_s[i],
            mass_sed_p_s,
            TP_Variables.K_ads_P,
            TP_Variables.Γ_inf,
        )

        # Calculate the sediment burial fluxes for the North and South regions (for mud, sand, rock, and peat)
        j_sedburial_m_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_m_n[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_N,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        j_sedburial_s_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_s_n[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_N,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        j_sedburial_r_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_r_n[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_N,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        j_sedburial_p_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_p_n[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_N,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )
        j_sedburial_m_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_m_s[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_S,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        j_sedburial_s_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_s_s[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_S,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        j_sedburial_r_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_r_s[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_S,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        j_sedburial_p_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_p_s[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_S,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )

        # Calculate the sediment resuspension for the North and South regions (for mud, sand, rock, and peat)
        sed_resusp_m_n[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_m_n[i])
        sed_resusp_s_n[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_s_n[i])
        sed_resusp_r_n[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_r_n[i])
        sed_resusp_p_n[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_p_n[i])
        sed_resusp_m_s[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_m_s[i])
        sed_resusp_s_s[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_s_s[i])
        sed_resusp_r_s[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_r_s[i])
        sed_resusp_p_s[i] = _calculate_sediment_resuspension(E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_p_s[i])

        # Calculate the phosphorus concentration for the next time step for the North and South regions (for mud, sand, rock, and peat)
        p_sed_m_n[i + 1] = _calculate_next_sediment(lake_o_a_m_n[i], tp_lake_n[i], dip_lake_n[i], j_sedburial_m_n[i], p_sed_m_n[i], mass_sed_m_n, TP_Variables.K_decomp_M, v_settle_n[i], sed_resusp_m_n[i], lake_o_storage_n[i])
        p_sed_s_n[i + 1] = _calculate_next_sediment(lake_o_a_s_n[i], tp_lake_n[i], dip_lake_n[i], j_sedburial_s_n[i], p_sed_s_n[i], mass_sed_s_n, TP_Variables.K_decomp_S, v_settle_n[i], sed_resusp_s_n[i], lake_o_storage_n[i])
        p_sed_r_n[i + 1] = _calculate_next_sediment(lake_o_a_r_n[i], tp_lake_n[i], dip_lake_n[i], j_sedburial_r_n[i], p_sed_r_n[i], mass_sed_r_n, TP_Variables.K_decomp_R, v_settle_n[i], sed_resusp_r_n[i], lake_o_storage_n[i])
        p_sed_p_n[i + 1] = _calculate_next_sediment(lake_o_a_p_n[i], tp_lake_n[i], dip_lake_n[i], j_sedburial_p_n[i], p_sed_p_n[i], mass_sed_p_n, TP_Variables.K_decomp_P, v_settle_n[i], sed_resusp_p_n[i], lake_o_storage_n[i])
        p_sed_m_s[i + 1] = _calculate_next_sediment(lake_o_a_m_s[i], tp_lake_s[i], dip_lake_s[i], j_sedburial_m_s[i], p_sed_m_s[i], mass_sed_m_s, TP_Variables.K_decomp_M, v_settle_s[i], sed_resusp_m_s[i], lake_o_storage_s[i])
        p_sed_s_s[i + 1] = _calculate_next_sediment(lake_o_a_s_s[i], tp_lake_s[i], dip_lake_s[i], j_sedburial_s_s[i], p_sed_s_s[i], mass_sed_s_s, TP_Variables.K_decomp_S, v_settle_s[i], sed_resusp_s_s[i], lake_o_storage_s[i])
        p_sed_r_s[i + 1] = _calculate_next_sediment(lake_o_a_r_s[i], tp_lake_s[i], dip_lake_s[i], j_sedburial_r_s[i], p_sed_r_s[i], mass_sed_r_s, TP_Variables.K_decomp_R, v_settle_s[i], sed_resusp_r_s[i], lake_o_storage_s[i])
        p_sed_p_s[i + 1] = _calculate_next_sediment(lake_o_a_p_s[i], tp_lake_s[i], dip_lake_s[i], j_sedburial_p_s[i], p_sed_p_s[i], mass_sed_p_s, TP_Variables.K_decomp_P, v_settle_s[i], sed_resusp_p_s[i], lake_o_storage_s[i])

        # Calculate the burial fluxes for the North and South regions (for mud, sand, rock, and peat)
        j_Γburial_m_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_m_n[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_N,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        j_Γburial_s_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_s_n[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_N,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        j_Γburial_r_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_r_n[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_N,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        j_Γburial_p_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_p_n[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_N,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )
        j_Γburial_m_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_m_s[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_S,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        j_Γburial_s_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_s_s[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_S,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        j_Γburial_r_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_r_s[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_S,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        j_Γburial_p_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_p_s[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_S,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )

        # Calculate the phosphorus concentration for the next time step for the North and South regions (for mud, sand, rock, and peat)
        Γ_m_n[i + 1] = _calculate_next_concentration(j_ads_m_n[i], j_des_m_n[i], j_Γburial_m_n[i], Γ_m_n[i], mass_sed_m_n)
        Γ_s_n[i + 1] = _calculate_next_concentration(j_ads_s_n[i], j_des_s_n[i], j_Γburial_s_n[i], Γ_s_n[i], mass_sed_s_n)
        Γ_r_n[i + 1] = _calculate_next_concentration(j_ads_r_n[i], j_des_r_n[i], j_Γburial_r_n[i], Γ_r_n[i], mass_sed_r_n)
        Γ_p_n[i + 1] = _calculate_next_concentration(j_ads_p_n[i], j_des_p_n[i], j_Γburial_p_n[i], Γ_p_n[i], mass_sed_p_n)
        Γ_m_s[i + 1] = _calculate_next_concentration(j_ads_m_s[i], j_des_m_s[i], j_Γburial_m_s[i], Γ_m_s[i], mass_sed_m_s)
        Γ_s_s[i + 1] = _calculate_next_concentration(j_ads_s_s[i], j_des_s_s[i], j_Γburial_s_s[i], Γ_s_s[i], mass_sed_s_s)
        Γ_r_s[i + 1] = _calculate_next_concentration(j_ads_r_s[i], j_des_r_s[i], j_Γburial_r_s[i], Γ_r_s[i], mass_sed_r_s)
        Γ_p_s[i + 1] = _calculate_next_concentration(j_ads_p_s[i], j_des_p_s[i], j_Γburial_p_s[i], Γ_p_s[i], mass_sed_p_s)

        # Calculate the decomposition rates for the North and South regions (for mud, sand, rock, and peat)
        j_decomp_m_n[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_M, p_sed_m_n[i], mass_sed_m_n
        )
        j_decomp_s_n[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_S, p_sed_s_n[i], mass_sed_s_n
        )
        j_decomp_r_n[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_R, p_sed_r_n[i], mass_sed_r_n
        )
        j_decomp_p_n[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_P, p_sed_p_n[i], mass_sed_p_n
        )
        j_decomp_m_s[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_M, p_sed_m_s[i], mass_sed_m_s
        )
        j_decomp_s_s[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_S, p_sed_s_s[i], mass_sed_s_s
        )
        j_decomp_r_s[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_R, p_sed_r_s[i], mass_sed_r_s
        )
        j_decomp_p_s[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_P, p_sed_p_s[i], mass_sed_p_s
        )

        # Calculate the DIP pore concentration for the North and South regions (for mud, sand, rock, and peat)
        dip_pore_m_n[i + 1] = _calculate_DIP_pore(Θ_m, dip_pore_m_n[i], dip_lake_n[i], j_des_m_n[i], j_ads_m_n[i], p_sed_m_n[i], mass_sed_m_n, TP_Variables.v_diff_M, TP_Variables.A_Mud_N, TP_Variables.K_decomp_M, TP_Variables.v_burial_M)
        dip_pore_s_n[i + 1] = _calculate_DIP_pore(Θ_s, dip_pore_s_n[i], dip_lake_n[i], j_des_s_n[i], j_ads_s_n[i], p_sed_s_n[i], mass_sed_s_n, TP_Variables.v_diff_S, TP_Variables.A_Sand_N, TP_Variables.K_decomp_S, TP_Variables.v_burial_S)
        dip_pore_r_n[i + 1] = _calculate_DIP_pore(Θ_r, dip_pore_r_n[i], dip_lake_n[i], j_des_r_n[i], j_ads_r_n[i], p_sed_r_n[i], mass_sed_r_n, TP_Variables.v_diff_R, TP_Variables.A_Rock_N, TP_Variables.K_decomp_R, TP_Variables.v_burial_R)
        dip_pore_p_n[i + 1] = _calculate_DIP_pore(Θ_p, dip_pore_p_n[i], dip_lake_n[i], j_des_p_n[i], j_ads_p_n[i], p_sed_p_n[i], mass_sed_p_n, TP_Variables.v_diff_P, TP_Variables.A_Peat_N, TP_Variables.K_decomp_P, TP_Variables.v_burial_P)
        dip_pore_m_s[i + 1] = _calculate_DIP_pore(Θ_m, dip_pore_m_s[i], dip_lake_s[i], j_des_m_s[i], j_ads_m_s[i], p_sed_m_s[i], mass_sed_m_s, TP_Variables.v_diff_M, TP_Variables.A_Mud_S, TP_Variables.K_decomp_M, TP_Variables.v_burial_M)
        dip_pore_s_s[i + 1] = _calculate_DIP_pore(Θ_s, dip_pore_s_s[i], dip_lake_s[i], j_des_s_s[i], j_ads_s_s[i], p_sed_s_s[i], mass_sed_s_s, TP_Variables.v_diff_S, TP_Variables.A_Sand_S, TP_Variables.K_decomp_S, TP_Variables.v_burial_S)
        dip_pore_r_s[i + 1] = _calculate_DIP_pore(Θ_r, dip_pore_r_s[i], dip_lake_s[i], j_des_r_s[i], j_ads_r_s[i], p_sed_r_s[i], mass_sed_r_s, TP_Variables.v_diff_R, TP_Variables.A_Rock_S, TP_Variables.K_decomp_R, TP_Variables.v_burial_R)
        dip_pore_p_s[i + 1] = _calculate_DIP_pore(Θ_p, dip_pore_p_s[i], dip_lake_s[i], j_des_p_s[i], j_ads_p_s[i], p_sed_p_s[i], mass_sed_p_s, TP_Variables.v_diff_P, TP_Variables.A_Peat_S, TP_Variables.K_decomp_P, TP_Variables.v_burial_P)

        # Calculate the settling phosphorus for the North and South regions
        settling_p_n[i] = TP_MBFR.Sett_P(
            tp_lake_n[i],
            dip_lake_n[i],
            lake_o_a_n[i],
            lake_o_storage_n[i],
            v_settle_n[i],
        )
        settling_p_s[i] = TP_MBFR.Sett_P(
            tp_lake_s[i],
            dip_lake_s[i],
            lake_o_a_s[i],
            lake_o_storage_s[i],
            v_settle_s[i],
        )

        # Calculate the diffusion phosphorus for the North and South regions (for mud, sand, rock, and peat)
        p_diff_m_n[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_M,
            dip_pore_m_n[i],
            dip_lake_n[i],
            Θ_m,
            TP_Variables.A_Mud_N,
            lake_o_storage_n[i],
        )
        p_diff_s_n[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_S,
            dip_pore_s_n[i],
            dip_lake_n[i],
            Θ_s,
            TP_Variables.A_Sand_N,
            lake_o_storage_n[i],
        )
        p_diff_r_n[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_R,
            dip_pore_r_n[i],
            dip_lake_n[i],
            Θ_r,
            TP_Variables.A_Rock_N,
            lake_o_storage_n[i],
        )
        p_diff_p_n[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_P,
            dip_pore_p_n[i],
            dip_lake_n[i],
            Θ_p,
            TP_Variables.A_Peat_N,
            lake_o_storage_n[i],
        )
        p_diff_m_s[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_M,
            dip_pore_m_s[i],
            dip_lake_s[i],
            Θ_m,
            TP_Variables.A_Mud_S,
            lake_o_storage_s[i],
        )
        p_diff_s_s[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_S,
            dip_pore_s_s[i],
            dip_lake_s[i],
            Θ_s,
            TP_Variables.A_Sand_S,
            lake_o_storage_s[i],
        )
        p_diff_r_s[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_R,
            dip_pore_r_s[i],
            dip_lake_s[i],
            Θ_r,
            TP_Variables.A_Rock_S,
            lake_o_storage_s[i],
        )
        p_diff_p_s[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_P,
            dip_pore_p_s[i],
            dip_lake_s[i],
            Θ_p,
            TP_Variables.A_Peat_S,
            lake_o_storage_s[i],
        )

        # Calculate the total phosphorus concentration for the next time step for the North and South regions
        tp_lake_n[i + 1] = _calculate_next_total_phosphorus_concentration(
            l_ext_m[i],
            atm_dep_n[i],
            Θ_m,
            Θ_s,
            Θ_r,
            Θ_p,
            dip_pore_m_n[i],
            dip_pore_s_n[i],
            dip_pore_r_n[i],
            dip_pore_p_n[i],
            dip_lake_n[i],
            q_n2s[i],
            lake_o_a_n[i],
            tp_lake_n[i],
            lake_o_storage_n[i],
            TP_Variables.v_diff_M,
            TP_Variables.v_diff_S,
            TP_Variables.v_diff_R,
            TP_Variables.v_diff_P,
            v_settle_n[i],
            sed_resusp_m_n[i],
            sed_resusp_s_n[i],
            sed_resusp_r_n[i],
            sed_resusp_p_n[i],
        )
        
        tp_lake_s[i + 1] = _calculate_next_total_phosphorus_concentration(
            atm_dep_s[i],
            q_n2s[i],
            Θ_m,
            Θ_s,
            Θ_r,
            Θ_p,
            dip_pore_m_s[i],
            dip_pore_s_s[i],
            dip_pore_r_s[i],
            dip_pore_p_s[i],
            dip_lake_s[i],
            q_o_m[i],
            lake_o_a_s[i],
            tp_lake_s[i],
            lake_o_storage_s[i],
            TP_Variables.v_diff_M,
            TP_Variables.v_diff_S,
            TP_Variables.v_diff_R,
            TP_Variables.v_diff_P,
            v_settle_s[i],
            sed_resusp_m_s[i],
            sed_resusp_s_s[i],
            sed_resusp_r_s[i],
            sed_resusp_p_s[i],
        )
        
        # Calculate the average total phosphorus in the lake for the next time step
        tp_lake_mean[i + 1] = (tp_lake_n[i + 1] + tp_lake_s[i + 1]) / 2

        # Calculate the phosphorus loads for the Caloosahatchee river, St. Lucie river, and the southern region
        p_load_cal[i] = s77_q[i] * CUBIC_METERS_IN_CUBIC_FOOT * SECONDS_IN_HOUR * HOURS_IN_DAY * tp_lake_s[i]  # mg/d P
        p_load_stl[i] = s308_q[i] * CUBIC_METERS_IN_CUBIC_FOOT * SECONDS_IN_HOUR * HOURS_IN_DAY * tp_lake_s[i]  # mg/d P
        p_load_south[i] = tot_reg_so[i] * CUBIC_METERS_IN_ACRE_FOOT * tp_lake_s[i]  # mg/d P

    # Create the dataframes for the phosphorus loads and the total phosphorus in the lake
    p_loads_df = pd.DataFrame(date_rng_0, columns=["Date"])  # 1/1/2008-12/31/2018
    p_lake_df = pd.DataFrame(date_rng_0, columns=["Date"])  # 1/1/2008-12/31/2018

    # Add phosphorus loads as tons
    p_loads_df["P_Load_Cal"] = pd.to_numeric(p_load_cal) / MILLIGRAMS_IN_TON  # tons
    p_loads_df["P_Load_StL"] = pd.to_numeric(p_load_stl) / MILLIGRAMS_IN_TON  # tons
    p_loads_df["P_Load_South"] = pd.to_numeric(p_load_south) / MILLIGRAMS_IN_TON  # tons
    
    # Add total phosphorus (in southern region) and average phosphorus in the lake
    p_lake_df["P_Lake"] = pd.to_numeric(tp_lake_mean)
    p_lake_df["TP_Lake_S"] = pd.to_numeric(tp_lake_s)
    
    # Set the date as the index for the dataframe
    p_loads_df = p_loads_df.set_index("Date")
    p_loads_df.index = pd.to_datetime(p_loads_df.index, unit="ns")
    
    # Resample the dataframes to get monthly totals
    p_loads_m = p_loads_df.resample("M").sum()
    p_loads_m = p_loads_m.reset_index()
    
    # Set the date as the index for the dataframe
    p_lake_df = p_lake_df.set_index("Date")

    return p_lake_df

def _calculate_porosity(bulk_density, particle_density, water_content):
    """
    Calculate the porosity of a sediment.

    Args:
        bulk_density (float): The bulk density of the sediment (in grams per cubic centimeter).
        particle_density (float): The particle density of the sediment (in graps per cubic centimeter).
        water_content (float): The water content of the sediment (a percentage, 0 to 100).

    Returns:
        float: The porosity of the sediment.
    """
    return 1 - ((bulk_density / particle_density) * ((100 - water_content) / 100))

def _calculate_mass_sediment(area, thickness, water_content, bulk_density):
    """
    Calculate the mass of sediment.

    Args:
        area (float): Area of sediment (in square meters).
        thickness (float): Thickness of sediment (in meters).
        water_content (float): Percentage of water content in sediment.
        bulk_density (float): Bulk density of sediment (in grams per cubic centimeter).

    Returns:
        float: The mass of sediment (in kilograms).
    """
    return area * thickness * ((100 - water_content) / 100) * bulk_density * 1000


def _calculate_DIP_pore(
    porosity,
    initial_DIP_pore,
    DIP_lake,
    desorption_flux,
    adsorption_flux,
    total_phosphorus,
    sediment_mass,
    diffusion_coefficient,
    surface_area,
    decomposition_rate,
    burial_velocity,
):
    """
    Calculate the dissolved inorganic phosphorus (DIP) concentration in the pore water.

    Args:
        porosity: The porosity of the sediment.
        initial_DIP_pore: The initial DIP concentration in the pore water.
        DIP_lake: The DIP concentration in the lake water.
        desorption_flux: The desorption flux of DIP from the sediment to the pore water.
        adsorption_flux: The adsorption flux of DIP from the pore water to the sediment.
        total_phosphorus: The total phosphorus (TP) concentration in the sediment.
        sediment_mass: The mass of the sediment.
        diffusion_coefficient: The diffusion coefficient of DIP in the pore water.
        surface_area: The surface area of the sediment-water interface.
        decomposition_rate: The decomposition rate constant of organic matter in the sediment.
        burial_velocity: The burial velocity of the sediment.

    Returns:
        The calculated DIP concentration in the pore water. If the result is negative, it is returned as 0.
    """
    dip_concentration = TP_MBFR.DIP_pore(
        porosity,
        initial_DIP_pore,
        DIP_lake,
        desorption_flux,
        adsorption_flux,
        total_phosphorus,
        sediment_mass,
        diffusion_coefficient,
        surface_area,
        decomposition_rate,
        burial_velocity,
    )
    return dip_concentration if dip_concentration > 0 else 0


def _calculate_sediment_area(lake_overall_area, area_type, area_mud, area_sand, area_rock, area_peat):
    """
    Calculate the area of a specific type of sediment in a region of a lake.

    Args:
        lake_overall_area (float): The overall area of the lake.
        area_type (float): The area of the specific type of sediment (in square meters).
        area_mud (float): The area of mud sediment (in square meters).
        area_sand (float): The area of sand sediment (in square meters).
        area_rock (float): The area of rock sediment (in square meters).
        area_peat (float): The area of peat sediment (in square meters).

    Returns:
        float: The calculated area of the specific type of sediment in the region.
    """
    return lake_overall_area * area_type / (area_mud + area_sand + area_rock + area_peat)


def _calculate_settling_velocity(gas_constant, gravity, clay_diameter, sand_diameter, 
                                clay_drag_coefficient, sand_drag_coefficient, 
                                clay_lift_coefficient, sand_lift_coefficient, 
                                dynamic_viscosity, mud_area, peat_area, 
                                sand_area, rock_area, total_area):
    """
    Calculate the settling velocity for a region of a lake.

    Args:
        gas_constant (float): The universal gas constant.
        gravity (float): The acceleration due to gravity.
        clay_diameter (float): The diameter of the clay particles.
        sand_diameter (float): The diameter of the sand particles.
        clay_drag_coefficient (float): The drag coefficient for clay.
        sand_drag_coefficient (float): The drag coefficient for sand.
        clay_lift_coefficient (float): The lift coefficient for clay.
        sand_lift_coefficient (float): The lift coefficient for sand.
        dynamic_viscosity (float): The dynamic viscosity of the fluid.
        mud_area (float): The area of mud sediment.
        peat_area (float): The area of peat sediment.
        sand_area (float): The area of sand sediment.
        rock_area (float): The area of rock sediment.
        total_area (float): The overall area of the region.

    Returns:
        float: The settling velocity for the region.
    """
    v_settle_clay = (gas_constant * gravity * clay_diameter**2) / (clay_drag_coefficient * dynamic_viscosity + (0.75 * clay_lift_coefficient * gas_constant * gravity * clay_diameter**3) ** 0.5)
    v_settle_sand = (gas_constant * gravity * sand_diameter**2) / (sand_drag_coefficient * dynamic_viscosity + (0.75 * sand_lift_coefficient * gas_constant * gravity * sand_diameter**3) ** 0.5)
    v_settle = v_settle_clay * ((mud_area + peat_area) / total_area) + v_settle_sand * ((sand_area + rock_area) / total_area)
    return v_settle


def _calculate_sediment_resuspension(sedimentation_constant, time_delay, E_1, E_2, 
                                    wind_shear_stress, critical_shear_stress, 
                                    lake_oxygen_water_depth, sediment_proportion):
    """
    Calculate the sediment resuspension.

    Args:
        sedimentation_constant (float): The sedimentation constant (E_0).
        time_delay (float): The time delay (Td).
        E_1 (float): The first power constant.
        E_2 (float): The second power constant.
        wind_shear_stress (float): The wind shear stress value.
        critical_shear_stress (float): The critical shear stress.
        lake_oxygen_water_depth (float): The lake oxygen water depth.
        sediment_proportion (float): The sediment proportion.

    Returns:
        float: The sediment resuspension for the given parameters.
    """
    if wind_shear_stress > critical_shear_stress:
        return ((sedimentation_constant / time_delay**E_1) * ((wind_shear_stress - critical_shear_stress) / critical_shear_stress) ** E_2) * 10 / lake_oxygen_water_depth * sediment_proportion
    else:
        return 0


def _calculate_next_sediment(lake_oxygen_area, total_phosphorus_lake, dissolved_inorganic_phosphorus_lake, sediment_burial_flux, current_sediment_phosphorus, current_sediment_mass, decomposition_constant, settling_velocity, sediment_resuspension, lake_oxygen_storage):
    """
    Calculate the next sediment value based on various parameters.

    Args:
        lake_oxygen_area (float): The current lake oxygen area.
        total_phosphorus_lake (float): The total phosphorus in the lake.
        dissolved_inorganic_phosphorus_lake (float): The dissolved inorganic phosphorus in the lake.
        sediment_burial_flux (float): The sediment burial flux.
        current_sediment_phosphorus (float): The current sediment phosphorus.
        current_sediment_mass (float): The current sediment mass.
        decomposition_constant (float): The decomposition constant.
        settling_velocity (float): The settling velocity.
        sediment_resuspension (float): The sediment resuspension.
        lake_oxygen_storage (float): The current lake oxygen storage.

    Returns:
        float: The next sediment value. If the calculated value is negative, it returns 0.
    """
    sediment = TP_MBFR.P_sed(lake_oxygen_area, total_phosphorus_lake, dissolved_inorganic_phosphorus_lake, sediment_burial_flux, current_sediment_phosphorus, current_sediment_mass, decomposition_constant, settling_velocity)
    sediment -= sediment_resuspension * lake_oxygen_storage / current_sediment_mass
    return sediment if sediment > 0 else 0


def _calculate_next_concentration(adsorption_rate, desorption_rate, burial_rate, current_concentration, sediment_mass):
    """
    Calculate the concentration of phosphorus in the sediment at the next time step.

    This function uses the Sor_P_conc method of the TP_MBFR object to calculate the new concentration based on the 
    current concentration and the rates of adsorption, desorption, and burial. If the calculated new concentration is 
    greater than 0, it is returned. Otherwise, 0 is returned to ensure that the concentration does not become negative.

    Args:
        adsorption_rate (float): The rate of adsorption of phosphorus onto the sediment.
        desorption_rate (float): The rate of desorption of phosphorus from the sediment.
        burial_rate (float): The rate of burial of phosphorus in the sediment.
        current_concentration (float): The current concentration of phosphorus in the sediment.
        sediment_mass (float): The mass of the sediment.

    Returns:
        float: The concentration of phosphorus in the sediment at the next time step.
    """
    next_concentration = TP_MBFR.Sor_P_conc(adsorption_rate, desorption_rate, burial_rate, current_concentration, sediment_mass)
    return next_concentration if next_concentration > 0 else 0


def _calculate_next_total_phosphorus_concentration(
    external_loading,
    atmospheric_deposition,
    mixing_factor_M,
    mixing_factor_S,
    mixing_factor_R,
    mixing_factor_P,
    pore_water_DIP_M,
    pore_water_DIP_S,
    pore_water_DIP_R,
    pore_water_DIP_P,
    lake_DIP,
    flow_rate_N2S,
    lake_oxygen_concentration,
    current_tp_concentration,
    lake_oxygen_storage,
    vertical_diffusion_M,
    vertical_diffusion_S,
    vertical_diffusion_R,
    vertical_diffusion_P,
    vertical_settling,
    sediment_resuspension_M,
    sediment_resuspension_S,
    sediment_resuspension_R,
    sediment_resuspension_P,
):
    """
    Calculates the next total phosphorus (TP) concentration in the lake.

    Args:
        external_loading (float): External loading of phosphorus from the catchment to the lake (M).
        atmospheric_deposition (float): Atmospheric deposition of phosphorus to the lake.
        mixing_factor_M (float): Mixing factor for phosphorus in the lake (M).
        mixing_factor_S (float): Mixing factor for phosphorus in the lake (S).
        mixing_factor_R (float): Mixing factor for phosphorus in the lake (R).
        mixing_factor_P (float): Mixing factor for phosphorus in the lake (P).
        pore_water_DIP_M (float): Dissolved inorganic phosphorus (DIP) from pore water to the lake (M).
        pore_water_DIP_S (float): Dissolved inorganic phosphorus (DIP) from pore water to the lake (S).
        pore_water_DIP_R (float): Dissolved inorganic phosphorus (DIP) from pore water to the lake (R).
        pore_water_DIP_P (float): Dissolved inorganic phosphorus (DIP) from pore water to the lake (P).
        lake_DIP (float): DIP concentration in the lake.
        flow_rate_N2S (float): Flow rate from the lake to the sea.
        lake_oxygen_concentration (float): Oxygen concentration in the lake (A).
        current_tp_concentration (float): Total phosphorus (TP) concentration in the lake.
        lake_oxygen_storage (float): Oxygen storage in the lake.
        vertical_diffusion_M (float): Vertical diffusion of phosphorus in the lake (M).
        vertical_diffusion_S (float): Vertical diffusion of phosphorus in the lake (S).
        vertical_diffusion_R (float): Vertical diffusion of phosphorus in the lake (R).
        vertical_diffusion_P (float): Vertical diffusion of phosphorus in the lake (P).
        vertical_settling (float): Vertical settling of phosphorus in the lake.
        sediment_resuspension_M (float): Sediment resuspension of phosphorus in the lake (M).
        sediment_resuspension_S (float): Sediment resuspension of phosphorus in the lake (S).
        sediment_resuspension_R (float): Sediment resuspension of phosphorus in the lake (R).
        sediment_resuspension_P (float): Sediment resuspension of phosphorus in the lake (P).

    Returns:
        float: The calculated next total phosphorus (TP) concentration in the lake.
    """
    tp_lake_result = (
        TP_MBFR.TP_Lake_N(
            external_loading,
            atmospheric_deposition,
            mixing_factor_M,
            mixing_factor_S,
            mixing_factor_R,
            mixing_factor_P,
            pore_water_DIP_M,
            pore_water_DIP_S,
            pore_water_DIP_R,
            pore_water_DIP_P,
            lake_DIP,
            flow_rate_N2S,
            lake_oxygen_concentration,
            current_tp_concentration,
            lake_oxygen_storage,
            vertical_diffusion_M,
            vertical_diffusion_S,
            vertical_diffusion_R,
            vertical_diffusion_P,
            vertical_settling,
        )
        + sediment_resuspension_M
        + sediment_resuspension_S
        + sediment_resuspension_R
        + sediment_resuspension_P
    )

    return tp_lake_result if tp_lake_result > 0 else 0
