import os
import argparse
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from loone.utils import load_config, stg_sto_ar
from loone.utils import tp_mass_balance_functions_regions as TP_MBFR
from loone.data import Data as DClass
from loone.data.tp_variables_regions import TP_Variables as TPVarClass


SECONDS_IN_DAY = 86400
CUBIC_METERS_IN_ACRE_FOOT = 1233.48
SQUARE_METERS_IN_ACRE = 4046.85642
METERS_IN_FOOT = 0.3048
METERS_IN_FOOT_3 = 0.305
SECONDS_IN_HOUR = 3600
HOURS_IN_DAY = 24
CUBIC_METERS_IN_CUBIC_FOOT = 0.028316847
MILLIGRAMS_IN_TON = 1e9


def LOONE_NUT(
    workspace: str,
    out_file_name: str,
    loads_external_filename: str,
    flow_df_filename: str,
    loone_q_path: str | None = None,
    input_data: str | None = None,
    forecast_mode: bool = False,
    simulation_data: dict = {},
    ensemble: int = None,
) -> pd.DataFrame:
    """Simulates nutrient (phosphorus) dynamics in the water column.

    Args:
        workspace (str): Path to the workspace containing config file.
        out_file_name (str): Path to LOONE nutrient file to be created.
        loads_external_filename (str): The name of the file that holds the external loads for the model.
        flow_df_filename (str): The name of the file that holds the flow data for the model.
        loone_q_path (str | None, optional): Path to LOONE Q CSV file. Defaults to "<workspace>/LOONE_Q_Outputs.csv".
        input_data (str | None, optional): Path to the data directory. Defaults to <workspace>.
        forecast_mode (bool, optional): Whether to run the model in forecast mode. Defaults to False.
        simulation_data (dict, optional): Dictionary containing simulation data when sim_type is 3. Defaults to {}. Required keys are 's77_q', 's308_q', 'tot_reg_so', 'stage_lo'.

    Returns:
        pd.DataFrame: A Dataframe containing an estimate of the total phosphorus concentration in the lake for a certain time series.
    """
    os.chdir(workspace)
    config = load_config(workspace)

    TP_Variables = TPVarClass(workspace)

    if not loone_q_path:
        loone_q_path = os.path.join(workspace, "LOONE_Q_Outputs.csv")
    Data = DClass(workspace, forecast_mode, ensemble)
    print("LOONE Nut Module is Running!")

    # Read in config values

    data_dir = input_data if input_data else workspace
    loone_q = pd.read_csv(loone_q_path) #this is in cubic feet
    # Based on the defined Start and End year, month, and day on the
    # Pre_defined_Variables File, Startdate and enddate are defined.
    if forecast_mode:
        startdate = datetime.now().date()  # datetime(year, month, day).date()
        enddate = startdate + timedelta(
            days=15
        )  # datetime(year, month, day).date()
    else:
        year, month, day = map(int, config["start_date_entry"])
        startdate = datetime(year, month, day).date()
        year, month, day = map(int, config["end_date_entry"])
        enddate = datetime(year, month, day).date()

    date_rng_0 = pd.date_range(start=startdate, end=enddate, freq="D")
    load_ext = pd.read_csv(os.path.join(data_dir, loads_external_filename))
    if forecast_mode:
        q_in = pd.read_csv(os.path.join(data_dir, f"LO_Inflows_BK_forecast_{ensemble:02}.csv")) #cmd
    else:
        q_in = pd.read_csv(os.path.join(data_dir, config["lo_inflows_bk"]))
    flow_df = pd.read_csv(os.path.join(data_dir, flow_df_filename)) #cubic meters per day
    # q_o = flow_df["Outflows"].values #cmd - this is wrong, this should come from LOONE_Q, not from geoglows
    s77_q = loone_q["S77_Q"].values if 's77_q' not in simulation_data else simulation_data['s77_q'] #cfs
    s308_q = loone_q["S308_Q"].values if 's308_q' not in simulation_data else simulation_data['s308_q'] #cfs
    #tot_reg_so should be coming from LOONE_Q, not from flow_df - these units are acft/day
    if 'tot_reg_so' in simulation_data:
        tot_reg_so = simulation_data['tot_reg_so']
    #TODO check with Osama
    else:
        tot_reg_so = loone_q["TotRegSo"]
    #New way to calculate q_o  - add s77_q, s308_q, and tot_reg_so - this should be converted to cmd - convert all of them to cmd, then add them together
    q_o = s308_q * CUBIC_METERS_IN_CUBIC_FOOT * SECONDS_IN_DAY + s77_q * CUBIC_METERS_IN_CUBIC_FOOT * SECONDS_IN_DAY + tot_reg_so * CUBIC_METERS_IN_ACRE_FOOT #cmd
    #TODO: Should it read the loone_q outputs for historical data as well?
    if forecast_mode:
        sto_stage = pd.read_csv(os.path.join(data_dir, f"LOONE_Q_Outputs_{ensemble:02}.csv")) #acft
        stage_lo = sto_stage["Stage"].values if 'stage_lo' not in simulation_data else simulation_data['stage_lo'] #feet
        storage = sto_stage["Storage"].values
        # storage dev is 0 in forecast mode
        storage_dev = np.zeros(len(storage), dtype=float)  # ac-ft
    else:
        sto_stage = pd.read_csv(
            os.path.join(data_dir, config["sto_stage"])
        )
        #TODO - we should get the dates for this
        stage_lo = sto_stage["Stage_ft"].values if 'stage_lo' not in simulation_data else simulation_data['stage_lo'] #feet
        storage = sto_stage["Storage_acft"].values
        storage_dev = Data.Storage_dev_df["DS_dev"]
        start_storage = stg_sto_ar.stg2sto(config["start_stage"], 0) 
        stage_2_ar[1] = stg_sto_ar.stg2ar(stage_lo[1], 0)
        storage[0] = start_storage  # ac-ft
        storage[1] = stg_sto_ar.stg2sto(stage_lo[1], 0)  # ac-ft
    n_rows = len(q_in.index)
    l_ext = load_ext["TP_Loads_In_mg"]  # mg
    atm_dep_n = TP_Variables.northern_percentage * load_ext["Atm_Loading_mg"]
    atm_dep_s = TP_Variables.southern_percentage * load_ext["Atm_Loading_mg"]

    # Read Shear Stress driven by Wind Speed
    if forecast_mode:
        wind_shear_str = pd.read_csv(
            os.path.join(data_dir, f"WindShearStress_{ensemble:02}.csv")
        )
    else:
        wind_shear_str = pd.read_csv(
            os.path.join(data_dir, config["wind_shear_stress"])
        )
    w_ss = wind_shear_str["ShearStress"]  # Dyne/cm2
    if forecast_mode:
        nu_ts = pd.read_csv(os.path.join(data_dir, f"nu_predicted.csv"))
    else:
        nu_ts = pd.read_csv(os.path.join(data_dir, config["nu"]))
    LO_BL = 0.5  # m (Bed Elevation of LO)
    g = 9.8  # m/s2 gravitational acceleration
    cal_res = pd.read_csv(os.path.join(data_dir, "nondominated_Sol_var.csv"))
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
    # q_o = np.zeros(n_rows, dtype=object) - this should not be initialized here, it should come from LOONE_Q
    q_o_m = np.zeros(n_rows, dtype=object)

    p_load_cal = np.zeros(n_rows, dtype=object)
    p_load_stl = np.zeros(n_rows, dtype=object)
    p_load_south = np.zeros(n_rows, dtype=object)

    v_settle_n = np.zeros(n_rows, dtype=object)
    v_settle_s = np.zeros(n_rows, dtype=object)

    ##Initial Values##
    # S.A. is calculated based on the Lake's previous time step Stage, but for
    # the S.A. at i=0 I used same time step Stage!
    # TODO: Does this need to be fixed in forecast mode? - come from dbhydro to get the stage for that day
    # In forecast mode does not read from the config file - make sure this is in feet
    # start_storage = stg_sto_ar.stg2sto(config["start_stage"], 0) 
    stage_2_ar[1] = stg_sto_ar.stg2ar(stage_lo[1], 0)
    # This is not needed in this case
    # storage[0] = start_storage  # ac-ft
    # storage[1] = stg_sto_ar.stg2sto(stage_lo[1], 0)  # ac-ft
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
    Θ_m = _calculate_porosity(
        TP_Variables.bulk_density_mud,
        TP_Variables.particle_density_mud,
        TP_Variables.water_content_mud,
    )
    Θ_s = _calculate_porosity(
        TP_Variables.bulk_density_sand,
        TP_Variables.particle_density_sand,
        TP_Variables.water_content_sand,
    )
    Θ_r = _calculate_porosity(
        TP_Variables.bulk_density_rock,
        TP_Variables.particle_density_rock,
        TP_Variables.water_content_rock,
    )
    Θ_p = _calculate_porosity(
        TP_Variables.bulk_density_peat,
        TP_Variables.particle_density_peat,
        TP_Variables.water_content_peat,
    )

    # Mass of sediment in surfacial mix Mud layer in the North Region(kg)
    mass_sed_m_n = _calculate_mass_sediment(
        TP_Variables.area_mud_north,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_mud,
        TP_Variables.bulk_density_mud,
    )

    # Mass of sediment in surfacial mix Sand layer in the North Region(kg)
    mass_sed_s_n = _calculate_mass_sediment(
        TP_Variables.area_sand_north,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_sand,
        TP_Variables.bulk_density_sand,
    )

    # Mass of sediment in surfacial mix Rock layer in the North Region(kg)
    mass_sed_r_n = _calculate_mass_sediment(
        TP_Variables.area_rock_north,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_rock,
        TP_Variables.bulk_density_rock,
    )

    # Mass of sediment in surfacial mix Peat layer in the North Region(kg)
    mass_sed_p_n = _calculate_mass_sediment(
        TP_Variables.area_peat_north,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_peat,
        TP_Variables.bulk_density_peat,
    )

    # Mass of sediment in surfacial mix Mud layer in the South Region(kg)
    mass_sed_m_s = _calculate_mass_sediment(
        TP_Variables.area_mud_south,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_mud,
        TP_Variables.bulk_density_mud,
    )

    # Mass of sediment in surfacial mix Sand layer in the South Region(kg)
    mass_sed_s_s = _calculate_mass_sediment(
        TP_Variables.area_sand_south,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_sand,
        TP_Variables.bulk_density_sand,
    )

    # Mass of sediment in surfacial mix Rock layer in the South Region(kg)
    mass_sed_r_s = _calculate_mass_sediment(
        TP_Variables.area_rock_south,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_rock,
        TP_Variables.bulk_density_rock,
    )

    # Mass of sediment in surfacial mix Peat layer in the South Region(kg)
    mass_sed_p_s = _calculate_mass_sediment(
        TP_Variables.area_peat_south,
        TP_Variables.sediment_depth,
        TP_Variables.water_content_peat,
        TP_Variables.bulk_density_peat,
    )

    for i in range(n_rows - 2):
        if forecast_mode:
            storage_dev[i] = 0  # in forecast mode storage_dev is 0
        if storage_dev[i] >= 0:
            q_i_m[i] = (
                q_i[i] + storage_dev[i] * CUBIC_METERS_IN_ACRE_FOOT #in forecast mode storage_dev is 0
            )  # m3/d
            
            q_o_m[i] = q_o[i]
            l_ext_m[i] = l_ext[i] + q_i_m[i] * tp_lake_n[i]
        else:
            q_o_m[i] = (
                q_o[i] - storage_dev[i] * CUBIC_METERS_IN_ACRE_FOOT
            )  # m3/d
            q_i_m[i] = q_i[i]
            l_ext_m[i] = l_ext[i]
        q_n2s[i] = (q_i_m[i] + q_o_m[i]) / 2
        stage_2_ar[i + 2] = stg_sto_ar.stg2ar(stage_lo[i + 2], 0)
        lo_wd[i] = stage_lo[i] * METERS_IN_FOOT - LO_BL
        # TODO - are these the right values for volume in the lake
        lake_o_storage_n[i] = (
            storage[i] #reads in acft
            * TP_Variables.northern_percentage
            * SQUARE_METERS_IN_ACRE
            * METERS_IN_FOOT_3
        )  # m3
        lake_o_storage_s[i] = (
            storage[i]  # reads in acft
            * TP_Variables.southern_percentage
            * SQUARE_METERS_IN_ACRE
            * METERS_IN_FOOT_3
        )  # m3
        lake_o_a_n[i] = (
            stage_2_ar[i] * TP_Variables.northern_percentage * SQUARE_METERS_IN_ACRE
        )  # m2
        lake_o_a_s[i] = (
            stage_2_ar[i] * TP_Variables.southern_percentage * SQUARE_METERS_IN_ACRE
        )  # m2
        lake_o_a_m_n[i] = _calculate_sediment_area(
            lake_o_a_n[i],
            TP_Variables.area_mud_north,
            TP_Variables.area_mud_north,
            TP_Variables.area_sand_north,
            TP_Variables.area_rock_north,
            TP_Variables.area_peat_north,
        )
        lake_o_a_s_n[i] = _calculate_sediment_area(
            lake_o_a_n[i],
            TP_Variables.area_sand_north,
            TP_Variables.area_mud_north,
            TP_Variables.area_sand_north,
            TP_Variables.area_rock_north,
            TP_Variables.area_peat_north,
        )
        lake_o_a_r_n[i] = _calculate_sediment_area(
            lake_o_a_n[i],
            TP_Variables.area_rock_north,
            TP_Variables.area_mud_north,
            TP_Variables.area_sand_north,
            TP_Variables.area_rock_north,
            TP_Variables.area_peat_north,
        )
        lake_o_a_p_n[i] = _calculate_sediment_area(
            lake_o_a_n[i],
            TP_Variables.area_peat_north,
            TP_Variables.area_mud_north,
            TP_Variables.area_sand_north,
            TP_Variables.area_rock_north,
            TP_Variables.area_peat_north,
        )
        lake_o_a_m_s[i] = _calculate_sediment_area(
            lake_o_a_s[i],
            TP_Variables.area_mud_south,
            TP_Variables.area_mud_south,
            TP_Variables.area_sand_south,
            TP_Variables.area_rock_south,
            TP_Variables.area_peat_south,
        )
        lake_o_a_s_s[i] = _calculate_sediment_area(
            lake_o_a_s[i],
            TP_Variables.area_sand_south,
            TP_Variables.area_mud_south,
            TP_Variables.area_sand_south,
            TP_Variables.area_rock_south,
            TP_Variables.area_peat_south,
        )
        lake_o_a_r_s[i] = _calculate_sediment_area(
            lake_o_a_s[i],
            TP_Variables.area_rock_south,
            TP_Variables.area_mud_south,
            TP_Variables.area_sand_south,
            TP_Variables.area_rock_south,
            TP_Variables.area_peat_south,
        )
        lake_o_a_p_s[i] = _calculate_sediment_area(
            lake_o_a_s[i],
            TP_Variables.area_peat_south,
            TP_Variables.area_mud_south,
            TP_Variables.area_sand_south,
            TP_Variables.area_rock_south,
            TP_Variables.area_peat_south,
        )

        dip_lake_n[i] = TP_MBFR.DIP_Lake(tp_lake_n[i])
        dip_lake_s[i] = TP_MBFR.DIP_Lake(tp_lake_s[i])

        # Calculate the settling velocity for the North and South regions
        v_settle_n[i] = _calculate_settling_velocity(
            R,
            g,
            d_c,
            d_s,
            c_1_c,
            c_1_s,
            c_2_c,
            c_2_s,
            nu_d[i],
            TP_Variables.area_mud_north,
            TP_Variables.area_peat_north,
            TP_Variables.area_sand_north,
            TP_Variables.area_rock_north,
            TP_Variables.total_area_north,
        )

        v_settle_s[i] = _calculate_settling_velocity(
            R,
            g,
            d_c,
            d_s,
            c_1_c,
            c_1_s,
            c_2_c,
            c_2_s,
            nu_d[i],
            TP_Variables.area_mud_south,
            TP_Variables.area_peat_south,
            TP_Variables.area_sand_south,
            TP_Variables.area_rock_south,
            TP_Variables.total_area_south,
        )

        # Calculate desorption fluxes for the North and South regions (for mud, sand, rock, and peat)
        j_des_m_n[i] = TP_MBFR.Des_flux(
            Γ_m_n[i], mass_sed_m_n, TP_Variables.desorption_rate_mud
        )
        j_des_s_n[i] = TP_MBFR.Des_flux(
            Γ_s_n[i], mass_sed_s_n, TP_Variables.desorption_rate_sand
        )
        j_des_r_n[i] = TP_MBFR.Des_flux(
            Γ_r_n[i], mass_sed_r_n, TP_Variables.desorption_rate_rock
        )
        j_des_p_n[i] = TP_MBFR.Des_flux(
            Γ_p_n[i], mass_sed_p_n, TP_Variables.desorption_rate_peat
        )
        j_des_m_s[i] = TP_MBFR.Des_flux(
            Γ_m_s[i], mass_sed_m_s, TP_Variables.desorption_rate_mud
        )
        j_des_s_s[i] = TP_MBFR.Des_flux(
            Γ_s_s[i], mass_sed_s_s, TP_Variables.desorption_rate_sand
        )
        j_des_r_s[i] = TP_MBFR.Des_flux(
            Γ_r_s[i], mass_sed_r_s, TP_Variables.desorption_rate_rock
        )
        j_des_p_s[i] = TP_MBFR.Des_flux(
            Γ_p_s[i], mass_sed_p_s, TP_Variables.desorption_rate_peat
        )

        # Calculate adsorption fluxes for the North and South regions
        j_ads_m_n[i] = TP_MBFR.Ads_flux(
            dip_pore_m_n[i],
            Γ_m_n[i],
            mass_sed_m_n,
            TP_Variables.adsorption_rate_mud,
            TP_Variables.inorganic_fraction,
        )
        j_ads_s_n[i] = TP_MBFR.Ads_flux(
            dip_pore_s_n[i],
            Γ_s_n[i],
            mass_sed_s_n,
            TP_Variables.adsorption_rate_sand,
            TP_Variables.inorganic_fraction,
        )
        j_ads_r_n[i] = TP_MBFR.Ads_flux(
            dip_pore_r_n[i],
            Γ_r_n[i],
            mass_sed_r_n,
            TP_Variables.adsorption_rate_rock,
            TP_Variables.inorganic_fraction,
        )
        j_ads_p_n[i] = TP_MBFR.Ads_flux(
            dip_pore_p_n[i],
            Γ_p_n[i],
            mass_sed_p_n,
            TP_Variables.adsorption_rate_peat,
            TP_Variables.inorganic_fraction,
        )
        j_ads_m_s[i] = TP_MBFR.Ads_flux(
            dip_pore_m_s[i],
            Γ_m_s[i],
            mass_sed_m_s,
            TP_Variables.adsorption_rate_mud,
            TP_Variables.inorganic_fraction,
        )
        j_ads_s_s[i] = TP_MBFR.Ads_flux(
            dip_pore_s_s[i],
            Γ_s_s[i],
            mass_sed_s_s,
            TP_Variables.adsorption_rate_sand,
            TP_Variables.inorganic_fraction,
        )
        j_ads_r_s[i] = TP_MBFR.Ads_flux(
            dip_pore_r_s[i],
            Γ_r_s[i],
            mass_sed_r_s,
            TP_Variables.adsorption_rate_rock,
            TP_Variables.inorganic_fraction,
        )
        j_ads_p_s[i] = TP_MBFR.Ads_flux(
            dip_pore_p_s[i],
            Γ_p_s[i],
            mass_sed_p_s,
            TP_Variables.adsorption_rate_peat,
            TP_Variables.inorganic_fraction,
        )

        # Calculate the sediment burial fluxes for the North and South regions (for mud, sand, rock, and peat)
        j_sedburial_m_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_m_n[i],
            TP_Variables.bulk_density_mud,
            TP_Variables.area_mud_north,
            TP_Variables.burial_velocity_mud,
            TP_Variables.water_content_mud,
        )
        j_sedburial_s_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_s_n[i],
            TP_Variables.bulk_density_sand,
            TP_Variables.area_sand_north,
            TP_Variables.burial_velocity_sand,
            TP_Variables.water_content_sand,
        )
        j_sedburial_r_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_r_n[i],
            TP_Variables.bulk_density_rock,
            TP_Variables.area_rock_north,
            TP_Variables.burial_velocity_rock,
            TP_Variables.water_content_rock,
        )
        j_sedburial_p_n[i] = TP_MBFR.Sed_burial_flux(
            p_sed_p_n[i],
            TP_Variables.bulk_density_peat,
            TP_Variables.area_peat_north,
            TP_Variables.burial_velocity_peat,
            TP_Variables.water_content_peat,
        )
        j_sedburial_m_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_m_s[i],
            TP_Variables.bulk_density_mud,
            TP_Variables.area_mud_south,
            TP_Variables.burial_velocity_mud,
            TP_Variables.water_content_mud,
        )
        j_sedburial_s_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_s_s[i],
            TP_Variables.bulk_density_sand,
            TP_Variables.area_sand_south,
            TP_Variables.burial_velocity_sand,
            TP_Variables.water_content_sand,
        )
        j_sedburial_r_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_r_s[i],
            TP_Variables.bulk_density_rock,
            TP_Variables.area_rock_south,
            TP_Variables.burial_velocity_rock,
            TP_Variables.water_content_rock,
        )
        j_sedburial_p_s[i] = TP_MBFR.Sed_burial_flux(
            p_sed_p_s[i],
            TP_Variables.bulk_density_peat,
            TP_Variables.area_peat_south,
            TP_Variables.burial_velocity_peat,
            TP_Variables.water_content_peat,
        )

        # Calculate the sediment resuspension for the North and South regions (for mud, sand, rock, and peat)
        sed_resusp_m_n[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_m_n[i]
        )
        sed_resusp_s_n[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_s_n[i]
        )
        sed_resusp_r_n[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_r_n[i]
        )
        sed_resusp_p_n[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_p_n[i]
        )
        sed_resusp_m_s[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_m_s[i]
        ) #mg
        sed_resusp_s_s[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_s_s[i]
        )
        sed_resusp_r_s[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_r_s[i]
        )
        sed_resusp_p_s[i] = _calculate_sediment_resuspension(
            E_0, td, E_1, E_2, w_ss[i], crtcl_sh_str, lo_wd[i], p_sed_p_s[i]
        )

        # Calculate the phosphorus concentration for the next time step for the North and South regions (for mud, sand, rock, and peat)
        p_sed_m_n[i + 1] = _calculate_next_sediment(
            lake_o_a_m_n[i],
            tp_lake_n[i],
            dip_lake_n[i],
            j_sedburial_m_n[i],
            p_sed_m_n[i],
            mass_sed_m_n,
            TP_Variables.decomposition_rate_mud,
            v_settle_n[i],
            sed_resusp_m_n[i],
            lake_o_storage_n[i],
        ) #units are correct. It is okay staying as just mg
        p_sed_s_n[i + 1] = _calculate_next_sediment(
            lake_o_a_s_n[i],
            tp_lake_n[i],
            dip_lake_n[i],
            j_sedburial_s_n[i],
            p_sed_s_n[i],
            mass_sed_s_n,
            TP_Variables.decomposition_rate_sand,
            v_settle_n[i],
            sed_resusp_s_n[i],
            lake_o_storage_n[i],
        )
        p_sed_r_n[i + 1] = _calculate_next_sediment(
            lake_o_a_r_n[i],
            tp_lake_n[i],
            dip_lake_n[i],
            j_sedburial_r_n[i],
            p_sed_r_n[i],
            mass_sed_r_n,
            TP_Variables.decomposition_rate_rock,
            v_settle_n[i],
            sed_resusp_r_n[i],
            lake_o_storage_n[i],
        )
        p_sed_p_n[i + 1] = _calculate_next_sediment(
            lake_o_a_p_n[i],
            tp_lake_n[i],
            dip_lake_n[i],
            j_sedburial_p_n[i],
            p_sed_p_n[i],
            mass_sed_p_n,
            TP_Variables.decomposition_rate_peat,
            v_settle_n[i],
            sed_resusp_p_n[i],
            lake_o_storage_n[i],
        )
        p_sed_m_s[i + 1] = _calculate_next_sediment(
            lake_o_a_m_s[i],
            tp_lake_s[i],
            dip_lake_s[i],
            j_sedburial_m_s[i],
            p_sed_m_s[i],
            mass_sed_m_s,
            TP_Variables.decomposition_rate_mud,
            v_settle_s[i],
            sed_resusp_m_s[i],
            lake_o_storage_s[i],
        )
        p_sed_s_s[i + 1] = _calculate_next_sediment(
            lake_o_a_s_s[i],
            tp_lake_s[i],
            dip_lake_s[i],
            j_sedburial_s_s[i],
            p_sed_s_s[i],
            mass_sed_s_s,
            TP_Variables.decomposition_rate_sand,
            v_settle_s[i],
            sed_resusp_s_s[i],
            lake_o_storage_s[i],
        )
        p_sed_r_s[i + 1] = _calculate_next_sediment(
            lake_o_a_r_s[i],
            tp_lake_s[i],
            dip_lake_s[i],
            j_sedburial_r_s[i],
            p_sed_r_s[i],
            mass_sed_r_s,
            TP_Variables.decomposition_rate_rock,
            v_settle_s[i],
            sed_resusp_r_s[i],
            lake_o_storage_s[i],
        )
        p_sed_p_s[i + 1] = _calculate_next_sediment(
            lake_o_a_p_s[i],
            tp_lake_s[i],
            dip_lake_s[i],
            j_sedburial_p_s[i],
            p_sed_p_s[i],
            mass_sed_p_s,
            TP_Variables.decomposition_rate_peat,
            v_settle_s[i],
            sed_resusp_p_s[i],
            lake_o_storage_s[i],
        )

        # Calculate the burial fluxes for the North and South regions (for mud, sand, rock, and peat)
        j_Γburial_m_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_m_n[i],
            TP_Variables.bulk_density_mud,
            TP_Variables.area_mud_north,
            TP_Variables.burial_velocity_mud,
            TP_Variables.water_content_mud,
        )
        j_Γburial_s_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_s_n[i],
            TP_Variables.bulk_density_sand,
            TP_Variables.area_sand_north,
            TP_Variables.burial_velocity_sand,
            TP_Variables.water_content_sand,
        )
        j_Γburial_r_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_r_n[i],
            TP_Variables.bulk_density_rock,
            TP_Variables.area_rock_north,
            TP_Variables.burial_velocity_rock,
            TP_Variables.water_content_rock,
        )
        j_Γburial_p_n[i] = TP_MBFR.Sor_P_burialflux(
            Γ_p_n[i],
            TP_Variables.bulk_density_peat,
            TP_Variables.area_peat_north,
            TP_Variables.burial_velocity_peat,
            TP_Variables.water_content_peat,
        )
        j_Γburial_m_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_m_s[i],
            TP_Variables.bulk_density_mud,
            TP_Variables.area_mud_south,
            TP_Variables.burial_velocity_mud,
            TP_Variables.water_content_mud,
        )
        j_Γburial_s_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_s_s[i],
            TP_Variables.bulk_density_sand,
            TP_Variables.area_sand_south,
            TP_Variables.burial_velocity_sand,
            TP_Variables.water_content_sand,
        )
        j_Γburial_r_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_r_s[i],
            TP_Variables.bulk_density_rock,
            TP_Variables.area_rock_south,
            TP_Variables.burial_velocity_rock,
            TP_Variables.water_content_rock,
        )
        j_Γburial_p_s[i] = TP_MBFR.Sor_P_burialflux(
            Γ_p_s[i],
            TP_Variables.bulk_density_peat,
            TP_Variables.area_peat_south,
            TP_Variables.burial_velocity_peat,
            TP_Variables.water_content_peat,
        )

        # Calculate the phosphorus concentration for the next time step for the North and South regions (for mud, sand, rock, and peat)
        Γ_m_n[i + 1] = _calculate_next_concentration(
            j_ads_m_n[i],
            j_des_m_n[i],
            j_Γburial_m_n[i],
            Γ_m_n[i],
            mass_sed_m_n,
        )
        Γ_s_n[i + 1] = _calculate_next_concentration(
            j_ads_s_n[i],
            j_des_s_n[i],
            j_Γburial_s_n[i],
            Γ_s_n[i],
            mass_sed_s_n,
        )
        Γ_r_n[i + 1] = _calculate_next_concentration(
            j_ads_r_n[i],
            j_des_r_n[i],
            j_Γburial_r_n[i],
            Γ_r_n[i],
            mass_sed_r_n,
        )
        Γ_p_n[i + 1] = _calculate_next_concentration(
            j_ads_p_n[i],
            j_des_p_n[i],
            j_Γburial_p_n[i],
            Γ_p_n[i],
            mass_sed_p_n,
        )
        Γ_m_s[i + 1] = _calculate_next_concentration(
            j_ads_m_s[i],
            j_des_m_s[i],
            j_Γburial_m_s[i],
            Γ_m_s[i],
            mass_sed_m_s,
        )
        Γ_s_s[i + 1] = _calculate_next_concentration(
            j_ads_s_s[i],
            j_des_s_s[i],
            j_Γburial_s_s[i],
            Γ_s_s[i],
            mass_sed_s_s,
        )
        Γ_r_s[i + 1] = _calculate_next_concentration(
            j_ads_r_s[i],
            j_des_r_s[i],
            j_Γburial_r_s[i],
            Γ_r_s[i],
            mass_sed_r_s,
        )
        Γ_p_s[i + 1] = _calculate_next_concentration(
            j_ads_p_s[i],
            j_des_p_s[i],
            j_Γburial_p_s[i],
            Γ_p_s[i],
            mass_sed_p_s,
        )

        # Calculate the decomposition rates for the North and South regions (for mud, sand, rock, and peat)
        j_decomp_m_n[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_mud, p_sed_m_n[i], mass_sed_m_n
        )
        j_decomp_s_n[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_sand, p_sed_s_n[i], mass_sed_s_n
        )
        j_decomp_r_n[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_rock, p_sed_r_n[i], mass_sed_r_n
        )
        j_decomp_p_n[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_peat, p_sed_p_n[i], mass_sed_p_n
        )
        j_decomp_m_s[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_mud, p_sed_m_s[i], mass_sed_m_s
        )
        j_decomp_s_s[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_sand, p_sed_s_s[i], mass_sed_s_s
        )
        j_decomp_r_s[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_rock, p_sed_r_s[i], mass_sed_r_s
        )
        j_decomp_p_s[i] = TP_MBFR.J_decomp(
            TP_Variables.decomposition_rate_peat, p_sed_p_s[i], mass_sed_p_s
        )

        # Calculate the DIP pore concentration for the North and South regions (for mud, sand, rock, and peat)
        dip_pore_m_n[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_m,
            dip_pore_m_n[i],
            dip_lake_n[i],
            j_des_m_n[i],
            j_ads_m_n[i],
            p_sed_m_n[i],
            mass_sed_m_n,
            TP_Variables.diffusion_velocity_mud,
            TP_Variables.area_mud_north,
            TP_Variables.decomposition_rate_mud,
            TP_Variables.burial_velocity_mud,
        )
        dip_pore_s_n[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_s,
            dip_pore_s_n[i],
            dip_lake_n[i],
            j_des_s_n[i],
            j_ads_s_n[i],
            p_sed_s_n[i],
            mass_sed_s_n,
            TP_Variables.diffusion_velocity_sand,
            TP_Variables.area_sand_north,
            TP_Variables.decomposition_rate_sand,
            TP_Variables.burial_velocity_sand,
        )
        dip_pore_r_n[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_r,
            dip_pore_r_n[i],
            dip_lake_n[i],
            j_des_r_n[i],
            j_ads_r_n[i],
            p_sed_r_n[i],
            mass_sed_r_n,
            TP_Variables.diffusion_velocity_rock,
            TP_Variables.area_rock_north,
            TP_Variables.decomposition_rate_rock,
            TP_Variables.burial_velocity_rock,
        )
        dip_pore_p_n[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_p,
            dip_pore_p_n[i],
            dip_lake_n[i],
            j_des_p_n[i],
            j_ads_p_n[i],
            p_sed_p_n[i],
            mass_sed_p_n,
            TP_Variables.diffusion_velocity_peat,
            TP_Variables.area_peat_north,
            TP_Variables.decomposition_rate_peat,
            TP_Variables.burial_velocity_peat,
        )
        dip_pore_m_s[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_m,
            dip_pore_m_s[i],
            dip_lake_s[i],
            j_des_m_s[i],
            j_ads_m_s[i],
            p_sed_m_s[i],
            mass_sed_m_s,
            TP_Variables.diffusion_velocity_mud,
            TP_Variables.area_mud_south,
            TP_Variables.decomposition_rate_mud,
            TP_Variables.burial_velocity_mud,
        )
        dip_pore_s_s[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_s,
            dip_pore_s_s[i],
            dip_lake_s[i],
            j_des_s_s[i],
            j_ads_s_s[i],
            p_sed_s_s[i],
            mass_sed_s_s,
            TP_Variables.diffusion_velocity_sand,
            TP_Variables.area_sand_south,
            TP_Variables.decomposition_rate_sand,
            TP_Variables.burial_velocity_sand,
        )
        dip_pore_r_s[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_r,
            dip_pore_r_s[i],
            dip_lake_s[i],
            j_des_r_s[i],
            j_ads_r_s[i],
            p_sed_r_s[i],
            mass_sed_r_s,
            TP_Variables.diffusion_velocity_rock,
            TP_Variables.area_rock_south,
            TP_Variables.decomposition_rate_rock,
            TP_Variables.burial_velocity_rock,
        )
        dip_pore_p_s[i + 1] = _calculate_DIP_pore(
            workspace,
            Θ_p,
            dip_pore_p_s[i],
            dip_lake_s[i],
            j_des_p_s[i],
            j_ads_p_s[i],
            p_sed_p_s[i],
            mass_sed_p_s,
            TP_Variables.diffusion_velocity_peat,
            TP_Variables.area_peat_south,
            TP_Variables.decomposition_rate_peat,
            TP_Variables.burial_velocity_peat,
        )

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
            TP_Variables.diffusion_velocity_mud,
            dip_pore_m_n[i],
            dip_lake_n[i],
            Θ_m,
            TP_Variables.area_mud_north,
            lake_o_storage_n[i],
        )
        p_diff_s_n[i] = TP_MBFR.Diff_P(
            TP_Variables.diffusion_velocity_sand,
            dip_pore_s_n[i],
            dip_lake_n[i],
            Θ_s,
            TP_Variables.area_sand_north,
            lake_o_storage_n[i],
        )
        p_diff_r_n[i] = TP_MBFR.Diff_P(
            TP_Variables.diffusion_velocity_rock,
            dip_pore_r_n[i],
            dip_lake_n[i],
            Θ_r,
            TP_Variables.area_rock_north,
            lake_o_storage_n[i],
        )
        p_diff_p_n[i] = TP_MBFR.Diff_P(
            TP_Variables.diffusion_velocity_peat,
            dip_pore_p_n[i],
            dip_lake_n[i],
            Θ_p,
            TP_Variables.area_peat_north,
            lake_o_storage_n[i],
        )
        p_diff_m_s[i] = TP_MBFR.Diff_P(
            TP_Variables.diffusion_velocity_mud,
            dip_pore_m_s[i],
            dip_lake_s[i],
            Θ_m,
            TP_Variables.area_mud_south,
            lake_o_storage_s[i],
        )
        p_diff_s_s[i] = TP_MBFR.Diff_P(
            TP_Variables.diffusion_velocity_sand,
            dip_pore_s_s[i],
            dip_lake_s[i],
            Θ_s,
            TP_Variables.area_sand_south,
            lake_o_storage_s[i],
        )
        p_diff_r_s[i] = TP_MBFR.Diff_P(
            TP_Variables.diffusion_velocity_rock,
            dip_pore_r_s[i],
            dip_lake_s[i],
            Θ_r,
            TP_Variables.area_rock_south,
            lake_o_storage_s[i],
        )
        p_diff_p_s[i] = TP_MBFR.Diff_P(
            TP_Variables.diffusion_velocity_peat,
            dip_pore_p_s[i],
            dip_lake_s[i],
            Θ_p,
            TP_Variables.area_peat_south,
            lake_o_storage_s[i],
        )

        # Calculate the total phosphorus concentration for the next time step for the North and South regions
        tp_lake_n[i + 1] = _calculate_next_total_phosphorus_concentration(
            workspace,
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
            TP_Variables.diffusion_velocity_mud,
            TP_Variables.diffusion_velocity_sand,
            TP_Variables.diffusion_velocity_rock,
            TP_Variables.diffusion_velocity_peat,
            v_settle_n[i],
            sed_resusp_m_n[i]/lake_o_storage_n[i],
            sed_resusp_s_n[i]/lake_o_storage_n[i],
            sed_resusp_r_n[i]/lake_o_storage_n[i],
            sed_resusp_p_n[i]/lake_o_storage_n[i],
        )

        tp_lake_s[i + 1] = _calculate_next_total_phosphorus_concentration(
            workspace,
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
            TP_Variables.diffusion_velocity_mud,
            TP_Variables.diffusion_velocity_sand,
            TP_Variables.diffusion_velocity_rock,
            TP_Variables.diffusion_velocity_peat,
            v_settle_s[i],
            sed_resusp_m_s[i]/lake_o_storage_s[i], #mg - in the equation make sure this gets converted to mg/m3
            sed_resusp_s_s[i]/lake_o_storage_s[i],
            sed_resusp_r_s[i]/lake_o_storage_s[i],
            sed_resusp_p_s[i]/lake_o_storage_s[i],
        )

        # Calculate the average total phosphorus in the lake for the next time step
        tp_lake_mean[i + 1] = (tp_lake_n[i + 1] + tp_lake_s[i + 1]) / 2

        # Calculate the phosphorus loads for the Caloosahatchee river, St. Lucie river, and the southern region
        p_load_cal[i] = (
            s77_q[i]
            * CUBIC_METERS_IN_CUBIC_FOOT
            * SECONDS_IN_HOUR
            * HOURS_IN_DAY
            * tp_lake_s[i]
        )  # mg/d P
        p_load_stl[i] = (
            s308_q[i]
            * CUBIC_METERS_IN_CUBIC_FOOT
            * SECONDS_IN_HOUR
            * HOURS_IN_DAY
            * tp_lake_s[i]
        )  # mg/d P
        p_load_south[i] = (
            tot_reg_so[i] * CUBIC_METERS_IN_ACRE_FOOT * tp_lake_s[i]
        )  # mg/d P

    # Create the dataframes for the phosphorus loads and the total phosphorus in the lake
    p_loads_df = pd.DataFrame(
        date_rng_0, columns=["Date"]
    )  # 1/1/2008-12/31/2018
    p_lake_df = pd.DataFrame(
        date_rng_0, columns=["Date"]
    )  # 1/1/2008-12/31/2018

    # Add phosphorus loads as tons
    p_loads_df["P_Load_Cal"] = (
        pd.to_numeric(p_load_cal) / MILLIGRAMS_IN_TON
    )  # tons
    p_loads_df["P_Load_StL"] = (
        pd.to_numeric(p_load_stl) / MILLIGRAMS_IN_TON
    )  # tons
    p_loads_df["P_Load_South"] = (
        pd.to_numeric(p_load_south) / MILLIGRAMS_IN_TON
    )  # tons

    # Add total phosphorus (in north and south regions) and average phosphorus in the lake
    p_lake_df["P_Lake"] = pd.to_numeric(tp_lake_mean)
    p_lake_df["TP_Lake_N"] = pd.to_numeric(tp_lake_n)
    p_lake_df["TP_Lake_S"] = pd.to_numeric(tp_lake_s)

    # Set the date as the index for the dataframe
    p_loads_df = p_loads_df.set_index("Date")
    p_loads_df.index = pd.to_datetime(p_loads_df.index, unit="ns")

    # Resample the dataframes to get monthly totals
    p_loads_m = p_loads_df.resample("ME").sum()
    p_loads_m = p_loads_m.reset_index()

    # Set the date as the index for the dataframe
    p_lake_df = p_lake_df.set_index("Date")
    p_lake_df.to_csv(out_file_name)

    if len(simulation_data) > 0:
        return [p_lake_df, p_loads_df]
    
    return p_lake_df


def _calculate_porosity(bulk_density: float, particle_density: float, water_content: float) -> float:
    """
    Calculate the porosity of a sediment.

    Args:
        bulk_density (float): The bulk density of the sediment (in grams per cubic centimeter).
        particle_density (float): The particle density of the sediment (in graps per cubic centimeter).
        water_content (float): The water content of the sediment (a percentage, 0 to 100).

    Returns:
        float: The porosity of the sediment.
    """
    return 1 - (
        (bulk_density / particle_density) * ((100 - water_content) / 100)
    )


def _calculate_mass_sediment(area: float, thickness: float, water_content: float, bulk_density: float) -> float:
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
    return (
        area * thickness * ((100 - water_content) / 100) * bulk_density * 1000
    )


def _calculate_DIP_pore(
    workspace,
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
        workspace,
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


def _calculate_sediment_area(
    lake_overall_area: float, area_type: float, area_mud: float, area_sand: float, area_rock: float, area_peat: float
) -> float:
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
    return (
        lake_overall_area
        * area_type
        / (area_mud + area_sand + area_rock + area_peat)
    )


def _calculate_settling_velocity(
    gas_constant: float,
    gravity: float,
    clay_diameter: float,
    sand_diameter: float,
    clay_drag_coefficient: float,
    sand_drag_coefficient: float,
    clay_lift_coefficient: float,
    sand_lift_coefficient: float,
    dynamic_viscosity: float,
    mud_area: float,
    peat_area: float,
    sand_area: float,
    rock_area: float,
    total_area: float,
) -> float:
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
    v_settle_clay = (gas_constant * gravity * clay_diameter**2) / (
        clay_drag_coefficient * dynamic_viscosity
        + (
            0.75
            * clay_lift_coefficient
            * gas_constant
            * gravity
            * clay_diameter**3
        )
        ** 0.5
    )
    v_settle_sand = (gas_constant * gravity * sand_diameter**2) / (
        sand_drag_coefficient * dynamic_viscosity
        + (
            0.75
            * sand_lift_coefficient
            * gas_constant
            * gravity
            * sand_diameter**3
        )
        ** 0.5
    )
    v_settle = v_settle_clay * (
        (mud_area + peat_area) / total_area
    ) + v_settle_sand * ((sand_area + rock_area) / total_area)
    return v_settle


def _calculate_sediment_resuspension(
    sedimentation_constant: float,
    time_delay: float,
    E_1: float,
    E_2: float,
    wind_shear_stress: float,
    critical_shear_stress: float,
    lake_oxygen_water_depth: float,
    sediment_proportion: float,
) -> float:
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
        return (
            (
                (sedimentation_constant / time_delay**E_1)
                * (
                    (wind_shear_stress - critical_shear_stress)
                    / critical_shear_stress
                )
                ** E_2
            )
            * 10
            / lake_oxygen_water_depth
            * sediment_proportion
        )
    else:
        return 0


def _calculate_next_sediment(
    lake_oxygen_area: float,
    total_phosphorus_lake: float,
    dissolved_inorganic_phosphorus_lake: float,
    sediment_burial_flux: float,
    current_sediment_phosphorus: float,
    current_sediment_mass: float,
    decomposition_constant: float,
    settling_velocity: float,
    sediment_resuspension: float,
    lake_oxygen_storage: float, #we do not need this value
) -> float:
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
        lake_oxygen_storage (float): The current lake oxygen storage. We do not need this here because we do not need to multiply by the volume

    Returns:
        float: The next sediment value. If the calculated value is negative, it returns 0.
    """
    sediment = TP_MBFR.P_sed(
        lake_oxygen_area,
        total_phosphorus_lake,
        dissolved_inorganic_phosphorus_lake,
        sediment_burial_flux,
        current_sediment_phosphorus,
        current_sediment_mass,
        decomposition_constant,
        settling_velocity,
    )
    sediment -= (
        # sediment_resuspension * lake_oxygen_storage / current_sediment_mass
        sediment_resuspension / current_sediment_mass
    )
    return sediment if sediment > 0 else 0


def _calculate_next_concentration(
    adsorption_rate: float,
    desorption_rate: float,
    burial_rate: float,
    current_concentration: float,
    sediment_mass: float,
) -> float:
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
    next_concentration = TP_MBFR.Sor_P_conc(
        adsorption_rate,
        desorption_rate,
        burial_rate,
        current_concentration,
        sediment_mass,
    )
    return next_concentration if next_concentration > 0 else 0


def _calculate_next_total_phosphorus_concentration(
    workspace: str,
    external_loading: float,
    atmospheric_deposition: float,
    mixing_factor_M: float,
    mixing_factor_S: float,
    mixing_factor_R: float,
    mixing_factor_P: float,
    pore_water_DIP_M: float,
    pore_water_DIP_S: float,
    pore_water_DIP_R: float,
    pore_water_DIP_P: float,
    lake_DIP: float,
    flow_rate_N2S: float,
    lake_oxygen_concentration: float,
    current_tp_concentration: float,
    lake_oxygen_storage: float,
    vertical_diffusion_M: float,
    vertical_diffusion_S: float,
    vertical_diffusion_R: float,
    vertical_diffusion_P: float,
    vertical_settling: float,
    sediment_resuspension_M: float,
    sediment_resuspension_S: float,
    sediment_resuspension_R: float,
    sediment_resuspension_P: float,
) -> float:
    """
    Calculates the next total phosphorus (TP) concentration in the lake.

    Args:
        workspace (str): The path to the working directory.
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
            workspace,
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
        + sediment_resuspension_M # these need to be divided by the volume of the lake in cubic meters
        + sediment_resuspension_S # these need to be divided by the volume of the lake in cubic meters
        + sediment_resuspension_R # these need to be divided by the volume of the lake in cubic meters
        + sediment_resuspension_P # these need to be divided by the volume of the lake in cubic meters
    )

    return tp_lake_result if tp_lake_result > 0 else 0


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "workspace",
        nargs=1,
        help="The path to the working directory.",
    )
    argparser.add_argument(
        "out_file_name",
        nargs=1,
        help="Path to LOONE nutrient file to be created.",
    )
    argparser.add_argument(
        "loads_external_filename",
        nargs=1,
        help="File name that holds the external loads for the model.",
    )
    argparser.add_argument(
        "flow_df_filename",
        nargs=1,
        help="The name of the file that holds the flow data for the model.",
    )
    argparser.add_argument(
        "--loone_q_path",
        nargs=1,
        help="Path to LOONE Q CSV file.",
    )
    argparser.add_argument(
        "--data_input",
        nargs=1,
        help="Path to data directory.",
    )
    argparser.add_argument(
        "--forecast_mode",
        action="store_true",
        help="Flag to indicate that the model is running in forecast mode.",
    )
    args = argparser.parse_args()
    workspace = args.workspace[0]
    out_file_name = args.out_file_name[0]
    loads_external_filename = args.loads_external_filename[0]
    flow_df_filename = args.flow_df_filename[0]
    loone_q_path = args.loone_q_path[0] if args.loone_q_path else None
    data_input = args.data_input[0] if args.data_input else None
    forecast_mode = args.forecast_mode if args.forecast_mode else False

    LOONE_NUT(
        workspace,
        out_file_name,
        loads_external_filename,
        flow_df_filename,
        loone_q_path,
        data_input,
        forecast_mode,
    )
