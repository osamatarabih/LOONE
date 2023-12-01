import os
import pandas as pd
from datetime import datetime, timedelta
from glob import glob
import numpy as np

from Model_Config import Model_Config
from Pre_defined_Variables import Pre_defined_Variables
from Stg_Sto_Ar import Stg_Sto_Ar
from Data import Data
from TP_Variables_Regions import TP_Variables
import TP_Mass_Balance_Functions_Regions as TP_MBFR


def LOONE_Nut(
    loone_q_path: str, loads_external: str, data_dir: str | None = None,
) -> pd.DataFrame:
    print("LOONE Nut Module is Running!")
    data_dir = data_dir if data_dir else Model_Config.Working_Path
    loone_q = pd.read_csv(loone_q_path)
    # Based on the defined Start and End year, month, and day on the
    # Pre_defined_Variables File, Startdate and enddate are defined.
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    startdate = datetime.now() #datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    year, month, day = map(int, Pre_defined_Variables.enddate_entry)
    enddate = datetime.now() + timedelta(days=15) #datetime(year, month, day).date()
    schedule = Pre_defined_Variables.Schedule

    date_rng_0 = pd.date_range(start=startdate, end=enddate, freq="D")
    Load_ext = pd.read_csv(
        os.path.join(
            data_dir,
            loads_external
        )
    )
    Q_in = pd.read_csv(
        os.path.join(
            data_dir,
            f"LO_Inflows_BK.csv"
        )
    )
    Flow_df = pd.read_csv(
        glob(os.path.join(
            data_dir,
            f"geoglows_flow_df*predicted.csv"
        ))[0]
    )
    Q_O = Flow_df["Outflows"].values
    S77_Q = loone_q["S77_Q"].values
    S308_Q = loone_q["S308_Q"].values
    TotRegSo = Flow_df[["S351_Out", "S352_Out", "S354_Out"]].sum(axis=1) * (
        70.0456 / 86400
    )
    Sto_Stage = pd.read_csv(
        os.path.join(
            data_dir,
            f"Average_LO_Storage_3MLag.csv"
        )
    )
    Stage_LO = Sto_Stage["Stage_ft"].values
    Storage = Sto_Stage["Storage_acft"].values
    n_rows = len(Q_in.index)
    Storage_dev = Data.Stroage_dev_df["DS_dev"]
    L_ext = Load_ext["TP_Loads_In_mg"]  # mg
    Atm_Dep_N = TP_Variables.N_Per * Load_ext["Atm_Loading_mg"]
    Atm_Dep_S = TP_Variables.S_Per * Load_ext["Atm_Loading_mg"]

    # Read Shear Stress driven by Wind Speed
    Wind_ShearStr = pd.read_csv(
        os.path.join(
            data_dir,
            f"WindShearStress.csv"
        )
    )
    W_SS = Wind_ShearStr["ShearStress"]  # Dyne/cm2
    nu_ts = pd.read_csv(
        os.path.join(
            data_dir,
            "nu.csv"
        )
    )
    LO_BL = 0.5  # m (Bed Elevation of LO)
    g = 9.8  # m/s2 gravitational acceleration
    Cal_Res = pd.read_csv(
        os.path.join(
            data_dir,
            f"nondominated_Sol_var.csv"
        )
    )
    Par = Cal_Res["Par"]
    d_c = Par[20]  # m (particle diameter 10 microm /1E6 to convert to m) clay
    d_s = Par[21]  # m sand
    nu_d = nu_ts["nu"]

    R = 1.65  # submerged specific gravity (1.65 for quartz in water)
    C_1_c = Par[16]
    C_2_c = Par[17]
    C_1_s = Par[18]
    C_2_s = Par[19]

    # Parameters associated with sediment resuspension
    E_0 = 1e-4
    E_1 = 2
    E_2 = 3
    Crtcl_ShStr = Par[22]  # 0.32 #Dyne/cm2
    Td = Par[23]  # days
    L_ext_M = np.zeros(n_rows, dtype=object)
    Q_N2S = np.zeros(n_rows, dtype=object)
    Stage2ar = np.zeros(n_rows, dtype=object)
    LO_WD = np.zeros(n_rows, dtype=object)
    Lake_O_Storage_N = np.zeros(n_rows, dtype=object)
    Lake_O_Storage_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_M_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_S_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_R_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_P_N = np.zeros(n_rows, dtype=object)
    Lake_O_A_M_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_S_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_R_S = np.zeros(n_rows, dtype=object)
    Lake_O_A_P_S = np.zeros(n_rows, dtype=object)
    DIP_Lake_N = np.zeros(n_rows, dtype=object)
    DIP_Lake_S = np.zeros(n_rows, dtype=object)
    TP_Lake_Mean = np.zeros(n_rows, dtype=object)
    J_des_M_N = np.zeros(n_rows, dtype=object)
    J_des_S_N = np.zeros(n_rows, dtype=object)
    J_des_R_N = np.zeros(n_rows, dtype=object)
    J_des_P_N = np.zeros(n_rows, dtype=object)
    J_des_M_S = np.zeros(n_rows, dtype=object)
    J_des_S_S = np.zeros(n_rows, dtype=object)
    J_des_R_S = np.zeros(n_rows, dtype=object)
    J_des_P_S = np.zeros(n_rows, dtype=object)
    J_ads_M_N = np.zeros(n_rows, dtype=object)
    J_ads_S_N = np.zeros(n_rows, dtype=object)
    J_ads_R_N = np.zeros(n_rows, dtype=object)
    J_ads_P_N = np.zeros(n_rows, dtype=object)
    J_ads_M_S = np.zeros(n_rows, dtype=object)
    J_ads_S_S = np.zeros(n_rows, dtype=object)
    J_ads_R_S = np.zeros(n_rows, dtype=object)
    J_ads_P_S = np.zeros(n_rows, dtype=object)
    P_sed_M_N = np.zeros(n_rows, dtype=object)
    P_sed_S_N = np.zeros(n_rows, dtype=object)
    P_sed_R_N = np.zeros(n_rows, dtype=object)
    P_sed_P_N = np.zeros(n_rows, dtype=object)
    P_sed_M_S = np.zeros(n_rows, dtype=object)
    P_sed_S_S = np.zeros(n_rows, dtype=object)
    P_sed_R_S = np.zeros(n_rows, dtype=object)
    P_sed_P_S = np.zeros(n_rows, dtype=object)
    J_sedburial_M_N = np.zeros(n_rows, dtype=object)
    J_sedburial_S_N = np.zeros(n_rows, dtype=object)
    J_sedburial_R_N = np.zeros(n_rows, dtype=object)
    J_sedburial_P_N = np.zeros(n_rows, dtype=object)
    J_sedburial_M_S = np.zeros(n_rows, dtype=object)
    J_sedburial_S_S = np.zeros(n_rows, dtype=object)
    J_sedburial_R_S = np.zeros(n_rows, dtype=object)
    J_sedburial_P_S = np.zeros(n_rows, dtype=object)
    J_Γburial_M_N = np.zeros(n_rows, dtype=object)
    J_Γburial_S_N = np.zeros(n_rows, dtype=object)
    J_Γburial_R_N = np.zeros(n_rows, dtype=object)
    J_Γburial_P_N = np.zeros(n_rows, dtype=object)
    J_Γburial_M_S = np.zeros(n_rows, dtype=object)
    J_Γburial_S_S = np.zeros(n_rows, dtype=object)
    J_Γburial_R_S = np.zeros(n_rows, dtype=object)
    J_Γburial_P_S = np.zeros(n_rows, dtype=object)
    Γ_M_N = np.zeros(n_rows, dtype=object)
    Γ_S_N = np.zeros(n_rows, dtype=object)
    Γ_R_N = np.zeros(n_rows, dtype=object)
    Γ_P_N = np.zeros(n_rows, dtype=object)
    Γ_M_S = np.zeros(n_rows, dtype=object)
    Γ_S_S = np.zeros(n_rows, dtype=object)
    Γ_R_S = np.zeros(n_rows, dtype=object)
    Γ_P_S = np.zeros(n_rows, dtype=object)
    DIP_pore_M_N = np.zeros(n_rows, dtype=object)
    DIP_pore_S_N = np.zeros(n_rows, dtype=object)
    DIP_pore_R_N = np.zeros(n_rows, dtype=object)
    DIP_pore_P_N = np.zeros(n_rows, dtype=object)
    DIP_pore_M_S = np.zeros(n_rows, dtype=object)
    DIP_pore_S_S = np.zeros(n_rows, dtype=object)
    DIP_pore_R_S = np.zeros(n_rows, dtype=object)
    DIP_pore_P_S = np.zeros(n_rows, dtype=object)
    TP_Lake_N = np.zeros(n_rows, dtype=object)
    TP_Lake_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_M_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_S_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_R_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_P_N = np.zeros(n_rows, dtype=object)
    Sed_Resusp_M_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_S_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_R_S = np.zeros(n_rows, dtype=object)
    Sed_Resusp_P_S = np.zeros(n_rows, dtype=object)

    J_decomp_M_N = np.zeros(n_rows, dtype=object)
    J_decomp_S_N = np.zeros(n_rows, dtype=object)
    J_decomp_R_N = np.zeros(n_rows, dtype=object)
    J_decomp_P_N = np.zeros(n_rows, dtype=object)
    J_decomp_M_S = np.zeros(n_rows, dtype=object)
    J_decomp_S_S = np.zeros(n_rows, dtype=object)
    J_decomp_R_S = np.zeros(n_rows, dtype=object)
    J_decomp_P_S = np.zeros(n_rows, dtype=object)

    Settling_P_N = np.zeros(n_rows, dtype=object)
    Settling_P_S = np.zeros(n_rows, dtype=object)

    P_diff_M_N = np.zeros(n_rows, dtype=object)
    P_diff_S_N = np.zeros(n_rows, dtype=object)
    P_diff_R_N = np.zeros(n_rows, dtype=object)
    P_diff_P_N = np.zeros(n_rows, dtype=object)
    P_diff_M_S = np.zeros(n_rows, dtype=object)
    P_diff_S_S = np.zeros(n_rows, dtype=object)
    P_diff_R_S = np.zeros(n_rows, dtype=object)
    P_diff_P_S = np.zeros(n_rows, dtype=object)

    Q_I = Q_in["Inflows_cmd"]
    Q_I_M = np.zeros(n_rows, dtype=object)
    Q_O = np.zeros(n_rows, dtype=object)
    Q_O_M = np.zeros(n_rows, dtype=object)

    P_Load_Cal = np.zeros(n_rows, dtype=object)
    P_Load_StL = np.zeros(n_rows, dtype=object)
    P_Load_South = np.zeros(n_rows, dtype=object)

    v_settle_N_c = np.zeros(n_rows, dtype=object)
    v_settle_N_s = np.zeros(n_rows, dtype=object)
    v_settle_N = np.zeros(n_rows, dtype=object)
    v_settle_S_c = np.zeros(n_rows, dtype=object)
    v_settle_S_s = np.zeros(n_rows, dtype=object)
    v_settle_S = np.zeros(n_rows, dtype=object)

    ##Initial Values##
    # S.A. is calculated based on the Lake's previous time step Stage, but for
    # the S.A. at i=0 I used same time step Stage!
    StartStorage = Stg_Sto_Ar.stg2sto(Pre_defined_Variables.startstage, 0)
    Stage2ar[0] = Stg_Sto_Ar.stg2ar(Stage_LO[0], 0)
    Stage2ar[1] = Stg_Sto_Ar.stg2ar(Stage_LO[1], 0)
    Storage[0] = StartStorage  # ac-ft
    Storage[1] = Stg_Sto_Ar.stg2sto(Stage_LO[1], 0)  # ac-ft
    # TP_MassBalanceModel Initial Values.
    TP_Lake_N[0] = 225  # mg/m3
    TP_Lake_S[0] = 275  # mg/m3
    TP_Lake_Mean[0] = (TP_Lake_N[0] + TP_Lake_S[0]) / 2
    Γ_M_N[0] = 25  # mg/kg
    Γ_S_N[0] = 25  # mg/kg
    Γ_R_N[0] = 25  # mg/kg
    Γ_P_N[0] = 25  # mg/kg
    Γ_M_S[0] = 25  # mg/kg
    Γ_S_S[0] = 25  # mg/kg
    Γ_R_S[0] = 25  # mg/kg
    Γ_P_S[0] = 25  # mg/kg
    DIP_pore_M_N[0] = 700  # 760 #mg/m3
    DIP_pore_S_N[0] = 240  # 205 #mg/m3
    DIP_pore_R_N[0] = 240  # 205 #mg/m3
    DIP_pore_P_N[0] = 160  # 160 #mg/m3
    DIP_pore_M_S[0] = 700  # 760 #mg/m3
    DIP_pore_S_S[0] = 240  # 205 #mg/m3
    DIP_pore_R_S[0] = 240  # 205 #mg/m3
    DIP_pore_P_S[0] = 160  # 160 #mg/m3
    P_sed_M_N[0] = 1100  # mg/kg
    P_sed_S_N[0] = 300  # mg/kg
    P_sed_R_N[0] = 300  # mg/kg
    P_sed_P_N[0] = 200  # mg/kg
    P_sed_M_S[0] = 1100  # mg/kg
    P_sed_S_S[0] = 300  # mg/kg
    P_sed_R_S[0] = 300  # mg/kg
    P_sed_P_S[0] = 200  # mg/kg
    Θ_M = 1 - (
        (TP_Variables.Bulk_density_M / TP_Variables.Particle_density_M)
        * ((100 - TP_Variables.Per_H2O_M) / 100)
    )
    Θ_S = 1 - (
        (TP_Variables.Bulk_density_S / TP_Variables.Particle_density_S)
        * ((100 - TP_Variables.Per_H2O_S) / 100)
    )
    Θ_R = 1 - (
        (TP_Variables.Bulk_density_R / TP_Variables.Particle_density_R)
        * ((100 - TP_Variables.Per_H2O_R) / 100)
    )
    Θ_P = 1 - (
        (TP_Variables.Bulk_density_P / TP_Variables.Particle_density_P)
        * ((100 - TP_Variables.Per_H2O_P) / 100)
    )
    # Mass of sediment in surfacial mix Mud layer in the North Region(kg)
    Mass_sed_M_N = (
        TP_Variables.A_Mud_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_M) / 100)
        * TP_Variables.Bulk_density_M
        * 1000
    )
    # Mass of sediment in surfacial mix Sand layer in the North Region(kg)
    Mass_sed_S_N = (
        TP_Variables.A_Sand_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_S) / 100)
        * TP_Variables.Bulk_density_S
        * 1000
    )
    # Mass of sediment in surfacial mix Rock layer in the North Region(kg)
    Mass_sed_R_N = (
        TP_Variables.A_Rock_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_R) / 100)
        * TP_Variables.Bulk_density_R
        * 1000
    )
    # Mass of sediment in surfacial mix Peat layer in the North Region(kg)
    Mass_sed_P_N = (
        TP_Variables.A_Peat_N
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_P) / 100)
        * TP_Variables.Bulk_density_P
        * 1000
    )
    # Mass of sediment in surfacial mix Mud layer in the South Region(kg)
    Mass_sed_M_S = (
        TP_Variables.A_Mud_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_M) / 100)
        * TP_Variables.Bulk_density_M
        * 1000
    )
    # Mass of sediment in surfacial mix Sand layer in the South Region(kg)
    Mass_sed_S_S = (
        TP_Variables.A_Sand_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_S) / 100)
        * TP_Variables.Bulk_density_S
        * 1000
    )
    # Mass of sediment in surfacial mix Rock layer in the South Region(kg)
    Mass_sed_R_S = (
        TP_Variables.A_Rock_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_R) / 100)
        * TP_Variables.Bulk_density_R
        * 1000
    )
    # Mass of sediment in surfacial mix Peat layer in the South Region(kg)
    Mass_sed_P_S = (
        TP_Variables.A_Peat_S
        * TP_Variables.Z_sed
        * ((100 - TP_Variables.Per_H2O_P) / 100)
        * TP_Variables.Bulk_density_P
        * 1000
    )

    for i in range(n_rows - 2):
        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48  # m3/d
            Q_O_M[i] = Q_O[i]
            L_ext_M[i] = L_ext[i] + Q_I_M[i] * TP_Lake_N[i]
        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48  # m3/d
            Q_I_M[i] = Q_I[i]
            L_ext_M[i] = L_ext[i]
        Q_N2S[i] = (Q_I_M[i] + Q_O_M[i]) / 2
        Stage2ar[i + 2] = Stg_Sto_Ar.stg2ar(Stage_LO[i + 2], 0)
        LO_WD[i] = Stage_LO[i] * 0.3048 - LO_BL
        Lake_O_Storage_N[i] = (
            Storage[i] * TP_Variables.N_Per * 4046.85642 * 0.305
        )  # m3
        Lake_O_Storage_S[i] = (
            Storage[i] * TP_Variables.S_Per * 4046.85642 * 0.305
        )  # m3
        Lake_O_A_N[i] = Stage2ar[i] * TP_Variables.N_Per * 4046.85642  # m2
        Lake_O_A_S[i] = Stage2ar[i] * TP_Variables.S_Per * 4046.85642  # m2
        Lake_O_A_M_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Mud_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_S_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Sand_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_R_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Rock_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_P_N[i] = (
            Lake_O_A_N[i]
            * TP_Variables.A_Peat_N
            / (
                TP_Variables.A_Mud_N
                + TP_Variables.A_Sand_N
                + TP_Variables.A_Rock_N
                + TP_Variables.A_Peat_N
            )
        )
        Lake_O_A_M_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Mud_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )
        Lake_O_A_S_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Sand_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )
        Lake_O_A_R_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Rock_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )
        Lake_O_A_P_S[i] = (
            Lake_O_A_S[i]
            * TP_Variables.A_Peat_S
            / (
                TP_Variables.A_Mud_S
                + TP_Variables.A_Sand_S
                + TP_Variables.A_Rock_S
                + TP_Variables.A_Peat_S
            )
        )

        DIP_Lake_N[i] = TP_MBFR.DIP_Lake(TP_Lake_N[i])
        DIP_Lake_S[i] = TP_MBFR.DIP_Lake(TP_Lake_S[i])

        v_settle_N_c[i] = (R * g * d_c**2) / (
            C_1_c * nu_d[i] + (0.75 * C_2_c * R * g * d_c**3) ** 0.5
        )
        v_settle_N_s[i] = (R * g * d_s**2) / (
            C_1_s * nu_d[i] + (0.75 * C_2_s * R * g * d_s**3) ** 0.5
        )
        v_settle_N[i] = v_settle_N_c[i] * (
            (TP_Variables.A_Mud_N + TP_Variables.A_Peat_N) / TP_Variables.A_N
        ) + v_settle_N_s[i] * (
            (TP_Variables.A_Sand_N + TP_Variables.A_Rock_N) / TP_Variables.A_N
        )

        v_settle_S_c[i] = (R * g * d_c**2) / (
            C_1_c * nu_d[i] + (0.75 * C_2_c * R * g * d_c**3) ** 0.5
        )
        v_settle_S_s[i] = (R * g * d_s**2) / (
            C_1_s * nu_d[i] + (0.75 * C_2_s * R * g * d_s**3) ** 0.5
        )
        v_settle_S[i] = v_settle_S_c[i] * (
            (TP_Variables.A_Mud_S + TP_Variables.A_Peat_S) / TP_Variables.A_S
        ) + v_settle_S_s[i] * (
            (TP_Variables.A_Sand_S + TP_Variables.A_Rock_S) / TP_Variables.A_S
        )

        J_des_M_N[i] = TP_MBFR.Des_flux(
            Γ_M_N[i], Mass_sed_M_N, TP_Variables.K_des_M
        )
        J_des_S_N[i] = TP_MBFR.Des_flux(
            Γ_S_N[i], Mass_sed_S_N, TP_Variables.K_des_S
        )
        J_des_R_N[i] = TP_MBFR.Des_flux(
            Γ_R_N[i], Mass_sed_R_N, TP_Variables.K_des_R
        )
        J_des_P_N[i] = TP_MBFR.Des_flux(
            Γ_P_N[i], Mass_sed_P_N, TP_Variables.K_des_P
        )
        J_des_M_S[i] = TP_MBFR.Des_flux(
            Γ_M_S[i], Mass_sed_M_S, TP_Variables.K_des_M
        )
        J_des_S_S[i] = TP_MBFR.Des_flux(
            Γ_S_S[i], Mass_sed_S_S, TP_Variables.K_des_S
        )
        J_des_R_S[i] = TP_MBFR.Des_flux(
            Γ_R_S[i], Mass_sed_R_S, TP_Variables.K_des_R
        )
        J_des_P_S[i] = TP_MBFR.Des_flux(
            Γ_P_S[i], Mass_sed_P_S, TP_Variables.K_des_P
        )

        J_ads_M_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_M_N[i],
            Γ_M_N[i],
            Mass_sed_M_N,
            TP_Variables.K_ads_M,
            TP_Variables.Γ_inf,
        )
        J_ads_S_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_S_N[i],
            Γ_S_N[i],
            Mass_sed_S_N,
            TP_Variables.K_ads_S,
            TP_Variables.Γ_inf,
        )
        J_ads_R_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_R_N[i],
            Γ_R_N[i],
            Mass_sed_R_N,
            TP_Variables.K_ads_R,
            TP_Variables.Γ_inf,
        )
        J_ads_P_N[i] = TP_MBFR.Ads_flux(
            DIP_pore_P_N[i],
            Γ_P_N[i],
            Mass_sed_P_N,
            TP_Variables.K_ads_P,
            TP_Variables.Γ_inf,
        )
        J_ads_M_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_M_S[i],
            Γ_M_S[i],
            Mass_sed_M_S,
            TP_Variables.K_ads_M,
            TP_Variables.Γ_inf,
        )
        J_ads_S_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_S_S[i],
            Γ_S_S[i],
            Mass_sed_S_S,
            TP_Variables.K_ads_S,
            TP_Variables.Γ_inf,
        )
        J_ads_R_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_R_S[i],
            Γ_R_S[i],
            Mass_sed_R_S,
            TP_Variables.K_ads_R,
            TP_Variables.Γ_inf,
        )
        J_ads_P_S[i] = TP_MBFR.Ads_flux(
            DIP_pore_P_S[i],
            Γ_P_S[i],
            Mass_sed_P_S,
            TP_Variables.K_ads_P,
            TP_Variables.Γ_inf,
        )

        J_sedburial_M_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_M_N[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_N,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_sedburial_S_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_S_N[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_N,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_sedburial_R_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_R_N[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_N,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_sedburial_P_N[i] = TP_MBFR.Sed_burial_flux(
            P_sed_P_N[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_N,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )
        J_sedburial_M_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_M_S[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_S,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_sedburial_S_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_S_S[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_S,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_sedburial_R_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_R_S[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_S,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_sedburial_P_S[i] = TP_MBFR.Sed_burial_flux(
            P_sed_P_S[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_S,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )

        Sed_Resusp_M_N[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_M_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_S_N[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_S_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_R_N[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_R_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_P_N[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_P_N[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_M_S[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_M_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_S_S[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_S_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_R_S[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_R_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )
        Sed_Resusp_P_S[i] = (
            (
                (E_0 / Td**E_1)
                * ((W_SS[i] - Crtcl_ShStr) / Crtcl_ShStr) ** E_2
            )
            * 10
            / LO_WD[i]
            * P_sed_P_S[i]
            if W_SS[i] > Crtcl_ShStr
            else 0
        )

        P_sed_M_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_M_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.K_decomp_M,
                v_settle_N[i],
            )
            - Sed_Resusp_M_N[i] * Lake_O_Storage_N[i] / Mass_sed_M_N
            if TP_MBFR.P_sed(
                Lake_O_A_M_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.K_decomp_M,
                v_settle_N[i],
            )
            - Sed_Resusp_M_N[i] * Lake_O_Storage_N[i] / Mass_sed_M_N
            > 0
            else 0
        )
        P_sed_S_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_S_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.K_decomp_S,
                v_settle_N[i],
            )
            - Sed_Resusp_S_N[i] * Lake_O_Storage_N[i] / Mass_sed_S_N
            if TP_MBFR.P_sed(
                Lake_O_A_S_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.K_decomp_S,
                v_settle_N[i],
            )
            - Sed_Resusp_S_N[i] * Lake_O_Storage_N[i] / Mass_sed_S_N
            > 0
            else 0
        )
        P_sed_R_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_R_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.K_decomp_R,
                v_settle_N[i],
            )
            - Sed_Resusp_R_N[i] * Lake_O_Storage_N[i] / Mass_sed_R_N
            if TP_MBFR.P_sed(
                Lake_O_A_R_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.K_decomp_R,
                v_settle_N[i],
            )
            - Sed_Resusp_R_N[i] * Lake_O_Storage_N[i] / Mass_sed_R_N
            > 0
            else 0
        )
        P_sed_P_N[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_P_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.K_decomp_P,
                v_settle_N[i],
            )
            - Sed_Resusp_P_N[i] * Lake_O_Storage_N[i] / Mass_sed_P_N
            if TP_MBFR.P_sed(
                Lake_O_A_P_N[i],
                TP_Lake_N[i],
                DIP_Lake_N[i],
                J_sedburial_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.K_decomp_P,
                v_settle_N[i],
            )
            - Sed_Resusp_P_N[i] * Lake_O_Storage_N[i] / Mass_sed_P_N
            > 0
            else 0
        )
        P_sed_M_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_M_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.K_decomp_M,
                v_settle_S[i],
            )
            - Sed_Resusp_M_S[i] * Lake_O_Storage_S[i] / Mass_sed_M_S
            if TP_MBFR.P_sed(
                Lake_O_A_M_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.K_decomp_M,
                v_settle_S[i],
            )
            - Sed_Resusp_M_S[i] * Lake_O_Storage_S[i] / Mass_sed_M_S
            > 0
            else 0
        )
        P_sed_S_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_S_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.K_decomp_S,
                v_settle_S[i],
            )
            - Sed_Resusp_S_S[i] * Lake_O_Storage_S[i] / Mass_sed_S_S
            if TP_MBFR.P_sed(
                Lake_O_A_S_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.K_decomp_S,
                v_settle_S[i],
            )
            - Sed_Resusp_S_S[i] * Lake_O_Storage_S[i] / Mass_sed_S_S
            > 0
            else 0
        )
        P_sed_R_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_R_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.K_decomp_R,
                v_settle_S[i],
            )
            - Sed_Resusp_R_S[i] * Lake_O_Storage_S[i] / Mass_sed_R_S
            if TP_MBFR.P_sed(
                Lake_O_A_R_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.K_decomp_R,
                v_settle_S[i],
            )
            - Sed_Resusp_R_S[i] * Lake_O_Storage_S[i] / Mass_sed_R_S
            > 0
            else 0
        )
        P_sed_P_S[i + 1] = (
            TP_MBFR.P_sed(
                Lake_O_A_P_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.K_decomp_P,
                v_settle_S[i],
            )
            - Sed_Resusp_P_S[i] * Lake_O_Storage_S[i] / Mass_sed_P_S
            if TP_MBFR.P_sed(
                Lake_O_A_P_S[i],
                TP_Lake_S[i],
                DIP_Lake_S[i],
                J_sedburial_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.K_decomp_P,
                v_settle_S[i],
            )
            - Sed_Resusp_P_S[i] * Lake_O_Storage_S[i] / Mass_sed_P_S
            > 0
            else 0
        )

        J_Γburial_M_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_M_N[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_N,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_Γburial_S_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_S_N[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_N,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_Γburial_R_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_R_N[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_N,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_Γburial_P_N[i] = TP_MBFR.Sor_P_burialflux(
            Γ_P_N[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_N,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )
        J_Γburial_M_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_M_S[i],
            TP_Variables.Bulk_density_M,
            TP_Variables.A_Mud_S,
            TP_Variables.v_burial_M,
            TP_Variables.Per_H2O_M,
        )
        J_Γburial_S_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_S_S[i],
            TP_Variables.Bulk_density_S,
            TP_Variables.A_Sand_S,
            TP_Variables.v_burial_S,
            TP_Variables.Per_H2O_S,
        )
        J_Γburial_R_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_R_S[i],
            TP_Variables.Bulk_density_R,
            TP_Variables.A_Rock_S,
            TP_Variables.v_burial_R,
            TP_Variables.Per_H2O_R,
        )
        J_Γburial_P_S[i] = TP_MBFR.Sor_P_burialflux(
            Γ_P_S[i],
            TP_Variables.Bulk_density_P,
            TP_Variables.A_Peat_S,
            TP_Variables.v_burial_P,
            TP_Variables.Per_H2O_P,
        )

        Γ_M_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_M_N[i],
                J_des_M_N[i],
                J_Γburial_M_N[i],
                Γ_M_N[i],
                Mass_sed_M_N,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_M_N[i],
                J_des_M_N[i],
                J_Γburial_M_N[i],
                Γ_M_N[i],
                Mass_sed_M_N,
            )
            > 0
            else 0
        )
        Γ_S_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_S_N[i],
                J_des_S_N[i],
                J_Γburial_S_N[i],
                Γ_S_N[i],
                Mass_sed_S_N,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_S_N[i],
                J_des_S_N[i],
                J_Γburial_S_N[i],
                Γ_S_N[i],
                Mass_sed_S_N,
            )
            > 0
            else 0
        )
        Γ_R_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_R_N[i],
                J_des_R_N[i],
                J_Γburial_R_N[i],
                Γ_R_N[i],
                Mass_sed_R_N,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_R_N[i],
                J_des_R_N[i],
                J_Γburial_R_N[i],
                Γ_R_N[i],
                Mass_sed_R_N,
            )
            > 0
            else 0
        )
        Γ_P_N[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_P_N[i],
                J_des_P_N[i],
                J_Γburial_P_N[i],
                Γ_P_N[i],
                Mass_sed_P_N,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_P_N[i],
                J_des_P_N[i],
                J_Γburial_P_N[i],
                Γ_P_N[i],
                Mass_sed_P_N,
            )
            > 0
            else 0
        )
        Γ_M_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_M_S[i],
                J_des_M_S[i],
                J_Γburial_M_S[i],
                Γ_M_S[i],
                Mass_sed_M_S,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_M_S[i],
                J_des_M_S[i],
                J_Γburial_M_S[i],
                Γ_M_S[i],
                Mass_sed_M_S,
            )
            > 0
            else 0
        )
        Γ_S_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_S_S[i],
                J_des_S_S[i],
                J_Γburial_S_S[i],
                Γ_S_S[i],
                Mass_sed_S_S,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_S_S[i],
                J_des_S_S[i],
                J_Γburial_S_S[i],
                Γ_S_S[i],
                Mass_sed_S_S,
            )
            > 0
            else 0
        )
        Γ_R_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_R_S[i],
                J_des_R_S[i],
                J_Γburial_R_S[i],
                Γ_R_S[i],
                Mass_sed_R_S,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_R_S[i],
                J_des_R_S[i],
                J_Γburial_R_S[i],
                Γ_R_S[i],
                Mass_sed_R_S,
            )
            > 0
            else 0
        )
        Γ_P_S[i + 1] = (
            TP_MBFR.Sor_P_conc(
                J_ads_P_S[i],
                J_des_P_S[i],
                J_Γburial_P_S[i],
                Γ_P_S[i],
                Mass_sed_P_S,
            )
            if TP_MBFR.Sor_P_conc(
                J_ads_P_S[i],
                J_des_P_S[i],
                J_Γburial_P_S[i],
                Γ_P_S[i],
                Mass_sed_P_S,
            )
            > 0
            else 0
        )

        J_decomp_M_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_M, P_sed_M_N[i], Mass_sed_M_N
        )
        J_decomp_S_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_S, P_sed_S_N[i], Mass_sed_S_N
        )
        J_decomp_R_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_R, P_sed_R_N[i], Mass_sed_R_N
        )
        J_decomp_P_N[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_P, P_sed_P_N[i], Mass_sed_P_N
        )
        J_decomp_M_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_M, P_sed_M_S[i], Mass_sed_M_S
        )
        J_decomp_S_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_S, P_sed_S_S[i], Mass_sed_S_S
        )
        J_decomp_R_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_R, P_sed_R_S[i], Mass_sed_R_S
        )
        J_decomp_P_S[i] = TP_MBFR.J_decomp(
            TP_Variables.K_decomp_P, P_sed_P_S[i], Mass_sed_P_S
        )

        DIP_pore_M_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_N[i],
                DIP_Lake_N[i],
                J_des_M_N[i],
                J_ads_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_N,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            if TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_N[i],
                DIP_Lake_N[i],
                J_des_M_N[i],
                J_ads_M_N[i],
                P_sed_M_N[i],
                Mass_sed_M_N,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_N,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            > 0
            else 0
        )
        DIP_pore_S_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_N[i],
                DIP_Lake_N[i],
                J_des_S_N[i],
                J_ads_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_N,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            if TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_N[i],
                DIP_Lake_N[i],
                J_des_S_N[i],
                J_ads_S_N[i],
                P_sed_S_N[i],
                Mass_sed_S_N,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_N,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            > 0
            else 0
        )
        DIP_pore_R_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_N[i],
                DIP_Lake_N[i],
                J_des_R_N[i],
                J_ads_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_N,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            if TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_N[i],
                DIP_Lake_N[i],
                J_des_R_N[i],
                J_ads_R_N[i],
                P_sed_R_N[i],
                Mass_sed_R_N,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_N,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            > 0
            else 0
        )
        DIP_pore_P_N[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                J_des_P_N[i],
                J_ads_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_N,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            if TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                J_des_P_N[i],
                J_ads_P_N[i],
                P_sed_P_N[i],
                Mass_sed_P_N,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_N,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            > 0
            else 0
        )
        DIP_pore_M_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_S[i],
                DIP_Lake_S[i],
                J_des_M_S[i],
                J_ads_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_S,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            if TP_MBFR.DIP_pore(
                Θ_M,
                DIP_pore_M_S[i],
                DIP_Lake_S[i],
                J_des_M_S[i],
                J_ads_M_S[i],
                P_sed_M_S[i],
                Mass_sed_M_S,
                TP_Variables.v_diff_M,
                TP_Variables.A_Mud_S,
                TP_Variables.K_decomp_M,
                TP_Variables.v_burial_M,
            )
            > 0
            else 0
        )
        DIP_pore_S_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_S[i],
                DIP_Lake_S[i],
                J_des_S_S[i],
                J_ads_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_S,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            if TP_MBFR.DIP_pore(
                Θ_S,
                DIP_pore_S_S[i],
                DIP_Lake_S[i],
                J_des_S_S[i],
                J_ads_S_S[i],
                P_sed_S_S[i],
                Mass_sed_S_S,
                TP_Variables.v_diff_S,
                TP_Variables.A_Sand_S,
                TP_Variables.K_decomp_S,
                TP_Variables.v_burial_S,
            )
            > 0
            else 0
        )
        DIP_pore_R_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_S[i],
                DIP_Lake_S[i],
                J_des_R_S[i],
                J_ads_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_S,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            if TP_MBFR.DIP_pore(
                Θ_R,
                DIP_pore_R_S[i],
                DIP_Lake_S[i],
                J_des_R_S[i],
                J_ads_R_S[i],
                P_sed_R_S[i],
                Mass_sed_R_S,
                TP_Variables.v_diff_R,
                TP_Variables.A_Rock_S,
                TP_Variables.K_decomp_R,
                TP_Variables.v_burial_R,
            )
            > 0
            else 0
        )
        DIP_pore_P_S[i + 1] = (
            TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                J_des_P_S[i],
                J_ads_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_S,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            if TP_MBFR.DIP_pore(
                Θ_P,
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                J_des_P_S[i],
                J_ads_P_S[i],
                P_sed_P_S[i],
                Mass_sed_P_S,
                TP_Variables.v_diff_P,
                TP_Variables.A_Peat_S,
                TP_Variables.K_decomp_P,
                TP_Variables.v_burial_P,
            )
            > 0
            else 0
        )

        Settling_P_N[i] = TP_MBFR.Sett_P(
            TP_Lake_N[i],
            DIP_Lake_N[i],
            Lake_O_A_N[i],
            Lake_O_Storage_N[i],
            v_settle_N[i],
        )
        Settling_P_S[i] = TP_MBFR.Sett_P(
            TP_Lake_S[i],
            DIP_Lake_S[i],
            Lake_O_A_S[i],
            Lake_O_Storage_S[i],
            v_settle_S[i],
        )

        P_diff_M_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_M,
            DIP_pore_M_N[i],
            DIP_Lake_N[i],
            Θ_M,
            TP_Variables.A_Mud_N,
            Lake_O_Storage_N[i],
        )
        P_diff_S_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_S,
            DIP_pore_S_N[i],
            DIP_Lake_N[i],
            Θ_S,
            TP_Variables.A_Sand_N,
            Lake_O_Storage_N[i],
        )
        P_diff_R_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_R,
            DIP_pore_R_N[i],
            DIP_Lake_N[i],
            Θ_R,
            TP_Variables.A_Rock_N,
            Lake_O_Storage_N[i],
        )
        P_diff_P_N[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_P,
            DIP_pore_P_N[i],
            DIP_Lake_N[i],
            Θ_P,
            TP_Variables.A_Peat_N,
            Lake_O_Storage_N[i],
        )
        P_diff_M_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_M,
            DIP_pore_M_S[i],
            DIP_Lake_S[i],
            Θ_M,
            TP_Variables.A_Mud_S,
            Lake_O_Storage_S[i],
        )
        P_diff_S_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_S,
            DIP_pore_S_S[i],
            DIP_Lake_S[i],
            Θ_S,
            TP_Variables.A_Sand_S,
            Lake_O_Storage_S[i],
        )
        P_diff_R_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_R,
            DIP_pore_R_S[i],
            DIP_Lake_S[i],
            Θ_R,
            TP_Variables.A_Rock_S,
            Lake_O_Storage_S[i],
        )
        P_diff_P_S[i] = TP_MBFR.Diff_P(
            TP_Variables.v_diff_P,
            DIP_pore_P_S[i],
            DIP_Lake_S[i],
            Θ_P,
            TP_Variables.A_Peat_S,
            Lake_O_Storage_S[i],
        )

        TP_Lake_N[i + 1] = (
            TP_MBFR.TP_Lake_N(
                L_ext_M[i],
                Atm_Dep_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_N[i],
                DIP_pore_S_N[i],
                DIP_pore_R_N[i],
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                Q_N2S[i],
                Lake_O_A_N[i],
                TP_Lake_N[i],
                Lake_O_Storage_N[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_N[i],
            )
            + (
                Sed_Resusp_M_N[i]
                + Sed_Resusp_S_N[i]
                + Sed_Resusp_R_N[i]
                + Sed_Resusp_P_N[i]
            )
            if TP_MBFR.TP_Lake_N(
                L_ext_M[i],
                Atm_Dep_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_N[i],
                DIP_pore_S_N[i],
                DIP_pore_R_N[i],
                DIP_pore_P_N[i],
                DIP_Lake_N[i],
                Q_N2S[i],
                Lake_O_A_N[i],
                TP_Lake_N[i],
                Lake_O_Storage_N[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_N[i],
            )
            + (
                Sed_Resusp_M_N[i]
                + Sed_Resusp_S_N[i]
                + Sed_Resusp_R_N[i]
                + Sed_Resusp_P_N[i]
            )
            > 0
            else 0
        )
        TP_Lake_S[i + 1] = (
            TP_MBFR.TP_Lake_S(
                Atm_Dep_S[i],
                Q_N2S[i],
                TP_Lake_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_S[i],
                DIP_pore_S_S[i],
                DIP_pore_R_S[i],
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                Q_O_M[i],
                Lake_O_A_S[i],
                TP_Lake_S[i],
                Lake_O_Storage_S[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_S[i],
            )
            + (
                Sed_Resusp_M_S[i]
                + Sed_Resusp_S_S[i]
                + Sed_Resusp_R_S[i]
                + Sed_Resusp_P_S[i]
            )
            if TP_MBFR.TP_Lake_S(
                Atm_Dep_S[i],
                Q_N2S[i],
                TP_Lake_N[i],
                Θ_M,
                Θ_S,
                Θ_R,
                Θ_P,
                DIP_pore_M_S[i],
                DIP_pore_S_S[i],
                DIP_pore_R_S[i],
                DIP_pore_P_S[i],
                DIP_Lake_S[i],
                Q_O_M[i],
                Lake_O_A_S[i],
                TP_Lake_S[i],
                Lake_O_Storage_S[i],
                TP_Variables.v_diff_M,
                TP_Variables.v_diff_S,
                TP_Variables.v_diff_R,
                TP_Variables.v_diff_P,
                v_settle_S[i],
            )
            + (
                Sed_Resusp_M_S[i]
                + Sed_Resusp_S_S[i]
                + Sed_Resusp_R_S[i]
                + Sed_Resusp_P_S[i]
            )
            > 0
            else 0
        )
        TP_Lake_Mean[i + 1] = (TP_Lake_N[i + 1] + TP_Lake_S[i + 1]) / 2

        P_Load_Cal[i] = (
            S77_Q[i] * 0.028316847 * 3600 * 24 * TP_Lake_S[i]
        )  # mg/d P
        P_Load_StL[i] = (
            S308_Q[i] * 0.028316847 * 3600 * 24 * TP_Lake_S[i]
        )  # mg/d P
        #breakpoint()
        P_Load_South[i] = TotRegSo[i] * 1233.48 * TP_Lake_S[i]  # mg/d P

    P_Loads_df = pd.DataFrame(
        date_rng_0, columns=["Date"]
    )  # 1/1/2008-12/31/2018
    P_Lake_df = pd.DataFrame(
        date_rng_0, columns=["Date"]
    )  # 1/1/2008-12/31/2018

    P_Loads_df["P_Load_Cal"] = pd.to_numeric(P_Load_Cal) / 1e9  # tons
    P_Loads_df["P_Load_StL"] = pd.to_numeric(P_Load_StL) / 1e9  # tons
    P_Loads_df["P_Load_South"] = pd.to_numeric(P_Load_South) / 1e9  # tons
    P_Lake_df["P_Lake"] = pd.to_numeric(TP_Lake_Mean)
    P_Lake_df["TP_Lake_S"] = pd.to_numeric(TP_Lake_S)
    P_Loads_df = P_Loads_df.set_index("Date")
    P_Loads_df.index = pd.to_datetime(P_Loads_df.index, unit="ns")
    P_Loads_M = P_Loads_df.resample("M").sum()
    P_Loads_M = P_Loads_M.reset_index()
    P_Lake_df = P_Lake_df.set_index("Date")

    return P_Lake_df
