from loone.data.tp_variables_regions import TP_Variables as TPVarClass


def DIP_Lake(TP_Lake: float) -> float:
    """
    Calculate the Dissolved Inorganic Phosphorus (DIP) concentration in the Lake Water column.

    Args:
        TP_Lake (float): Total phosphorus concentration in the lake (mg/m3).

    Returns:
        float: Calculated DIP concentration in the lake (mg/m3).
    """
    # Linear Regression done for 1973-2019 monthly data
    DIP_L = 8.68935118641328 + 0.244095176900591 * TP_Lake
    return DIP_L


def Des_flux(
    Γ: float,  # Phosphorus concentration in sediment (mg/m3)
    Mass_sed: float,  # Sediment mass in a given area (kg/m2)
    K_des: float,  # Desorption rate constant (yr-1)
) -> float:
    """
    Calculate the total desorptive flux from the sediment (mg/yr)

    Args:
        Γ (float): Phosphorus concentration in sediment (mg/m3)
        Mass_sed (float): Sediment mass in a given area (kg/m2)
        K_des (float): Desorption rate constant (yr-1)

    Returns:
        float: Total desorptive flux (mg/yr)
    """
    # Total desorptive flux (mg/yr)
    # for value of i
    J_des = K_des * Γ * Mass_sed
    return J_des


def Ads_flux(
    DIP_pore: float,  # DIP concentration in the pore water (mg/m3)
    Γ: float,  # Phosphorus concentration in sediment (mg/m3)
    Mass_sed: float,  # Sediment mass in a given area (kg/m2)
    K_ads: float,  # Adsorption rate constant (yr-1)
    Γ_inf: float,  # Equilibrium phosphorus concentration in sediment (mg/m3)
) -> float:
    """
    Calculate the total adsorptive flux from the pore water into the sediment (mg/yr)

    Args:
        DIP_pore (float): DIP concentration in the pore water (mg/m3)
        Γ (float): Phosphorus concentration in sediment (mg/m3)
        Mass_sed (float): Sediment mass in a given area (kg/m2)
        K_ads (float): Adsorption rate constant (yr-1)
        Γ_inf (float): Equilibrium phosphorus concentration in sediment (mg/m3)

    Returns:
        float: Total adsorptive flux (mg/yr)
    """
    # Total adsorptive flux (mg/yr)
    # for value of i
    J_ads = K_ads * DIP_pore * (Γ_inf - Γ) * Mass_sed
    return J_ads


# def P_sed(
#     Lake_O_A: float,  # Lake oxygen area (m2)
#     TP_Lake: float,  # Total phosphorus concentration in the lake (mg/m3)
#     Sed_burial_flux: float,  # Sediment burial flux (mg/yr)
#     P_sed: float,  # Phosphorus concentration in sediment (mg/m3)
#     Mass_sed: float,  # Sediment mass in a given area (kg/m2)
#     K_decomp: float,  # Decomposition rate constant (yr-1)
#     v_settle: float,  # Settling velocity (m/yr)
# ) -> float:
#     """
#     Calculate the next phosphorus concentration in sediment (mg/m3) based on various parameters.

#     Args:
#         Lake_O_A (float): Lake oxygen area (m2)
#         TP_Lake (float): Total phosphorus concentration in the lake (mg/m3)
#         Sed_burial_flux (float): Sediment burial flux (mg/yr)
#         P_sed (float): Phosphorus concentration in sediment (mg/m3)
#         Mass_sed (float): Sediment mass in a given area (kg/m2)
#         K_decomp (float): Decomposition rate constant (yr-1)
#         v_settle (float): Settling velocity (m/yr)

#     Returns:
#         float: Phosphorus concentration in sediment (mg/m3) at the next time step
#     """
#     #for value of i - 1
#     P_sed_Nxt = ((v_settle * Lake_O_A * TP_Lake - Sed_burial_flux - K_decomp * P_sed * Mass_sed)/Mass_sed) + P_sed
#     return(P_sed_Nxt)


## Multiply v_settle * (TP-DIP)
def P_sed(
    Lake_O_A: float,
    TP_Lake: float,
    DIP_Lake: float,
    Sed_burial_flux: float,
    P_sed: float,
    Mass_sed: float,
    K_decomp: float,
    v_settle: float,
) -> float:
    """
    Calculate the next phosphorus concentration in sediment (mg/m3) based on various parameters.

    Args:
        Lake_O_A (float): Lake oxygen area (m2)
        TP_Lake (float): Total phosphorus concentration in the lake (mg/m3)
        DIP_Lake (float): Dissolved inorganic phosphorus concentration in the lake (mg/m3)
        Sed_burial_flux (float): Sediment burial flux (mg/yr)
        P_sed (float): Phosphorus concentration in sediment (mg/m3)
        Mass_sed (float): Sediment mass in a given area (kg/m2)
        K_decomp (float): Decomposition rate constant (yr-1)
        v_settle (float): Settling velocity (m/yr)

    Returns:
        float: Phosphorus concentration in sediment (mg/m3) at the next time step
    """
    # for value of i - 1
    P_sed_Nxt = (
        (
            v_settle * Lake_O_A * (TP_Lake - DIP_Lake)
            - Sed_burial_flux
            - K_decomp * P_sed * Mass_sed
        )
        / Mass_sed
    ) + P_sed
    return P_sed_Nxt


def Sed_burial_flux(
    P_sed: float,  # Phosphorus concentration in sediment (mg/m3)
    Bulk_density: float,  # Bulk density of sediment (g/cm3)
    A_sed: float,  # Sediment area (m2)
    v_burial: float,  # Burial velocity (m/yr)
    Per_H2O: float,  # Percentage of water in sediment (0-100)
) -> float:
    """
    Calculate the sediment burial flux (mg/yr) based on various parameters.

    Args:
        P_sed (float): Phosphorus concentration in sediment (mg/m3)
        Bulk_density (float): Bulk density of sediment (g/cm3)
        A_sed (float): Sediment area (m2)
        v_burial (float): Burial velocity (m/yr)
        Per_H2O (float): Percentage of water in sediment (0-100)

    Returns:
        float: Sediment burial flux (mg/yr)
    """
    # for value of i
    J_sedburial = (
        Bulk_density
        * ((100 - Per_H2O) / 100)
        * 1000
        * P_sed
        * A_sed
        * v_burial
    )
    return J_sedburial


def Sor_P_burialflux(
    Γ: float,  # Sorbed phosphorus concentration (mg/kg)
    Bulk_density: float,  # Bulk density of sediment (g/cm3)
    A_sed: float,  # Sediment area (m2)
    v_burial: float,  # Burial velocity (m/yr)
    Per_H2O: float,  # Percentage of water in sediment (0-100)
) -> float:
    """
    Calculate the sorbed phosphorus burial flux (mg/yr) based on various parameters.

    Args:
        Γ (float): Sorbed phosphorus concentration (mg/kg)
        Bulk_density (float): Bulk density of sediment (g/cm3)
        A_sed (float): Sediment area (m2)
        v_burial (float): Burial velocity (m/yr)
        Per_H2O (float): Percentage of water in sediment (0-100)

    Returns:
        float: Sorbed phosphorus burial flux (mg/yr)
    """
    # for value of i
    J_Γburial = (
        Bulk_density
        * ((100 - Per_H2O) / 100)
        * 1000
        * Γ
        * A_sed
        * v_burial
    )
    return J_Γburial


def Sor_P_conc(
    Ads_flux: float,  # Adsorption flux (mg/yr)
    Des_flux: float,  # Desorption flux (mg/yr)
    Sor_P_burialflux: float,  # Sorbed phosphorus burial flux (mg/yr)
    Γ: float,  # Sorbed phosphorus concentration (mg/kg)
    Mass_sed: float,  # Sediment mass (kg)
) -> float:
    """
    Calculate the sorbed phosphorus concentration (mg/kg) based on various parameters.

    Args:
        Ads_flux (float): Adsorption flux (mg/yr)
        Des_flux (float): Desorption flux (mg/yr)
        Sor_P_burialflux (float): Sorbed phosphorus burial flux (mg/yr)
        Γ (float): Sorbed phosphorus concentration (mg/kg)
        Mass_sed (float): Sediment mass (kg)

    Returns:
        float: Sorbed phosphorus concentration (mg/kg) at the next time step
    """
    # for value of i - 1
    Γ_Nxt = ((Ads_flux - Des_flux - Sor_P_burialflux) / Mass_sed) + Γ
    return Γ_Nxt


def J_decomp(
    K_decomp: float,  # Decomposition rate constant (yr-1)
    P_sed: float,  # Phosphorus concentration in sediment (mg/m3)
    Mass_sed: float,  # Sediment mass (kg)
) -> float:
    """
    Calculate the decomposition flux (mg/yr) based on various parameters.

    Args:
        K_decomp (float): Decomposition rate constant (yr-1)
        P_sed (float): Phosphorus concentration in sediment (mg/m3)
        Mass_sed (float): Sediment mass (kg)

    Returns:
        float: Decomposition flux (mg/yr)
    """
    J_decomp = K_decomp * P_sed * Mass_sed
    return J_decomp


def DIP_pore(
    workspace: str,
    Θ: float,
    DIP_pore: float,
    DIP_Lake: float,
    Des_flux: float,
    Ads_flux: float,
    P_sed: float,
    Mass_sed: float,
    v_diff: float,
    A_sed: float,
    K_decomp: float,
    v_burial: float,
) -> float:
    """
    Calculate the next dissolved inorganic phosphorus (DIP) concentration in the pore water.

    Args:
        workspace (str): The path to the working directory.
        Θ (float): Porosity of the sediment.
        DIP_pore (float): Current DIP concentration in the pore water (mg/m3).
        DIP_Lake (float): DIP concentration in the lake water (mg/m3).
        Des_flux (float): Desorption flux of DIP from the sediment to the pore water (mg/yr).
        Ads_flux (float): Adsorption flux of DIP from the pore water to the sediment (mg/yr).
        P_sed (float): Phosphorus concentration in the sediment (mg/m3).
        Mass_sed (float): Sediment mass (kg).
        v_diff (float): Diffusion coefficient (m2/yr).
        A_sed (float): Surface area of the sediment-water interface (m2).
        K_decomp (float): Decomposition rate constant (yr-1).
        v_burial (float): Burial velocity of the sediment (m/yr).

    Returns:
        float: The next DIP concentration in the pore water (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    # for value of i - 1
    DIP_p_Nxt = (
        (
            -v_diff * Θ * A_sed * (DIP_pore - DIP_Lake)
            + Des_flux
            - Ads_flux
            + K_decomp * P_sed * Mass_sed
            - v_burial * Θ * A_sed * DIP_pore
        )
        / (Θ * TP_Variables.sediment_depth * A_sed)
    ) + DIP_pore
    return DIP_p_Nxt


# def TP_Lake_N(
#     L_ext: float,
#     Atm_Dep_N: float,
#     Θ_M: float,
#     Θ_S: float,
#     Θ_R: float,
#     Θ_P: float,
#     DIP_pore_M_N: float,
#     DIP_pore_S_N: float,
#     DIP_pore_R_N: float,
#     DIP_pore_P_N: float,
#     DIP_Lake_N: float,
#     Q_N2S: float,
#     Lake_O_A_N: float,
#     TP_Lake_N: float,
#     Lake_V_N: float,
#     v_diff_M: float,
#     v_diff_S: float,
#     v_diff_R: float,
#     v_diff_P: float,
#     v_settle: float
# ) -> float:
#     """
#     Calculate the next total phosphorus concentration for the northern lake segment.

#     Args:
#         L_ext (float): External phosphorus loading.
#         Atm_Dep_N (float): Atmospheric deposition of phosphorus.
#         Θ_M (float): Mixing factor for mud.
#         Θ_S (float): Mixing factor for sand.
#         Θ_R (float): Mixing factor for rock.
#         Θ_P (float): Mixing factor for peat.
#         DIP_pore_M_N (float): DIP concentration in the mud pore water.
#         DIP_pore_S_N (float): DIP concentration in the sand pore water.
#         DIP_pore_R_N (float): DIP concentration in the rock pore water.
#         DIP_pore_P_N (float): DIP concentration in the peat pore water.
#         DIP_Lake_N (float): Current DIP concentration in the lake.
#         Q_N2S (float): Flow rate from north to south.
#         Lake_O_A_N (float): Oxygen concentration in the lake.
#         TP_Lake_N (float): Current total phosphorus concentration in the lake.
#         Lake_V_N (float): Volume of the northern lake segment.
#         v_diff_M (float): Vertical diffusion rate for mud.
#         v_diff_S (float): Vertical diffusion rate for sand.
#         v_diff_R (float): Vertical diffusion rate for rock.
#         v_diff_P (float): Vertical diffusion rate for peat.
#         v_settle (float): Settling velocity of phosphorus.

#     Returns:
#         float: The next total phosphorus concentration in the northern lake segment.
#     """
#     # for value of i - 1
#     TP_L_N_Nxt = (
#         (
#             L_ext
#             + Atm_Dep_N
#             + v_diff_M * (DIP_pore_M_N - DIP_Lake_N) * TP_Variables.A_Mud_N * Θ_M
#             + v_diff_S * (DIP_pore_S_N - DIP_Lake_N) * TP_Variables.A_Sand_N * Θ_S
#             + v_diff_R * (DIP_pore_R_N - DIP_Lake_N) * TP_Variables.A_Rock_N * Θ_R
#             + v_diff_P * (DIP_pore_P_N - DIP_Lake_N) * TP_Variables.A_Peat_N * Θ_P
#             - (Q_N2S + v_settle * Lake_O_A_N) * TP_Lake_N
#         ) / Lake_V_N
#     ) + TP_Lake_N
#     return TP_L_N_Nxt


# def TP_Lake_S(
#     Atm_Dep_S: float,
#     Q_N2S: float,
#     TP_Lake_N: float,
#     Θ_M: float,
#     Θ_S: float,
#     Θ_R: float,
#     Θ_P: float,
#     DIP_pore_M_S: float,
#     DIP_pore_S_S: float,
#     DIP_pore_R_S: float,
#     DIP_pore_P_S: float,
#     DIP_Lake_S: float,
#     Q_O: float,
#     Lake_O_A_S: float,
#     TP_Lake_S: float,
#     Lake_V_S: float,
#     v_diff_M: float,
#     v_diff_S: float,
#     v_diff_R: float,
#     v_diff_P: float,
#     v_settle: float
# ) -> float:
#     """
#     Calculate the next total phosphorus concentration for the southern lake segment.

#     Args:
#         Atm_Dep_S (float): Atmospheric deposition of phosphorus.
#         Q_N2S (float): Flow rate from north to south.
#         TP_Lake_N (float): Current total phosphorus concentration in the northern lake segment.
#         Θ_M (float): Mixing factor for mud.
#         Θ_S (float): Mixing factor for sand.
#         Θ_R (float): Mixing factor for rock.
#         Θ_P (float): Mixing factor for peat.
#         DIP_pore_M_S (float): DIP concentration in the mud pore water.
#         DIP_pore_S_S (float): DIP concentration in the sand pore water.
#         DIP_pore_R_S (float): DIP concentration in the rock pore water.
#         DIP_pore_P_S (float): DIP concentration in the peat pore water.
#         DIP_Lake_S (float): Current DIP concentration in the lake.
#         Q_O (float): Flow rate from south to ocean.
#         Lake_O_A_S (float): Oxygen concentration in the lake.
#         TP_Lake_S (float): Current total phosphorus concentration in the lake.
#         Lake_V_S (float): Volume of the southern lake segment.
#         v_diff_M (float): Vertical diffusion rate for mud.
#         v_diff_S (float): Vertical diffusion rate for sand.
#         v_diff_R (float): Vertical diffusion rate for rock.
#         v_diff_P (float): Vertical diffusion rate for peat.
#         v_settle (float): Settling velocity of phosphorus.

#     Returns:
#         float: The next total phosphorus concentration in the southern lake segment.
#     """
#     # for value of i - 1
#     TP_L_S_Nxt = (
#         (
#             Atm_Dep_S
#             + Q_N2S * TP_Lake_N
#             + v_diff_M * (DIP_pore_M_S - DIP_Lake_S) * TP_Variables.A_Mud_S * Θ_M
#             + v_diff_S * (DIP_pore_S_S - DIP_Lake_S) * TP_Variables.A_Sand_S * Θ_S
#             + v_diff_R * (DIP_pore_R_S - DIP_Lake_S) * TP_Variables.A_Rock_S * Θ_R
#             + v_diff_P * (DIP_pore_P_S - DIP_Lake_S) * TP_Variables.A_Peat_S * Θ_P
#             - (Q_O + v_settle * Lake_O_A_S) * TP_Lake_S
#         ) / Lake_V_S
#     ) + TP_Lake_S
#     return TP_L_S_Nxt


# Calculate Settling Phosphorus (mg/m3) explicitly for analysis purposes
def Sett_P(
    TP_Lake: float,  # Total phosphorus concentration in the lake (mg/m3)
    DIP_Lake: float,  # Dissolved inorganic phosphorus concentration in the lake (mg/m3)
    Lake_O_A: float,  # Oxygen concentration in the lake (mg/m3)
    Lake_V: float,  # Volume of the lake (m3)
    v_settle: float,  # Settling velocity of phosphorus (m/yr)
) -> float:
    """
    Calculate the settling phosphorus (mg/m3) explicitly for analysis purposes.

    Args:
        TP_Lake (float): Total phosphorus concentration in the lake (mg/m3).
        DIP_Lake (float): Dissolved inorganic phosphorus concentration in the lake (mg/m3).
        Lake_O_A (float): Oxygen concentration in the lake (mg/m3).
        Lake_V (float): Volume of the lake (m3).
        v_settle (float): Settling velocity of phosphorus (m/yr).

    Returns:
        float: The settling phosphorus concentration (mg/m3) in the lake.
    """
    Sett_P = v_settle * Lake_O_A * (TP_Lake - DIP_Lake) / Lake_V
    return Sett_P


def P_N_to_S(
    Q_N2S: float,  # Flow rate from North to South lake (m3/s)
    TP_Lake_N: float,  # Total phosphorus concentration in the North lake (mg/m3)
    Lake_V_N: float,  # Volume of the North lake (m3)
) -> float:
    """
    Calculate phosphorus flux from North to South lake.

    Args:
        Q_N2S (float): Flow rate from North to South lake (m3/s).
        TP_Lake_N (float): Total phosphorus concentration in the North lake (mg/m3).
        Lake_V_N (float): Volume of the North lake (m3).

    Returns:
        float: Phosphorus flux from North to South lake (mg/s).
    """
    P_N_t_S = Q_N2S * TP_Lake_N / Lake_V_N
    return P_N_t_S


def P_Out(
    Q_O: float,  # Flow rate out of the lake (m3/s)
    TP_Lake_S: float,  # Total phosphorus concentration in the South lake (mg/m3)
    Lake_V_S: float,  # Volume of the South lake (m3)
) -> float:
    """
    Calculate phosphorus flux out of the South lake.

    Args:
        Q_O (float): Flow rate out of the lake (m3/s).
        TP_Lake_S (float): Total phosphorus concentration in the South lake (mg/m3).
        Lake_V_S (float): Volume of the South lake (m3).

    Returns:
        float: Phosphorus flux out of the South lake (mg/s).
    """
    P_Out = Q_O * TP_Lake_S / Lake_V_S
    return P_Out


def Diff_P(
    v_diff: float,  # Diffusion coefficient (m2/s)
    DIP_pore: float,  # Dissolved inorganic phosphorus concentration in the pore water (mg/m3)
    DIP_Lake: float,  # Dissolved inorganic phosphorus concentration in the lake water (mg/m3)
    Θ: float,  # Porosity of the sediment
    A_sed: float,  # Surface area of the sediment-water interface (m2)
    Lake_V: float,  # Volume of the lake (m3)
) -> float:
    """
    Calculate phosphorus flux between the pore water and the lake water due to diffusion.

    Args:
        v_diff (float): Diffusion coefficient (m2/s).
        DIP_pore (float): Dissolved inorganic phosphorus concentration in the pore water (mg/m3).
        DIP_Lake (float): Dissolved inorganic phosphorus concentration in the lake water (mg/m3).
        Θ (float): Porosity of the sediment.
        A_sed (float): Surface area of the sediment-water interface (m2).
        Lake_V (float): Volume of the lake (m3).

    Returns:
        float: Phosphorus flux between the pore water and the lake water due to diffusion (mg/s).
    """
    Diff_P = v_diff * (DIP_pore - DIP_Lake) * A_sed * Θ / Lake_V
    return Diff_P


#### Multiply V_settle * (TP-DIP_Lake)
def TP_Lake_N(
    workspace: str,  # Path to the workspace directory
    L_ext: float,  # External loading of phosphorus from the catchment to the lake (mg/s)
    Atm_Dep_N: float,  # Atmospheric deposition of phosphorus to the lake (mg/s)
    Θ_M: float,  # Porosity of the mud sediment
    Θ_S: float,  # Porosity of the sand sediment
    Θ_R: float,  # Porosity of the rock sediment
    Θ_P: float,  # Porosity of the peat sediment
    DIP_pore_M_N: float,  # Dissolved inorganic phosphorus concentration in the pore water of mud sediment (mg/m3)
    DIP_pore_S_N: float,  # Dissolved inorganic phosphorus concentration in the pore water of sand sediment (mg/m3)
    DIP_pore_R_N: float,  # Dissolved inorganic phosphorus concentration in the pore water of rock sediment (mg/m3)
    DIP_pore_P_N: float,  # Dissolved inorganic phosphorus concentration in the pore water of peat sediment (mg/m3)
    DIP_Lake_N: float,  # Dissolved inorganic phosphorus concentration in the lake water (mg/m3)
    Q_N2S: float,  # Flow rate from the North lake to the South lake (m3/s)
    Lake_O_A_N: float,  # Oxygen concentration in the North lake (mg/m3)
    TP_Lake_N: float,  # Total phosphorus concentration in the North lake (mg/m3)
    Lake_V_N: float,  # Volume of the North lake (m3)
    v_diff_M: float,  # Diffusion coefficient of mud sediment (m2/s)
    v_diff_S: float,  # Diffusion coefficient of sand sediment (m2/s)
    v_diff_R: float,  # Diffusion coefficient of rock sediment (m2/s)
    v_diff_P: float,  # Diffusion coefficient of peat sediment (m2/s)
    v_settle: float,  # Settling velocity of phosphorus in the lake water (m/s)
) -> float:
    """
    Calculate the next total phosphorus concentration in the North lake.

    Parameters:
        workspace (str): Path to the workspace directory.
        L_ext (float): External loading of phosphorus from the catchment to the lake (mg/s).
        Atm_Dep_N (float): Atmospheric deposition of phosphorus to the lake (mg/s).
        Θ_M (float): Porosity of the mud sediment.
        Θ_S (float): Porosity of the sand sediment.
        Θ_R (float): Porosity of the rock sediment.
        Θ_P (float): Porosity of the peat sediment.
        DIP_pore_M_N (float): Dissolved inorganic phosphorus concentration in the pore water of mud sediment (mg/m3).
        DIP_pore_S_N (float): Dissolved inorganic phosphorus concentration in the pore water of sand sediment (mg/m3).
        DIP_pore_R_N (float): Dissolved inorganic phosphorus concentration in the pore water of rock sediment (mg/m3).
        DIP_pore_P_N (float): Dissolved inorganic phosphorus concentration in the pore water of peat sediment (mg/m3).
        DIP_Lake_N (float): Dissolved inorganic phosphorus concentration in the lake water (mg/m3).
        Q_N2S (float): Flow rate from the North lake to the South lake (m3/s).
        Lake_O_A_N (float): Oxygen concentration in the North lake (mg/m3).
        TP_Lake_N (float): Total phosphorus concentration in the North lake (mg/m3).
        Lake_V_N (float): Volume of the North lake (m3).
        v_diff_M (float): Diffusion coefficient of mud sediment (m2/s).
        v_diff_S (float): Diffusion coefficient of sand sediment (m2/s).
        v_diff_R (float): Diffusion coefficient of rock sediment (m2/s).
        v_diff_P (float): Diffusion coefficient of peat sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: Total phosphorus concentration in the North lake (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    # for value of i - 1
    TP_L_N_Nxt = (
        (
            L_ext
            + Atm_Dep_N
            + v_diff_M
            * (DIP_pore_M_N - DIP_Lake_N)
            * TP_Variables.area_mud_north
            * Θ_M
            + v_diff_S
            * (DIP_pore_S_N - DIP_Lake_N)
            * TP_Variables.area_sand_north
            * Θ_S
            + v_diff_R
            * (DIP_pore_R_N - DIP_Lake_N)
            * TP_Variables.area_rock_north
            * Θ_R
            + v_diff_P
            * (DIP_pore_P_N - DIP_Lake_N)
            * TP_Variables.area_peat_north
            * Θ_P
            - (
                Q_N2S * TP_Lake_N
                + v_settle * Lake_O_A_N * (TP_Lake_N - DIP_Lake_N)
            )
        )
        / Lake_V_N
    ) + TP_Lake_N
    return TP_L_N_Nxt


def TP_Lake_S(
    workspace: str,
    Atm_Dep_S: float,
    Q_N2S: float,
    TP_Lake_N: float,
    Θ_M: float,
    Θ_S: float,
    Θ_R: float,
    Θ_P: float,
    DIP_pore_M_S: float,
    DIP_pore_S_S: float,
    DIP_pore_R_S: float,
    DIP_pore_P_S: float,
    DIP_Lake_S: float,
    Q_O: float,
    Lake_O_A_S: float,
    TP_Lake_S: float,
    Lake_V_S: float,
    v_diff_M: float,
    v_diff_S: float,
    v_diff_R: float,
    v_diff_P: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the South lake.

    Args:
        workspace (str): The path to the working directory.
        Atm_Dep_S (float): Atmospheric deposition of phosphorus to the South lake (kg/yr).
        Q_N2S (float): Flow rate from the North lake to the South lake (m3/s).
        TP_Lake_N (float): Total phosphorus concentration in the North lake (mg/m3).
        Θ_M (float): Porosity of the mud sediment.
        Θ_S (float): Porosity of the sand sediment.
        Θ_R (float): Porosity of the rock sediment.
        Θ_P (float): Porosity of the peat sediment.
        DIP_pore_M_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the mud sediment in the South lake (mg/m3).
        DIP_pore_S_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the sand sediment in the South lake (mg/m3).
        DIP_pore_R_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the rock sediment in the South lake (mg/m3).
        DIP_pore_P_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the peat sediment in the South lake (mg/m3).
        DIP_Lake_S (float): Dissolved inorganic phosphorus (DIP) concentration in the South lake (mg/m3).
        Q_O (float): Flow rate out of the South lake (m3/s).
        Lake_O_A_S (float): Oxygen concentration in the South lake (mg/m3).
        TP_Lake_S (float): Total phosphorus concentration in the South lake (mg/m3).
        Lake_V_S (float): Volume of the South lake (m3).
        v_diff_M (float): Diffusion coefficient of mud sediment (m2/s).
        v_diff_S (float): Diffusion coefficient of sand sediment (m2/s).
        v_diff_R (float): Diffusion coefficient of rock sediment (m2/s).
        v_diff_P (float): Diffusion coefficient of peat sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns
        float: The next total phosphorus concentration in the South lake (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    # for value of i - 1
    TP_L_S_Nxt = (
        (
            Atm_Dep_S
            + Q_N2S * TP_Lake_N
            + v_diff_M
            * (DIP_pore_M_S - DIP_Lake_S)
            * TP_Variables.area_mud_south
            * Θ_M
            + v_diff_S
            * (DIP_pore_S_S - DIP_Lake_S)
            * TP_Variables.area_sand_south
            * Θ_S
            + v_diff_R
            * (DIP_pore_R_S - DIP_Lake_S)
            * TP_Variables.area_rock_south
            * Θ_R
            + v_diff_P
            * (DIP_pore_P_S - DIP_Lake_S)
            * TP_Variables.area_peat_south
            * Θ_P
            - (
                Q_O * TP_Lake_S
                + v_settle * Lake_O_A_S * (TP_Lake_S - DIP_Lake_S)
            )
        )
        / Lake_V_S
    ) + TP_Lake_S
    return TP_L_S_Nxt


# TP_Lake Function of almost all parameters (i.e. I substituted for some parameters in the main function to their basic parameters (e.g. DIP_Pore, P_Sed, etc.))
# def TP_Lake_4Cal(v_diff,K_decomp,v_settle,K_des,K_ads,v_burial,Z_sed,Θ,A_mud,Mass_sed,L_ext_i_1,Γ_i_2, DIP_Lake_i_2, DIP_Lake_i_1, Lake_O_A_i_3, Lake_O_A_i_1, TP_Lake_i_3, TP_Lake_i_1, J_sedburial_i_3, P_sed_i_3, DIP_pore_i_2, Q_o_i_1, Lake_V_i_1):
#     #TP_Lake_i
#     model = ((L_ext_i_1 + v_diff * ((((-v_diff * Θ * A_mud * (DIP_pore_i_2 - DIP_Lake_i_2) + (K_des*Γ_i_2*Mass_sed) - (K_ads*DIP_pore_i_2*(TP_Variables.Γ_inf - Γ_i_2)*Mass_sed) + K_decomp * (((v_settle * Lake_O_A_i_3 * TP_Lake_i_3 - J_sedburial_i_3 - K_decomp * P_sed_i_3 * Mass_sed)/Mass_sed) + P_sed_i_3) * Mass_sed - v_burial * Θ * A_mud * DIP_pore_i_2)/(Θ * Z_sed * A_mud)) + DIP_pore_i_2) - DIP_Lake_i_1) * A_mud * Θ - (Q_o_i_1 + v_settle * Lake_O_A_i_1)*TP_Lake_i_1)/Lake_V_i_1) + TP_Lake_i_1
#     return (model)
# def TP_Lake_4Cal(v_diff,K_decomp,v_settle,K_des,K_ads,v_burial,Z_sed,Θ,A_mud,Mass_sed,L_ext_i_1,Γ_i_2, DIP_Lake_i_2, DIP_Lake_i_1, Lake_O_A_i_3, Lake_O_A_i_1, TP_Lake_i_3, TP_Lake_i_1, J_sedburial_i_3, P_sed_i_3, DIP_pore_i_2, Q_o_i_1, Lake_V_i_1):
#     DIP_p_i_1 = (((-v_diff * Θ * A_mud * (DIP_pore_i_2 - DIP_Lake_i_2) + (K_des*Γ_i_2*Mass_sed) - (K_ads*DIP_pore_i_2*(TP_Variables.Γ_inf - Γ_i_2)*Mass_sed) + K_decomp * (((v_settle * Lake_O_A_i_3 * TP_Lake_i_3 - J_sedburial_i_3 - K_decomp * P_sed_i_3 * Mass_sed)/Mass_sed) + P_sed_i_3) * Mass_sed - v_burial * Θ * A_mud * DIP_pore_i_2)/(Θ * Z_sed * A_mud)) + DIP_pore_i_2) if (((-v_diff * Θ * A_mud * (DIP_pore_i_2 - DIP_Lake_i_2) + (K_des*Γ_i_2*Mass_sed) - (K_ads*DIP_pore_i_2*(TP_Variables.Γ_inf - Γ_i_2)*Mass_sed) + K_decomp * (((v_settle * Lake_O_A_i_3 * TP_Lake_i_3 - J_sedburial_i_3 - K_decomp * P_sed_i_3 * Mass_sed)/Mass_sed) + P_sed_i_3) * Mass_sed - v_burial * Θ * A_mud * DIP_pore_i_2)/(Θ * Z_sed * A_mud)) + DIP_pore_i_2) >0 else 0
#     model = ((L_ext_i_1 + v_diff * (DIP_p_i_1 - DIP_Lake_i_1) * A_mud * Θ - (Q_o_i_1 + v_settle * Lake_O_A_i_1)*TP_Lake_i_1)/Lake_V_i_1) + TP_Lake_i_1
#     return (model)


# Determine TP in the 8 regions
def TP_L_M_N(
    workspace: str,
    L_ext: float,
    Atm_Dep_N: float,
    Θ_M: float,
    DIP_pore_M_N: float,
    DIP_Lake_M_N: float,
    Q_N2S: float,
    Lake_O_A_M_N: float,
    TP_Lake_M_N: float,
    Lake_V_M_N: float,
    v_diff_M: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the North mud region.

    Args:
        workspace (str): The path to the working directory.
        L_ext (float): External loading of phosphorus to the North mud region (kg/yr).
        Atm_Dep_N (float): Atmospheric deposition of phosphorus to the North mud region (kg/yr).
        Θ_M (float): Porosity of the mud sediment.
        DIP_pore_M_N (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the mud sediment in the North mud region (mg/m3).
        DIP_Lake_M_N (float): Dissolved inorganic phosphorus (DIP) concentration in the North mud region (mg/m3).
        Q_N2S (float): Flow rate from the North mud region to the South mud region (m3/s).
        Lake_O_A_M_N (float): Oxygen concentration in the North mud region (mg/m3).
        TP_Lake_M_N (float): Total phosphorus concentration in the North mud region (mg/m3).
        Lake_V_M_N (float): Volume of the North mud region (m3).
        v_diff_M (float): Diffusion coefficient of mud sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the North mud region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    # for value of i - 1
    TP_M_N_Nxt = (
        (
            L_ext
            + Atm_Dep_N
            + v_diff_M
            * (DIP_pore_M_N - DIP_Lake_M_N)
            * TP_Variables.area_mud_north
            * Θ_M
            - (
                Q_N2S * TP_Lake_M_N
                + v_settle * Lake_O_A_M_N * (TP_Lake_M_N - DIP_Lake_M_N)
            )
        )
        / Lake_V_M_N
    ) + TP_Lake_M_N
    return TP_M_N_Nxt


def TP_L_S_N(
    workspace: str,
    L_ext: float,
    Atm_Dep_N: float,
    Θ_S: float,
    DIP_pore_S_N: float,
    DIP_Lake_S_N: float,
    Q_N2S: float,
    Lake_O_A_S_N: float,
    TP_Lake_S_N: float,
    Lake_V_S_N: float,
    v_diff_S: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the North sand region.

    Args:
        workspace (str): The path to the working directory.
        L_ext (float): External loading of phosphorus to the North sand region (kg/yr).
        Atm_Dep_N (float): Atmospheric deposition of phosphorus to the North sand region (kg/yr).
        Θ_S (float): Porosity of the sand sediment.
        DIP_pore_S_N (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the sand sediment in the North sand region (mg/m3).
        DIP_Lake_S_N (float): Dissolved inorganic phosphorus (DIP) concentration in the North sand region (mg/m3).
        Q_N2S (float): Flow rate from the North sand region to the South sand region (m3/s).
        Lake_O_A_S_N (float): Oxygen concentration in the North sand region (mg/m3).
        TP_Lake_S_N (float): Total phosphorus concentration in the North sand region (mg/m3).
        Lake_V_S_N (float): Volume of the North sand region (m3).
        v_diff_S (float): Diffusion coefficient of sand sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the North sand region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    # for value of i - 1
    TP_S_N_Nxt = (
        (
            L_ext
            + Atm_Dep_N
            + v_diff_S
            * (DIP_pore_S_N - DIP_Lake_S_N)
            * TP_Variables.area_sand_north
            * Θ_S
            - (
                Q_N2S * TP_Lake_S_N
                + v_settle * Lake_O_A_S_N * (TP_Lake_S_N - DIP_Lake_S_N)
            )
        )
        / Lake_V_S_N
    ) + TP_Lake_S_N
    return TP_S_N_Nxt


def TP_L_R_N(
    workspace: str,
    L_ext: float,
    Atm_Dep_N: float,
    Θ_R: float,
    DIP_pore_R_N: float,
    DIP_Lake_R_N: float,
    Q_N2S: float,
    Lake_O_A_R_N: float,
    TP_Lake_R_N: float,
    Lake_V_R_N: float,
    v_diff_R: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the North rock region.

    Args:
        workspace (str): The path to the working directory.
        L_ext (float): External loading of phosphorus to the North rock region (kg/yr).
        Atm_Dep_N (float): Atmospheric deposition of phosphorus to the North rock region (kg/yr).
        Θ_R (float): Porosity of the rock sediment.
        DIP_pore_R_N (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the rock sediment in the North rock region (mg/m3).
        DIP_Lake_R_N (float): Dissolved inorganic phosphorus (DIP) concentration in the North rock region (mg/m3).
        Q_N2S (float): Flow rate from the North rock region to the South rock region (m3/s).
        Lake_O_A_R_N (float): Oxygen concentration in the North rock region (mg/m3).
        TP_Lake_R_N (float): Total phosphorus concentration in the North rock region (mg/m3).
        Lake_V_R_N (float): Volume of the North rock region (m3).
        v_diff_R (float): Diffusion coefficient of rock sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the North rock region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    # for value of i - 1
    TP_R_N_Nxt = (
        (
            L_ext
            + Atm_Dep_N
            + v_diff_R
            * (DIP_pore_R_N - DIP_Lake_R_N)
            * TP_Variables.area_rock_north
            * Θ_R
            - (
                Q_N2S * TP_Lake_R_N
                + v_settle * Lake_O_A_R_N * (TP_Lake_R_N - DIP_Lake_R_N)
            )
        )
        / Lake_V_R_N
    ) + TP_Lake_R_N
    return TP_R_N_Nxt


def TP_L_P_N(
    workspace: str,
    L_ext: float,
    Atm_Dep_N: float,
    Θ_P: float,
    DIP_pore_P_N: float,
    DIP_Lake_P_N: float,
    Q_N2S: float,
    Lake_O_A_P_N: float,
    TP_Lake_P_N: float,
    Lake_V_P_N: float,
    v_diff_P: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the North peat region.

    Args:
        workspace (str): The path to the working directory.
        L_ext (float): External loading of phosphorus to the North peat region (kg/yr).
        Atm_Dep_N (float): Atmospheric deposition of phosphorus to the North peat region (kg/yr).
        Θ_P (float): Porosity of the peat sediment.
        DIP_pore_P_N (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the peat sediment in the North peat region (mg/m3).
        DIP_Lake_P_N (float): Dissolved inorganic phosphorus (DIP) concentration in the North peat region (mg/m3).
        Q_N2S (float): Flow rate from the North peat region to the South peat region (m3/s).
        Lake_O_A_P_N (float): Oxygen concentration in the North peat region (mg/m3).
        TP_Lake_P_N (float): Total phosphorus concentration in the North peat region (mg/m3).
        Lake_V_P_N (float): Volume of the North peat region (m3).
        v_diff_P (float): Diffusion coefficient of peat sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the North peat region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    # for value of i - 1
    TP_P_N_Nxt = (
        (
            L_ext
            + Atm_Dep_N
            + v_diff_P
            * (DIP_pore_P_N - DIP_Lake_P_N)
            * TP_Variables.area_peat_north
            * Θ_P
            - (
                Q_N2S * TP_Lake_P_N
                + v_settle * Lake_O_A_P_N * (TP_Lake_P_N - DIP_Lake_P_N)
            )
        )
        / Lake_V_P_N
    ) + TP_Lake_P_N
    return TP_P_N_Nxt


def TP_L_M_S(
    workspace: str,
    Atm_Dep_S: float,
    Q_N2S: float,
    TP_Lake_N: float,
    Θ_M: float,
    DIP_pore_M_S: float,
    DIP_Lake_M_S: float,
    Q_O: float,
    Lake_O_A_M_S: float,
    TP_Lake_M_S: float,
    Lake_V_M_S: float,
    v_diff_M: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the South mud region.

    Args:
        workspace (str): The path to the working directory.
        Atm_Dep_S (float): Atmospheric deposition of phosphorus to the South mud region (kg/yr).
        Q_N2S (float): Flow rate from the North mud region to the South mud region (m3/s).
        TP_Lake_N (float): Total phosphorus concentration in the North mud region (mg/m3).
        Θ_M (float): Porosity of the mud sediment.
        DIP_pore_M_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the mud sediment in the South mud region (mg/m3).
        DIP_Lake_M_S (float): Dissolved inorganic phosphorus (DIP) concentration in the South mud region (mg/m3).
        Q_O (float): Flow rate out of the South mud region (m3/s).
        Lake_O_A_M_S (float): Oxygen concentration in the South mud region (mg/m3).
        TP_Lake_M_S (float): Total phosphorus concentration in the South mud region (mg/m3).
        Lake_V_M_S (float): Volume of the South mud region (m3).
        v_diff_M (float): Diffusion coefficient of mud sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the South mud region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    TP_M_S_Nxt = (
        (
            Atm_Dep_S
            + Q_N2S * TP_Lake_N
            + v_diff_M
            * (DIP_pore_M_S - DIP_Lake_M_S)
            * TP_Variables.area_mud_south
            * Θ_M
            - (
                Q_O * TP_Lake_M_S
                + v_settle * Lake_O_A_M_S * (TP_Lake_M_S - DIP_Lake_M_S)
            )
        )
        / Lake_V_M_S
    ) + TP_Lake_M_S
    return TP_M_S_Nxt


def TP_L_S_S(
    workspace: str,
    Atm_Dep_S: float,
    Q_N2S: float,
    TP_Lake_N: float,
    Θ_S: float,
    DIP_pore_S_S: float,
    DIP_Lake_S_S: float,
    Q_O: float,
    Lake_O_A_S_S: float,
    TP_Lake_S_S: float,
    Lake_V_S_S: float,
    v_diff_S: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the South sand region.

    Args:
        workspace (str): The path to the working directory.
        Atm_Dep_S (float): Atmospheric deposition of phosphorus to the South sand region (kg/yr).
        Q_N2S (float): Flow rate from the North sand region to the South sand region (m3/s).
        TP_Lake_N (float): Total phosphorus concentration in the North sand region (mg/m3).
        Θ_S (float): Porosity of the sand sediment.
        DIP_pore_S_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the sand sediment in the South sand region (mg/m3).
        DIP_Lake_S_S (float): Dissolved inorganic phosphorus (DIP) concentration in the South sand region (mg/m3).
        Q_O (float): Flow rate out of the South sand region (m3/s).
        Lake_O_A_S_S (float): Oxygen concentration in the South sand region (mg/m3).
        TP_Lake_S_S (float): Total phosphorus concentration in the South sand region (mg/m3).
        Lake_V_S_S (float): Volume of the South sand region (m3).
        v_diff_S (float): Diffusion coefficient of sand sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the South sand region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    TP_S_S_Nxt = (
        (
            Atm_Dep_S
            + Q_N2S * TP_Lake_N
            + v_diff_S
            * (DIP_pore_S_S - DIP_Lake_S_S)
            * TP_Variables.area_sand_south
            * Θ_S
            - (
                Q_O * TP_Lake_S_S
                + v_settle * Lake_O_A_S_S * (TP_Lake_S_S - DIP_Lake_S_S)
            )
        )
        / Lake_V_S_S
    ) + TP_Lake_S_S
    return TP_S_S_Nxt


def TP_L_R_S(
    workspace: str,
    Atm_Dep_S: float,
    Q_N2S: float,
    TP_Lake_N: float,
    Θ_R: float,
    DIP_pore_R_S: float,
    DIP_Lake_R_S: float,
    Q_O: float,
    Lake_O_A_R_S: float,
    TP_Lake_R_S: float,
    Lake_V_R_S: float,
    v_diff_R: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the South rock region.

    Args:
        workspace (str): The path to the working directory.
        Atm_Dep_S (float): Atmospheric deposition of phosphorus to the South rock region (kg/yr).
        Q_N2S (float): Flow rate from the North rock region to the South rock region (m3/s).
        TP_Lake_N (float): Total phosphorus concentration in the North rock region (mg/m3).
        Θ_R (float): Porosity of the rock sediment.
        DIP_pore_R_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the rock sediment in the South rock region (mg/m3).
        DIP_Lake_R_S (float): Dissolved inorganic phosphorus (DIP) concentration in the South rock region (mg/m3).
        Q_O (float): Flow rate out of the South rock region (m3/s).
        Lake_O_A_R_S (float): Oxygen concentration in the South rock region (mg/m3).
        TP_Lake_R_S (float): Total phosphorus concentration in the South rock region (mg/m3).
        Lake_V_R_S (float): Volume of the South rock region (m3).
        v_diff_R (float): Diffusion coefficient of rock sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the South rock region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    TP_R_S_Nxt = (
        (
            Atm_Dep_S
            + Q_N2S * TP_Lake_N
            + v_diff_R
            * (DIP_pore_R_S - DIP_Lake_R_S)
            * TP_Variables.area_rock_south
            * Θ_R
            - (
                Q_O * TP_Lake_R_S
                + v_settle * Lake_O_A_R_S * (TP_Lake_R_S - DIP_Lake_R_S)
            )
        )
        / Lake_V_R_S
    ) + TP_Lake_R_S
    return TP_R_S_Nxt


def TP_L_P_S(
    workspace: str,
    Atm_Dep_S: float,
    Q_N2S: float,
    TP_Lake_N: float,
    Θ_P: float,
    DIP_pore_P_S: float,
    DIP_Lake_P_S: float,
    Q_O: float,
    Lake_O_A_P_S: float,
    TP_Lake_P_S: float,
    Lake_V_P_S: float,
    v_diff_P: float,
    v_settle: float,
) -> float:
    """
    Calculate the next total phosphorus concentration in the South peat region.

    Args:
        workspace (str): The path to the working directory.
        Atm_Dep_S (float): Atmospheric deposition of phosphorus to the South peat region (kg/yr).
        Q_N2S (float): Flow rate from the North peat region to the South peat region (m3/s).
        TP_Lake_N (float): Total phosphorus concentration in the North peat region (mg/m3).
        Θ_P (float): Porosity of the peat sediment.
        DIP_pore_P_S (float): Dissolved inorganic phosphorus (DIP) concentration in the pore water of the peat sediment in the South peat region (mg/m3).
        DIP_Lake_P_S (float): Dissolved inorganic phosphorus (DIP) concentration in the South peat region (mg/m3).
        Q_O (float): Flow rate out of the South peat region (m3/s).
        Lake_O_A_P_S (float): Oxygen concentration in the South peat region (mg/m3).
        TP_Lake_P_S (float): Total phosphorus concentration in the South peat region (mg/m3).
        Lake_V_P_S (float): Volume of the South peat region (m3).
        v_diff_P (float): Diffusion coefficient of peat sediment (m2/s).
        v_settle (float): Settling velocity of phosphorus in the lake water (m/s).

    Returns:
        float: The next total phosphorus concentration in the South peat region (mg/m3).
    """
    TP_Variables = TPVarClass(workspace)
    TP_P_S_Nxt = (
        (
            Atm_Dep_S
            + Q_N2S * TP_Lake_N
            + v_diff_P
            * (DIP_pore_P_S - DIP_Lake_P_S)
            * TP_Variables.area_peat_south
            * Θ_P
            - (
                Q_O * TP_Lake_P_S
                + v_settle * Lake_O_A_P_S * (TP_Lake_P_S - DIP_Lake_P_S)
            )
        )
        / Lake_V_P_S
    ) + TP_Lake_P_S
    return TP_P_S_Nxt
