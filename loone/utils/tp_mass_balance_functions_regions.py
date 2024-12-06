from loone.data.tp_variables_regions import TP_Variables as TPVarClass


def DIP_Lake(TP_Lake):
    # Dissolved Inorganic Phosphorus concentration in the Lake Water column (mg/m3)
    # DIP_L = -0.00376 + 0.31877 * TP_Lake #From Pollman and James (2011)
    DIP_L = (
        8.68935118641328 + 0.244095176900591 * TP_Lake
    )  # Linear Regression done for 1973-2019 monthly data!
    # DIP_L = -7.1e-09 * TP_Lake**4 + 1.02e-05 * TP_Lake**3 - 0.00667 * TP_Lake**2 + 2.825741 * TP_Lake - 27.7648 * TP_Lake**0.5 + 85.65048
    return DIP_L


def Des_flux(Γ, Mass_sed, K_des):
    # Total desorptive flux (mg/yr)
    # for value of i
    J_des = K_des * Γ * Mass_sed
    return J_des


def Ads_flux(DIP_pore, Γ, Mass_sed, K_ads, Γ_inf):
    # Total adsorptive flux (mg/yr)
    # for value of i
    J_ads = K_ads * DIP_pore * (Γ_inf - Γ) * Mass_sed
    return J_ads


# def P_sed(Lake_O_A,TP_Lake,Sed_burial_flux,P_sed,Mass_sed,K_decomp,v_settle):
#     #for value of i - 1
#     P_sed_Nxt = ((v_settle * Lake_O_A * TP_Lake - Sed_burial_flux - K_decomp * P_sed * Mass_sed)/Mass_sed) + P_sed
#     return(P_sed_Nxt)


## Multiply v_settle * (TP-DIP)
def P_sed(
    Lake_O_A,
    TP_Lake,
    DIP_Lake,
    Sed_burial_flux,
    P_sed,
    Mass_sed,
    K_decomp,
    v_settle,
):
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


def Sed_burial_flux(P_sed, Bulk_density, A_sed, v_burial, Per_H2O):
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


def Sor_P_burialflux(Γ, Bulk_density, A_sed, v_burial, Per_H2O):
    # for value of i
    J_Γburial = (
        Bulk_density * ((100 - Per_H2O) / 100) * 1000 * Γ * A_sed * v_burial
    )
    return J_Γburial


def Sor_P_conc(Ads_flux, Des_flux, Sor_P_burialflux, Γ, Mass_sed):
    # for value of i - 1
    Γ_Nxt = ((Ads_flux - Des_flux - Sor_P_burialflux) / Mass_sed) + Γ
    return Γ_Nxt


def J_decomp(K_decomp, P_sed, Mass_sed):
    J_decomp = K_decomp * P_sed * Mass_sed
    return J_decomp


def DIP_pore(
    workspace,
    Θ,
    DIP_pore,
    DIP_Lake,
    Des_flux,
    Ads_flux,
    P_sed,
    Mass_sed,
    v_diff,
    A_sed,
    K_decomp,
    v_burial,
):
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


# def TP_Lake_N(L_ext,Atm_Dep_N,Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_N,DIP_pore_S_N,DIP_pore_R_N,DIP_pore_P_N,DIP_Lake_N,Q_N2S,Lake_O_A_N,TP_Lake_N,Lake_V_N,v_diff_M,v_diff_S,v_diff_R,v_diff_P,v_settle):
#     #for value of i - 1
#     TP_L_N_Nxt = ((L_ext + Atm_Dep_N + v_diff_M * (DIP_pore_M_N - DIP_Lake_N) * TP_Variables.A_Mud_N * Θ_M + v_diff_S * (DIP_pore_S_N - DIP_Lake_N) * TP_Variables.A_Sand_N * Θ_S + v_diff_R * (DIP_pore_R_N - DIP_Lake_N) * TP_Variables.A_Rock_N * Θ_R + v_diff_P * (DIP_pore_P_N - DIP_Lake_N) * TP_Variables.A_Peat_N * Θ_P - (Q_N2S + v_settle * Lake_O_A_N)*TP_Lake_N)/Lake_V_N) + TP_Lake_N
#     return(TP_L_N_Nxt)
# def TP_Lake_S(Atm_Dep_S,Q_N2S,TP_Lake_N,Θ_M,Θ_S,Θ_R,Θ_P,DIP_pore_M_S,DIP_pore_S_S,DIP_pore_R_S,DIP_pore_P_S,DIP_Lake_S,Q_O,Lake_O_A_S,TP_Lake_S,Lake_V_S,v_diff_M,v_diff_S,v_diff_R,v_diff_P,v_settle):
#     #for value of i - 1
#     TP_L_S_Nxt = ((Atm_Dep_S + Q_N2S * TP_Lake_N + v_diff_M * (DIP_pore_M_S - DIP_Lake_S) * TP_Variables.A_Mud_S * Θ_M + v_diff_S * (DIP_pore_S_S - DIP_Lake_S) * TP_Variables.A_Sand_S * Θ_S + v_diff_R * (DIP_pore_R_S - DIP_Lake_S) * TP_Variables.A_Rock_S * Θ_R + v_diff_P * (DIP_pore_P_S - DIP_Lake_S) * TP_Variables.A_Peat_S * Θ_P - (Q_O + v_settle * Lake_O_A_S)*TP_Lake_S)/Lake_V_S) + TP_Lake_S
#     return(TP_L_S_Nxt)


# Calculate Settling Phosphorus (mg/m3) explicitly for analysis purposes
def Sett_P(TP_Lake, DIP_Lake, Lake_O_A, Lake_V, v_settle):
    Sett_P = v_settle * Lake_O_A * (TP_Lake - DIP_Lake) / Lake_V
    return Sett_P


def P_N_to_S(Q_N2S, TP_Lake_N, Lake_V_N):
    P_N_t_S = Q_N2S * TP_Lake_N / Lake_V_N
    return P_N_t_S


def P_Out(Q_O, TP_Lake_S, Lake_V_S):
    P_Out = Q_O * TP_Lake_S / Lake_V_S
    return P_Out


def Diff_P(v_diff, DIP_pore, DIP_Lake, Θ, A_sed, Lake_V):
    Diff_P = v_diff * (DIP_pore - DIP_Lake) * A_sed * Θ / Lake_V
    return Diff_P


#### Multiply V_settle * (TP-DIP_Lake)
def TP_Lake_N(
    workspace,
    L_ext,
    Atm_Dep_N,
    Θ_M,
    Θ_S,
    Θ_R,
    Θ_P,
    DIP_pore_M_N,
    DIP_pore_S_N,
    DIP_pore_R_N,
    DIP_pore_P_N,
    DIP_Lake_N,
    Q_N2S,
    Lake_O_A_N,
    TP_Lake_N,
    Lake_V_N,
    v_diff_M,
    v_diff_S,
    v_diff_R,
    v_diff_P,
    v_settle,
):
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
    workspace,
    Atm_Dep_S,
    Q_N2S,
    TP_Lake_N,
    Θ_M,
    Θ_S,
    Θ_R,
    Θ_P,
    DIP_pore_M_S,
    DIP_pore_S_S,
    DIP_pore_R_S,
    DIP_pore_P_S,
    DIP_Lake_S,
    Q_O,
    Lake_O_A_S,
    TP_Lake_S,
    Lake_V_S,
    v_diff_M,
    v_diff_S,
    v_diff_R,
    v_diff_P,
    v_settle,
):
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
    workspace,
    L_ext,
    Atm_Dep_N,
    Θ_M,
    DIP_pore_M_N,
    DIP_Lake_M_N,
    Q_N2S,
    Lake_O_A_M_N,
    TP_Lake_M_N,
    Lake_V_M_N,
    v_diff_M,
    v_settle,
):
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
    workspace,
    L_ext,
    Atm_Dep_N,
    Θ_S,
    DIP_pore_S_N,
    DIP_Lake_S_N,
    Q_N2S,
    Lake_O_A_S_N,
    TP_Lake_S_N,
    Lake_V_S_N,
    v_diff_S,
    v_settle,
):
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
    workspace,
    L_ext,
    Atm_Dep_N,
    Θ_R,
    DIP_pore_R_N,
    DIP_Lake_R_N,
    Q_N2S,
    Lake_O_A_R_N,
    TP_Lake_R_N,
    Lake_V_R_N,
    v_diff_R,
    v_settle,
):
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
    workspace,
    L_ext,
    Atm_Dep_N,
    Θ_P,
    DIP_pore_P_N,
    DIP_Lake_P_N,
    Q_N2S,
    Lake_O_A_P_N,
    TP_Lake_P_N,
    Lake_V_P_N,
    v_diff_P,
    v_settle,
):
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
    workspace,
    Atm_Dep_S,
    Q_N2S,
    TP_Lake_N,
    Θ_M,
    DIP_pore_M_S,
    DIP_Lake_M_S,
    Q_O,
    Lake_O_A_M_S,
    TP_Lake_M_S,
    Lake_V_M_S,
    v_diff_M,
    v_settle,
):
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
    workspace,
    Atm_Dep_S,
    Q_N2S,
    TP_Lake_N,
    Θ_S,
    DIP_pore_S_S,
    DIP_Lake_S_S,
    Q_O,
    Lake_O_A_S_S,
    TP_Lake_S_S,
    Lake_V_S_S,
    v_diff_S,
    v_settle,
):
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
    workspace,
    Atm_Dep_S,
    Q_N2S,
    TP_Lake_N,
    Θ_R,
    DIP_pore_R_S,
    DIP_Lake_R_S,
    Q_O,
    Lake_O_A_R_S,
    TP_Lake_R_S,
    Lake_V_R_S,
    v_diff_R,
    v_settle,
):
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
    workspace,
    Atm_Dep_S,
    Q_N2S,
    TP_Lake_N,
    Θ_P,
    DIP_pore_P_S,
    DIP_Lake_P_S,
    Q_O,
    Lake_O_A_P_S,
    TP_Lake_P_S,
    Lake_V_P_S,
    v_diff_P,
    v_settle,
):
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
