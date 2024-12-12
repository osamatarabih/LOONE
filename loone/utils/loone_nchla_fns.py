# Nitrogen Model Functions

# For NO, NH4
def f_T_alt1(T, T_opt, T_min, T_max):
    import math
    if T < T_opt:
        f_T = math.exp(-2.3 * ((T - T_opt) / (T_opt - T_min)) ** 2)
    elif T >= T_opt:
        f_T = math.exp(-2.3 * ((T - T_opt) / (T_max - T_opt)) ** 2)
    return f_T


# For Chla
def f_T__Chla_alt1(T, T_opt, T_min, T_max, month):
    import math
    if T < T_opt and month in (6, 7, 8, 9, 10):
        f_T = math.exp(-2.3 * ((T - T_opt) / (T_opt - T_min)) ** 2) * 1.2
    elif T < T_opt and month in (11, 12, 1, 2, 3, 4, 5):
        f_T = math.exp(-2.3 * ((T - T_opt) / (T_opt - T_min)) ** 2) * 0.8
    elif T >= T_opt and month in (6, 7, 8, 9, 10):
        f_T = math.exp(-2.3 * ((T - T_opt) / (T_max - T_opt)) ** 2) * 1.2
    elif T >= T_opt and month in (11, 12, 1, 2, 3, 4, 5):
        f_T = math.exp(-2.3 * ((T - T_opt) / (T_max - T_opt)) ** 2) * 0.8
    return f_T


def f_T_alt2(T, T_opt, KTg1, KTg2):
    import math
    if T <= T_opt:
        f_T = math.exp(-KTg1 * ((T - T_opt) ** 2))
    else:
        f_T = math.exp(-KTg2 * ((T_opt - T) ** 2))
    return f_T


def f_T_alt3(T, T_ref, T_min):
    if T <= T_min:
        f_T = 0
    else:
        f_T = (T - T_min) / (T_ref - T_min)
    return f_T


def f_T_alt4(T, T_opt, k1, k2):
    import math
    if T <= T_opt:
        f_T = math.exp(-k1 * (T - T_opt) ** 2)
    else:
        f_T = math.exp(-k2 * (T - T_opt) ** 2)
    return f_T


def f_L_alt1(TD, RAD, Kw, Kc, Chla, z, K1, K2):
    import math
    I_z = (27.25 / TD) * RAD * math.exp(-(Kw + Kc * Chla) * z)
    f_L = (I_z * (1 + 2 * math.sqrt(K1 / K2))) / (I_z + K1 + I_z ** 2 / K2)
    return f_L


def f_L_alt2(fd, Kw, Kc, Chla, z, Io, Im):
    import math
    Ke = Kw + Kc * Chla
    f_L = (2.72 * fd / (Ke * z)) * (math.exp(-Io * math.exp(-Ke * z) / Im) - math.exp(-Io / Im))
    return f_L


def f_L_alt4(fd, Kw, z, PAR, Io, K_LP):
    import math
    PAR_d = PAR * Io / fd
    f_L = (fd / (Kw * z)) * math.log1p((K_LP + PAR_d) / (K_LP + PAR_d * math.exp(-Kw * z)))
    return f_L


def Nit_Rate(Method, Tot_nit_R, kni0, kni, NH4, K_NH, DO, KDO, fox_min, DO_Cr_n, DO_Opt_n, a):
    if Method == 'Simple':
        Nit_R = Tot_nit_R * NH4
    elif Method == 'Opt1':
        Nit_R = kni0 + kni * (NH4 / (K_NH + NH4)) * (DO / (KDO + DO))
    elif Method == 'Opt2':
        Nit_R = kni0 + kni * NH4 * ((1 - fox_min) * ((DO - DO_Cr_n) / (DO_Opt_n - DO_Cr_n)) ** (10 ** a) + fox_min)
    return Nit_R


def Denit_Rate(Method, Tot_denit_R, kden0, kden, NOx, KN, DO, KDO, DO_Cr_d, DO_Opt_d):
    if Method == 'Simple':
        Denit_R = Tot_denit_R * NOx
    elif Method == 'Opt1':
        Denit_R = kden0 + kden * (NOx / (KN + NOx)) * (1 - (DO / (KDO + DO)))
    elif Method == 'Opt2':
        Denit_R = kden0 + kden * NOx * ((DO_Cr_d - DO) / (DO_Cr_d - DO_Opt_d))
    return Denit_R


def NOx_Uptake_Fn(G_max, f_T, f_L, f_P, f_N, K_NH, K_TN, YNOChla, NH4, NOx, Chla):
    Uptake = (G_max * f_T * min(f_L, f_P, f_N) * (1 - (NH4 / (K_NH + NH4))) * ((NH4 + NOx) / (K_TN + NH4 + NOx)) * YNOChla * Chla)
    return Uptake


def NO_Sed(S_NO, Area, NO_Temp_Adj, Vol):
    NO_Sed = (S_NO * Area * NO_Temp_Adj) / Vol
    return NO_Sed


def NOx_Dynamics_alt1(L_ext, Q_o, NOx, kn_an, k_deni_T, NH4, Vol, G_max, f_T, f_L, f_P, f_N, K_NH, K_TN, YNOChla, Chla, S_NO, Area, NO_Temp_Adj):
    NOx_Nxt = ((L_ext - Q_o * NOx + (714.6 * 1000 * 1000) + kn_an * Vol - k_deni_T * Vol - (G_max * f_T * min(f_L, f_P, f_N) * (1 - (NH4 / (K_NH + NH4))) * ((NH4 + NOx) / (K_TN + NH4 + NOx)) * YNOChla * Chla) * Vol + S_NO * Area * NO_Temp_Adj) / Vol) + NOx
    return NOx_Nxt if NOx_Nxt >= 0 else 0


def NOx_Dynamics_N(L_ext, Atm_dep_N, Q_NS, NO_N, kn_an, k_deni_T, NH4_N, Vol_N, G_max, f_T, f_L, f_P_N, f_N_N, K_NH, K_TN, YNOChla, Chla_N, S_NO, Area_N, NO_Temp_Adj):
    NOx_Nxt = ((L_ext - Q_NS * NO_N + Atm_dep_N + kn_an * Vol_N - k_deni_T * Vol_N - (G_max * f_T * min(f_L, f_P_N, f_N_N) * (1 - (NH4_N / (K_NH + NH4_N))) * ((NH4_N + NO_N) / (K_TN + NH4_N + NO_N)) * YNOChla * Chla_N) * Vol_N + S_NO * Area_N * NO_Temp_Adj) / Vol_N) + NO_N
    return NOx_Nxt if NOx_Nxt >= 0 else 0


def NOx_Dynamics_S(Atm_dep_S, Q_N2S, Q_o, NO_N, NO_S, kn_an, k_deni_T, NH4_S, Vol_S, G_max, f_T, f_L, f_P_S, f_N_S, K_NH, K_TN, YNOChla, Chla_S, S_NO, Area_S, NO_Temp_Adj):
    NOx_Nxt = ((Q_N2S * NO_N - Q_o * NO_S + Atm_dep_S + kn_an * Vol_S - k_deni_T * Vol_S - (G_max * f_T * min(f_L, f_P_S, f_N_S) * (1 - (NH4_S / (K_NH + NH4_S))) * ((NH4_S + NO_S) / (K_TN + NH4_S + NO_S)) * YNOChla * Chla_S) * Vol_S + S_NO * Area_S * NO_Temp_Adj) / Vol_S) + NO_S
    return NOx_Nxt if NOx_Nxt >= 0 else 0


def NOx_N_DiffEq(NO_N, t, L_ext, Atm_dep_N, Q_NS, kn_an, k_deni_T, NH4_N, Vol_N, G_max, f_T, f_L, f_P_N, f_N_N, K_NH, K_TN, YNOChla, Chla_N, S_NO, Area_N, NO_Temp_Adj):
    dNOxdt = ((L_ext - Q_NS * NO_N + Atm_dep_N + kn_an * Vol_N - k_deni_T * Vol_N - (G_max * f_T * min(f_L, f_P_N, f_N_N) * (1 - (NH4_N / (K_NH + NH4_N))) * ((NH4_N + NO_N) / (K_TN + NH4_N + NO_N)) * YNOChla * Chla_N) * Vol_N + S_NO * Area_N * NO_Temp_Adj) / Vol_N)
    return dNOxdt


def NH4_N_DiffEq(NH4_N, t, L_ext, Atm_dep_N, Q_NS, kn_an, NO_N, Vol, G_max, f_T, f_L, f_P, f_N, K_NH, K_TN, Chla_N, YNHChla, K_r_T, S_NH4, Area_N, NO_Temp_Adj):
    dNH4dt = ((L_ext - Q_NS * NH4_N + Atm_dep_N - kn_an * NH4_N * Vol - (G_max * f_T * min(f_L, f_P, f_N) * (1 - (NH4_N / (K_NH + NH4_N))) * ((NH4_N + NO_N) / (K_TN + NH4_N + NO_N)) * YNHChla * Chla_N) * Vol + YNHChla * K_r_T * Chla_N * Vol + S_NH4 * Area_N * NO_Temp_Adj) / Vol)
    return dNH4dt


def Chla_N_alt1(Load, Q_NS, Chla_N, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P_N, f_N_N, Vol_N):
    Chla_Nxt = ((Load - Q_NS * Chla_N - vc * Chla_N * Vol_N / z - K_m_T * Chla_N * Vol_N - K_r_T * Chla_N * Vol_N - Grazing * Vol_N + G_max * f_T * min(f_L, f_P_N, f_N_N) * Chla_N * Vol_N) / Vol_N) + Chla_N
    return Chla_Nxt if Chla_Nxt >= 5 else 5


def Chla_N_alt2(Load, Q_NS, Chla_N, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P_N, f_N_N, Vol_N):
    Chla_Nxt = ((Load - Q_NS * Chla_N - vc * Chla_N * Vol_N / z - K_m_T * Chla_N * Vol_N - K_r_T * Chla_N * Vol_N - Grazing * Vol_N + G_max * f_T * f_L * min(f_P_N, f_N_N) * Chla_N * Vol_N) / Vol_N) + Chla_N
    return Chla_Nxt if Chla_Nxt >= 5 else 5


def Chla_Growth_N(Chla_N, G_max, f_T, f_L, f_P_N, f_N_N, Vol_N):
    Growth = (G_max * f_T * f_L * min(f_P_N, f_N_N) * Chla_N * Vol_N) / Vol_N
    return Growth if Growth >= 5 else 5


def NOx_S_DiffEq(NO_S, t, Atm_dep_S, Q_N2S, Q_o, NO_N, kn_an, k_deni_T, NH4_S, Vol_S, G_max, f_T, f_L, f_P_S, f_N_S, K_NH, K_TN, YNOChla, Chla_S, S_NO, Area_S, NO_Temp_Adj):
    dNOxdt = ((Q_N2S * NO_N - Q_o * NO_S + Atm_dep_S + kn_an * Vol_S - k_deni_T * Vol_S - (G_max * f_T * min(f_L, f_P_S, f_N_S) * (1 - (NH4_S / (K_NH + NH4_S))) * ((NH4_S + NO_S) / (K_TN + NH4_S + NO_S)) * YNOChla * Chla_S) * Vol_S + S_NO * Area_S * NO_Temp_Adj) / Vol_S)
    return dNOxdt


def NH4_S_DiffEq(NH4_S, t, Atm_dep_S, Q_N2S, Q_o, NH4_N, kn_an, NO_S, Vol, G_max, f_T, f_L, f_P, f_N, K_NH, K_TN, Chla_S, YNHChla, K_r_T, S_NH4, Area_S, NO_Temp_Adj):
    dNH4dt = ((Q_N2S * NH4_N - Q_o * NH4_S + Atm_dep_S - kn_an * NH4_S * Vol - (G_max * f_T * min(f_L, f_P, f_N) * (1 - (NH4_S / (K_NH + NH4_S))) * ((NH4_S + NO_S) / (K_TN + NH4_S + NO_S)) * YNHChla * Chla_S) * Vol + YNHChla * K_r_T * Chla_S * Vol + + S_NH4 * Area_S * NO_Temp_Adj) / Vol)
    return dNH4dt


def Chla_S_alt1(Q_N2S, Q_o, Chla_N, Chla_S, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P_S, f_N_S, Vol_S):
    Chla_Nxt = ((Q_N2S * Chla_N - Q_o * Chla_S - vc * Chla_S * Vol_S / z - K_m_T * Chla_S * Vol_S - K_r_T * Chla_S * Vol_S - Grazing * Vol_S + G_max * f_T * min(f_L, f_P_S, f_N_S) * Chla_S * Vol_S) / Vol_S) + Chla_S
    return Chla_Nxt if Chla_Nxt >= 5 else 5


def Chla_S_alt2(Q_N2S, Q_o, Chla_N, Chla_S, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P_S, f_N_S, Vol_S):
    Chla_Nxt = ((Q_N2S * Chla_N - Q_o * Chla_S - vc * Chla_S * Vol_S / z - K_m_T * Chla_S * Vol_S - K_r_T * Chla_S * Vol_S - Grazing * Vol_S + G_max * f_T * f_L * min(f_P_S, f_N_S) * Chla_S * Vol_S) / Vol_S) + Chla_S
    return Chla_Nxt if Chla_Nxt >= 5 else 5


def Chla_Growth_S(Chla_S, G_max, f_T, f_L, f_P_S, f_N_S, Vol_S):
    Growth = (G_max * f_T * f_L * min(f_P_S, f_N_S) * Chla_S * Vol_S) / Vol_S
    return Growth if Growth >= 0 else 0


# def NOx_N_DiffEq(NO_N, t, L_ext, Atm_dep_N, Q_NS, kn_an, k_deni_T, NH4_N, Vol_N, G_max, f_T, f_L, f_P_N, f_N_N, K_NH, K_TN, YNOChla, Chla_N, S_NO, Area_N, NO_Temp_Adj):
#     dNOxdt = ((L_ext - Q_NS * NO_N + Atm_dep_N + kn_an * Vol_N - k_deni_T * Vol_N - (G_max * f_T * f_L * min(f_P_N, f_N_N) * (1 - (NH4_N / (K_NH + NH4_N))) * ((NH4_N + NO_N) / (K_TN + NH4_N + NO_N)) * YNOChla * Chla_N) * Vol_N + S_NO * Area_N * NO_Temp_Adj) / Vol_N)
#     return dNOxdt

# def NOx_S_DiffEq(NO_S, t, Atm_dep_S, Q_N2S, Q_o, NO_N, kn_an, k_deni_T, NH4_S, Vol_S, G_max, f_T, f_L, f_P_S, f_N_S, K_NH, K_TN, YNOChla, Chla_S, S_NO, Area_S, NO_Temp_Adj):
#     dNOxdt = ((Q_N2S * NO_N - Q_o * NO_S + Atm_dep_S + kn_an * Vol_S - k_deni_T * Vol_S - (G_max * f_T * f_L * min(f_P_S, f_N_S) * (1 - (NH4_S / (K_NH + NH4_S))) * ((NH4_S + NO_S) / (K_TN + NH4_S + NO_S)) * YNOChla * Chla_S) * Vol_S + S_NO * Area_S * NO_Temp_Adj) / Vol_S)
#     return dNOxdt


def NOx_Dynamics_alt2(L_ext, Q_o, NOx, kn_an, k_deni_T, NH4, Vol, G_max, f_T, f_L, f_P, f_N, K_NH, K_TN, YNOChla, Chla, S_NO, Area, NO_Temp_Adj):
    NOx_Nxt = ((L_ext - Q_o * NOx + (714.6 * 1000 * 1000) + kn_an * Vol - k_deni_T * Vol - (G_max * f_T * f_L * min(f_P, f_N) * (1 - (NH4 / (K_NH + NH4))) * ((NH4 + NOx) / (K_TN + NH4 + NOx)) * YNOChla * Chla) * Vol + S_NO * Area * NO_Temp_Adj) / Vol) + NOx
    return NOx_Nxt if NOx_Nxt >= 0 else 0


def NH4_Dynamics_v1(L_ext, Q_o, NOx, kn_an, NH4, Vol, G_max, f_T, f_L, f_P, f_N, K_NH, K_TN, YNOChla, Chla, YNHChla, K_r_T):
    NH4_Nxt = ((L_ext - Q_o * NH4 + (661.7 * 1000 * 1000) - kn_an * NH4 * Vol - (G_max * f_T * min(f_L, f_P, f_N) * (1 - (NH4 / (K_NH + NH4))) * ((NH4 + NOx) / (K_TN + NH4 + NOx)) * YNOChla * Chla) * Vol + YNHChla * K_r_T * Chla * Vol) / Vol) + NH4
    return NH4_Nxt if NH4_Nxt >= 0 else 0


def Chl_a_Dynamics_alt1(Load, Q_o, Chla, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P, f_N, Vol):
    Chla_Nxt = ((Load - Q_o * Chla - vc * Chla * Vol / z - K_m_T * Chla * Vol - K_r_T * Chla * Vol - Grazing * Vol + G_max * f_T * min(f_L, f_P, f_N) * Chla * Vol) / Vol) + Chla
    return Chla_Nxt if Chla_Nxt >= 0 else 0


def Chl_a_Dynamics_alt2(Load, Q_o, Chla, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P, f_N, Vol):
    Chla_Nxt = ((Load - Q_o * Chla - vc * Chla * Vol / z - K_m_T * Chla * Vol - K_r_T * Chla * Vol - Grazing * Vol + G_max * f_T * f_L * min(f_P, f_N) * Chla * Vol) / Vol) + Chla
    return Chla_Nxt if Chla_Nxt >= 0 else 0


def Chl_a_Dynamics_alt3(Load, Q_o, Chla, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P, f_N, Vol):
    Chla_Nxt = ((Load - Q_o * Chla - vc * Chla * Vol / z - K_m_T * Chla * Vol - K_r_T * Chla * Vol - Grazing * Vol + G_max * f_T * f_L * min(f_P, f_N) * Chla * Vol) / Vol) + Chla
    return Chla_Nxt if Chla_Nxt >= 0 else 0


def Chl_a_Dynamics_alt4(Load, Q_o, Chla, vc, z, K_m_T, K_r_T, Grazing, G_max, f_T, f_L, f_P, f_N, Vol):
    Chla_Nxt = ((Load - Q_o * Chla - vc * Chla * Vol / z - K_m_T * Chla * Vol - K_r_T * Chla * Vol - Grazing * Vol + G_max * f_T * f_L * f_P * f_N * Chla * Vol) / Vol) + Chla
    return Chla_Nxt if Chla_Nxt >= 0 else 0
