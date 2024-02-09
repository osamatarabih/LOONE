# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 08:05:45 2022

@author: osama
"""

# This script calculates wind-induced wave parameters and bottom shear stress due to waves!

# Import required PAckages
import pandas as pd
import numpy as np
from scipy.optimize import fsolve
import os
import math

Working_dir = "C:/Work/Research/LOONE/LOONE_Data_Pre"
os.chdir("%s" % Working_dir)
# Read Mean Wind Speed in LO
LO_WS = pd.read_csv("./LOWS.csv")
LO_WS["WS_mps"] = LO_WS["LO_Avg_WS_MPH"] * 0.44704  # MPH to m/s
# Read LO Stage to consider water depth changes
LO_Stage = pd.read_csv("./LO_Stg_Sto_SA_2008-2023.csv")
LO_Stage["Stage_m"] = LO_Stage["Stage_ft"] * 0.3048
Bottom_Elev = 0.5  # m (Karl E. Havens • Alan D. Steinman 2013)
LO_Wd = LO_Stage["Stage_m"] - Bottom_Elev
g = 9.81  # m/s2 gravitational acc
d = 1.5  # m  LO Mean water depth
F = 57500  # Fetch length of wind (m) !!!!!!
nu = 1.0034 / 1e6  # m2/s (kinematic viscosity of water at T = 20 C)
ru = 1000  # kg/m3

n = len(LO_WS.index)


class Wind_Func:
    def H(g, d, F, WS):
        H = (
            (
                0.283
                * np.tanh(0.53 * (g * d / WS**2) ** 0.75)
                * np.tanh(
                    0.00565
                    * (g * F / WS**2) ** 0.5
                    / np.tanh(0.53 * (g * d / WS**2) ** (3 / 8))
                )
            )
            * WS**2
            / g
        )
        return H

    def T(g, d, F, WS):
        T = (
            (
                7.54
                * np.tanh(0.833 * (g * d / WS**2) ** (3 / 8))
                * np.tanh(
                    0.0379
                    * (g * F / WS**2) ** 0.5
                    / np.tanh(0.833 * (g * d / WS**2) ** (3 / 8))
                )
            )
            * WS
            / g
        )
        return T

    def L(g, d, T):

        def func(x):
            return [(g * T**2 / 2 * np.pi) * np.tanh(2 * np.pi * d / x[0]) - x[0]]

        sol = fsolve(func, [1])
        L = sol[0]
        return L


W_H = np.zeros(n, dtype=object)
W_T = np.zeros(n, dtype=object)
W_L = np.zeros(n, dtype=object)
W_ShearStress = np.zeros(n, dtype=object)
for i in range(n):
    W_H[i] = Wind_Func.H(g, LO_Wd[i], F, LO_WS["WS_mps"].iloc[i])
    W_T[i] = Wind_Func.T(g, LO_Wd[i], F, LO_WS["WS_mps"].iloc[i])
    W_L[i] = Wind_Func.L(g, LO_Wd[i], W_T[i])
    W_ShearStress[i] = (
        W_H[i]
        * (ru * (nu * (2 * np.pi / W_T[i]) ** 3) ** 0.5)
        / (2 * np.sinh(2 * np.pi * LO_Wd[i] / W_L[i]))
    )


Wind_ShearStress = pd.DataFrame(LO_WS["date"], columns=["date"])
Wind_ShearStress["ShearStress"] = W_ShearStress * 10  # Convert N/m2 to Dyne/cm2
Wind_ShearStress.to_csv("./WindShearStress.csv")

# #Monthly
# Wind_ShearStress['Date'] = pd.to_datetime(Wind_ShearStress['Date'])
# Wind_ShearStress_df = pd.DataFrame()
# Wind_ShearStress_df['Date'] = Wind_ShearStress['Date'].dt.date
# Wind_ShearStress_df['ShearStress'] = pd.to_numeric(Wind_ShearStress['ShearStress'])
# Wind_ShearStress_df = Wind_ShearStress_df.set_index(['Date'])
# Wind_ShearStress_df.index = pd.to_datetime(Wind_ShearStress_df.index, unit = 'ns')
# Wind_ShearStress_df = Wind_ShearStress_df.resample('M').mean()
# Wind_ShearStress_df.to_csv('C:/Work/Research/Data Analysis/Lake_O_Weather_Data/WindSpeed_Processed/WindShearStress_M.csv')

# The drag coefficient
CD = 0.001 * (0.75 + 0.067 * LO_WS["WS_mps"])
air_ru = 1.293  # kg/m3


def tau_w(WS, CD, air_ru):
    tau_w = air_ru * CD * (WS**2)
    return tau_w


def Current_bottom_shear_stress(ru, tau_w):
    # Constants
    kappa = 0.41  # Von Karman constant
    # Calculate the bottom friction velocity
    u_b = math.sqrt(tau_w / ru)
    # Calculate the bottom shear stress
    tau_b = ru * kappa**2 * u_b**2
    return tau_b


Current_Stress = np.zeros(n, dtype=object)
Wind_Stress = np.zeros(n, dtype=object)
for i in range(n):
    Wind_Stress[i] = tau_w(LO_WS["WS_mps"].iloc[i], CD[i], air_ru)
    Current_Stress[i] = Current_bottom_shear_stress(ru, Wind_Stress[i])

Current_ShearStress_df = pd.DataFrame(LO_WS["date"], columns=["date"])
Current_ShearStress_df["Current_Stress"] = (
    Current_Stress * 10
)  # Convert N/m2 to Dyne/cm2
Current_ShearStress_df["Wind_Stress"] = Wind_Stress * 10  # Convert N/m2 to Dyne/cm2
Current_ShearStress_df["Wind_Speed_m/s"] = LO_WS["WS_mps"]


def Current_bottom_shear_stress_2(u, k, nu, ks, z, ru):
    def func1(u_str1):
        return [u_str1[0] - u * k * np.exp(z / (0.11 * nu / u_str1[0]))]

    sol1 = fsolve(func1, [1])

    def func2(u_str2):
        return [u_str2[0] - u * k * np.exp(z / (0.0333 * ks))]

    sol2 = fsolve(func2, [1])

    def func3(u_str3):
        return [u_str3[0] - u * k * np.exp(z / ((0.11 * nu / u_str3[0]) + 0.0333 * ks))]

    sol3 = fsolve(func3, [1])
    if sol1[0] * ks / nu <= 5:
        u_str = sol1[0]
    elif sol2[0] * ks / nu >= 70:
        u_str = sol2[0]
    elif sol3[0] * ks / nu > 5 and sol3[0] * ks / nu < 70:
        u_str = sol3[0]
    tau_c = ru * u_str**2
    return tau_c


def Current_bottom_shear_stress_3(u, k, nu, ks, z, ru):
    def func1(u_str1):
        return [u_str1[0] - u * k * (1 / np.log(z / (0.11 * nu / u_str1[0])))]

    sol1 = fsolve(func1, [1])

    def func2(u_str2):
        return [u_str2[0] - u * k * (1 / np.log(z / (0.0333 * ks)))]

    sol2 = fsolve(func2, [1])

    def func3(u_str3):
        return [
            u_str3[0]
            - u * k * (1 / np.log(z / ((0.11 * nu / u_str3[0]) + 0.0333 * ks)))
        ]

    sol3 = fsolve(func3, [1])
    if sol1[0] * ks / nu <= 5:
        u_str = sol1[0]
    elif sol2[0] * ks / nu >= 70:
        u_str = sol2[0]
    elif sol3[0] * ks / nu > 5 and sol3[0] * ks / nu < 70:
        u_str = sol3[0]
    else:
        u_str = 0
    tau_c = ru * u_str**2
    return tau_c


ks = 5.27e-4  # m
current_stress_3 = np.zeros(n, dtype=object)
for i in range(n):
    current_stress_3[i] = Current_bottom_shear_stress_3(
        0.05, 0.41, nu, ks, LO_Wd[i], ru
    )
Current_ShearStress_df["Current_Stress_3"] = (
    current_stress_3 * 10
)  # Convert N/m2 to Dyne/cm2
Current_ShearStress_df.to_csv("./Current_ShearStress.csv")
