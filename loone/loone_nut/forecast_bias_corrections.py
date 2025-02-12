# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 12:29:35 2024

@author: osama
"""
import pandas as pd
import numpy as np


# L001_Bias Corrected
obs_quantiles = [98, 122, 152]

correction_factor_1 = 0.5913678934015472
correction_factor_2 = 0.8353793753288427
correction_factor_3 = 1.0524001582697413
correction_factor_4 = 1.4605983806263183

Forecasted_TP_Lake_N = [100,90,65,150,168]
sim_BC_Qnt_EntPer = np.zeros(5)
Most_recent_TP_L001_Obs = 78
# Apply the correction factor to Forecasted TP_Lake_N
for i in range(len(Forecasted_TP_Lake_N)):
    if Most_recent_TP_L001_Obs <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_1 
    elif Most_recent_TP_L001_Obs > obs_quantiles[0] and Most_recent_TP_L001_Obs <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_2 
    elif Most_recent_TP_L001_Obs > obs_quantiles[1] and Most_recent_TP_L001_Obs <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_3
    elif Most_recent_TP_L001_Obs > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i]
        
print(sim_BC_Qnt_EntPer)


# L005_Bias Corrected
obs_quantiles = [56.25, 85, 115.25]

correction_factor_1 = 0.29828881937048657
correction_factor_2 = 0.5283271275855149
correction_factor_3 = 0.7632318707941915
correction_factor_4 = 1.188700366068547

Forecasted_TP_Lake_N = [100,90,65,150,168]
sim_BC_Qnt_EntPer = np.zeros(5)
Most_recent_TP_L005_Obs = 78
# Apply the correction factor to Forecasted TP_Lake_N
for i in range(len(Forecasted_TP_Lake_N)):
    if Most_recent_TP_L005_Obs <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_1 
    elif Most_recent_TP_L005_Obs > obs_quantiles[0] and Most_recent_TP_L005_Obs <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_2 
    elif Most_recent_TP_L005_Obs > obs_quantiles[1] and Most_recent_TP_L005_Obs <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_3
    elif Most_recent_TP_L005_Obs > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i]
        
print(sim_BC_Qnt_EntPer)

# L008_Bias Corrected
obs_quantiles = [98, 140.5, 174]

correction_factor_1 = 0.6062543848348515
correction_factor_2 = 0.9218176856485235
correction_factor_3 = 1.2001724608466868
correction_factor_4 = 1.7926071099770216

Forecasted_TP_Lake_N = [100,90,65,150,168]
sim_BC_Qnt_EntPer = np.zeros(5)
Most_recent_TP_L008_Obs = 78
# Apply the correction factor to Forecasted TP_Lake_N
for i in range(len(Forecasted_TP_Lake_N)):
    if Most_recent_TP_L008_Obs <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_1 
    elif Most_recent_TP_L008_Obs > obs_quantiles[0] and Most_recent_TP_L008_Obs <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_2 
    elif Most_recent_TP_L008_Obs > obs_quantiles[1] and Most_recent_TP_L008_Obs <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_3
    elif Most_recent_TP_L008_Obs > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_N[i]
        
print(sim_BC_Qnt_EntPer)

# L004_Bias Corrected
obs_quantiles = [122, 145, 177.5]

correction_factor_1 = 0.8955760993304536
correction_factor_2 = 1.1510189122660914
correction_factor_3 = 1.3689041532929282
correction_factor_4 = 1.9750371175582473

Forecasted_TP_Lake_S = [100,90,65,150,168]
sim_BC_Qnt_EntPer = np.zeros(5)
Most_recent_TP_L004_Obs = 78
# Apply the correction factor to Forecasted TP_Lake_N
for i in range(len(Forecasted_TP_Lake_N)):
    if Most_recent_TP_L004_Obs <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_1 
    elif Most_recent_TP_L004_Obs > obs_quantiles[0] and Most_recent_TP_L004_Obs <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_2 
    elif Most_recent_TP_L004_Obs > obs_quantiles[1] and Most_recent_TP_L004_Obs <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_3
    elif Most_recent_TP_L004_Obs > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i]
        
print(sim_BC_Qnt_EntPer)

# L006_Bias Corrected
obs_quantiles = [106.5, 135.5, 174]

correction_factor_1 = 0.7347463391896812
correction_factor_2 = 1.0422794056109943
correction_factor_3 = 1.3076781699380404
correction_factor_4 = 1.9930823718916406

Forecasted_TP_Lake_S = [100,90,65,150,168]
sim_BC_Qnt_EntPer = np.zeros(5)
Most_recent_TP_L006_Obs = 78
# Apply the correction factor to Forecasted TP_Lake_N
for i in range(len(Forecasted_TP_Lake_N)):
    if Most_recent_TP_L006_Obs <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_1 
    elif Most_recent_TP_L006_Obs > obs_quantiles[0] and Most_recent_TP_L006_Obs <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_2 
    elif Most_recent_TP_L006_Obs > obs_quantiles[1] and Most_recent_TP_L006_Obs <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_3
    elif Most_recent_TP_L006_Obs > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i]
        
print(sim_BC_Qnt_EntPer)

# L007_Bias Corrected
obs_quantiles = [86, 116, 151.125]

correction_factor_1 = 0.5927850757579236
correction_factor_2 = 0.9145550431423951
correction_factor_3 = 1.1594227707702824
correction_factor_4 = 1.6784629147274226

Forecasted_TP_Lake_S = [100,90,65,150,168]
sim_BC_Qnt_EntPer = np.zeros(5)
Most_recent_TP_L007_Obs = 78
# Apply the correction factor to Forecasted TP_Lake_N
for i in range(len(Forecasted_TP_Lake_N)):
    if Most_recent_TP_L007_Obs <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_1 
    elif Most_recent_TP_L007_Obs > obs_quantiles[0] and Most_recent_TP_L007_Obs <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_2 
    elif Most_recent_TP_L007_Obs > obs_quantiles[1] and Most_recent_TP_L007_Obs <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_3
    elif Most_recent_TP_L007_Obs > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i]
        
print(sim_BC_Qnt_EntPer)

# LZ40_Bias Corrected
obs_quantiles = [132.25, 158.75, 183.75]

correction_factor_1 = 0.9184276576796921
correction_factor_2 = 1.2443912322787118
correction_factor_3 = 1.4516529344806937
correction_factor_4 = 2.171551720078875

Forecasted_TP_Lake_S = [100,90,65,150,168]
sim_BC_Qnt_EntPer = np.zeros(5)
Most_recent_TP_LZ40_Obs = 78
# Apply the correction factor to Forecasted TP_Lake_N
for i in range(len(Forecasted_TP_Lake_N)):
    if Most_recent_TP_LZ40_Obs <= obs_quantiles[0]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_1 
    elif Most_recent_TP_LZ40_Obs > obs_quantiles[0] and Most_recent_TP_LZ40_Obs <= obs_quantiles[1]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_2 
    elif Most_recent_TP_LZ40_Obs > obs_quantiles[1] and Most_recent_TP_LZ40_Obs <= obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_3
    elif Most_recent_TP_LZ40_Obs > obs_quantiles[2]:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i] * correction_factor_4
    else:
        sim_BC_Qnt_EntPer[i] = Forecasted_TP_Lake_S[i]
        
print(sim_BC_Qnt_EntPer)
