# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:35:15 2022

@author: osamatarabih
"""

import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.stats as stats
import pandas as pd
from Model_Config import Model_Config 
Working_Path = Model_Config.Working_Path
os.chdir('%s'%Working_Path) 
# from LOONE_Q import LOONE_Q
from LOONE_Q_4NO import LOONE_Q
from LOONE_Nut_3M import LOONE_Nut_NS_3M, LOONE_Nut_NS
os.chdir('./Code/')
from LOONE_WQ import LOONE_NO_SimQ, LOONE_NO, LOONE_Chla, LOONE_Chla_SimQ
os.chdir('%s'%Working_Path) 

# from LOONE_Q_Sims import LOONE_Q_Sims
# from LOONE_Nut import LOONE_Nut
# Read Decision Variable Values out ot Optimization
# Dec_Var = pd.read_csv('./Outputs/Opt_Decision_Var.csv')
# Dec_Var_Val = Dec_Var['Value']
# LOONE_Q_Outputs = LOONE_Q(Dec_Var_Val[0],Dec_Var_Val[1],Dec_Var_Val[2:14],Dec_Var_Val[14:26],P_Lake_df['TP_Lake_S'])

LOONE_Q_Outputs = LOONE_Q(0,0,0,0,0) 
# LOONE_Q_Outputs[0].to_csv('./LOONE_Q_Outputs.csv')  

# LOONE_P = LOONE_Nut_NS_3M()
# LOONE_P.to_csv('./P_Outputs.csv')
# Q_Outputs = pd.read_csv('./P_Mod_Data_LORS20082023_3M/Outputs1_20072023.csv')
# LOONE_P = LOONE_Nut_NS(Q_Outputs)
# LOONE_P.to_csv('./P_Outputs_QSim1.csv')

# LOONE_NO_Outputs = LOONE_NO_SimQ(LOONE_Q_Outputs[0])
# LOONE_NO_Outputs[0].to_csv('./LOONE_NO_Loads_SimQ(LOSOM).csv')

# Q_Outputs = pd.read_csv('./Model_Data_Filled_20082023/Outputs1_20072023.csv')
# LOONE_NO_Outputs = LOONE_NO_SimQ(Q_Outputs)
# LOONE_NO_Outputs.to_csv('./LOONE_NO_Outputs_SimQ.csv')
# LOONE_NO_Outputs_df = pd.DataFrame(LOONE_NO_Outputs[0])
# LOONE_NO_Outputs_df.to_csv('./LOONE_NO_Loads_Sim_2.csv')
# LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
# LOONE_Q_Outputs_df.to_csv('./LOONE_Q_Nov2023.csv')

### Chla
# Chla_Outputs = LOONE_Chla()
# Chla_Outputs.to_csv('./Chla_ObsQ.csv')

# Q_Outputs = pd.read_csv('./Model_Data_Filled_20082023/Outputs1_20072023.csv')
Chla_Outputs = LOONE_Chla_SimQ(LOONE_Q_Outputs[0])
Chla_Outputs[0].to_csv('./Chla_SimQ(LORS).csv')


#### Loads
# Q_Outputs = pd.read_csv('./Model_Data_Filled_20082023/Outputs1_20072023.csv')
# LOONE_NO_Outputs = LOONE_NO_SimQ(Q_Outputs)
# LOONE_NO_Outputs[0].to_csv('./LOONE_NO_Loads_SimQ.csv')

#### TP Loads
# Q_Outputs = pd.read_csv('./P_Mod_Data_LORS20082023_3M/Outputs1_20072023.csv')
# LOONE_TP_Outputs = LOONE_Nut_NS(Q_Outputs)
# LOONE_TP_Outputs[0].to_csv('./LOONE_TP_Loads_SimQ.csv')


# Q_Outputs = pd.read_csv('./Model_Data_Filled_20082023/Outputs1_20072023.csv')
# LOONE_Chla_Outputs = LOONE_Chla_SimQ(Q_Outputs)
# LOONE_Chla_Outputs[0].to_csv('./LOONE_Chla_Loads_SimQ.csv')



### NO Opt Validation
# Read Decision Variable Values out ot Optimization
# Dec_Var = pd.read_csv('./Outputs/Dec_Var.csv')
# Dec_Var_Val = Dec_Var['Value']
# LOONE_Q_Outputs = LOONE_Q(Dec_Var_Val[0],Dec_Var_Val[1],Dec_Var_Val[2:14],Dec_Var_Val[14:26],P_Lake_df['TP_Lake_S'])


#S77_D1
# S77_RegRelRates_D1 = [0,1000,2000,3000]
# for i in range(len(S77_RegRelRates_D1)):
#     LOONE_Q_Outputs = LOONE_Q_Sims(0,0,0,0,0,S77_RegRelRates_D1[i])
#     LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
#     LOONE_Q_Outputs_df.to_csv('./Outputs/LOONE_Q_Sims_S77_RegRelRates_D1_%s.csv'%S77_RegRelRates_D1[i])
    
# #S77_D0    
# S77_RegRelRates_D0 = [0,100,200,300,400,450,500,600,700]
# for i in range(len(S77_RegRelRates_D0)):
#     LOONE_Q_Outputs = LOONE_Q_Sims(0,0,0,0,0,S77_RegRelRates_D0[i])
#     LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
#     LOONE_Q_Outputs_df.to_csv('./Outputs/LOONE_Q_Sims_S77_RegRelRates_D0_%s.csv'%S77_RegRelRates_D0[i])
  
# #MaxQstgTrigger    
# MaxQstgTrigger = [15,16,17,18,19,20]
# for i in range(len(MaxQstgTrigger)):
#     LOONE_Q_Outputs = LOONE_Q_Sims(0,0,0,0,0,MaxQstgTrigger[i])
#     LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
#     LOONE_Q_Outputs_df.to_csv('./Outputs/LOONE_Q_Sims_MaxQstgTrigger_%s.csv'%MaxQstgTrigger[i])
  
# S80_RegRelRates_Zone_A = [3000,4000,5000,6000,7200,8000,9000,10000]
# for i in range(len(S80_RegRelRates_Zone_A)):
#     LOONE_Q_Outputs = LOONE_Q_Sims(0,0,0,0,0,S80_RegRelRates_Zone_A[i])
#     LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
#     LOONE_Q_Outputs_df.to_csv('./Outputs/LOONE_Q_Sims_S80_RegRelRates_Zone_A_%s.csv'%S80_RegRelRates_Zone_A[i])    
# LOONE_Nut_out = LOONE_Nut(LOONE_Q_Outputs_df)
# LOONE_Nut_Lds_out_df = pd.DataFrame(LOONE_Nut_out)
# LOONE_Nut_Lds_out_df.to_csv('./Outputs/LOONE_Nut_Lds_Outputs_2023.csv')
# LOONE_Nut_Lake_out_df = pd.DataFrame(LOONE_Nut_out[1])
# LOONE_Nut_Lake_out_df.to_csv('./Outputs/LOONE_Nut_Lake_Outputs.csv')
# LOONE_Nut_Smr_StL_df = pd.DataFrame(LOONE_Nut_out[2])
# LOONE_Nut_Smr_StL_df.to_csv('./Outputs/LOONE_Nut_Smr_StL.csv')
# LOONE_Nut_Smr_Calo_df = pd.DataFrame(LOONE_Nut_out[3])
# LOONE_Nut_Smr_Calo_df.to_csv('./Outputs/LOONE_Nut_Smr_Calo.csv')
