# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:35:15 2022

@author: osamatarabih
"""

import os
import pandas as pd
from Model_Config import Model_Config 
Working_Path = Model_Config.Working_Path
os.chdir('%s'%Working_Path) 
from LOONE_Q import LOONE_Q
from LOONE_Nut import LOONE_Nut
# Read Decision Variable Values out ot Optimization
# Dec_Var = pd.read_csv('./Outputs/Opt_Decision_Var.csv')
# Dec_Var_Val = Dec_Var['Value']
# LOONE_Q_Outputs = LOONE_Q(Dec_Var_Val[0],Dec_Var_Val[1],Dec_Var_Val[2:14],Dec_Var_Val[14:26],P_Lake_df['TP_Lake_S'])

# LOONE_Q_Outputs = LOONE_Q(0,0,0,0,0)
# LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
# LOONE_Q_Outputs_df.to_csv('./Outputs/LOONE_Q_Outputs.csv')

# LOONE_Nut_out = LOONE_Nut(LOONE_Q_Outputs_df)
# LOONE_Nut_Lds_out_df = pd.DataFrame(LOONE_Nut_out[0])
# LOONE_Nut_Lds_out_df.to_csv('./Outputs/LOONE_Nut_Lds_Outputs.csv')
# LOONE_Nut_Lake_out_df = pd.DataFrame(LOONE_Nut_out[1])
# LOONE_Nut_Lake_out_df.to_csv('./Outputs/LOONE_Nut_Lake_Outputs.csv')
# LOONE_Nut_Smr_StL_df = pd.DataFrame(LOONE_Nut_out[2])
# LOONE_Nut_Smr_StL_df.to_csv('./Outputs/LOONE_Nut_Smr_StL.csv')
# LOONE_Nut_Smr_Calo_df = pd.DataFrame(LOONE_Nut_out[3])
# LOONE_Nut_Smr_Calo_df.to_csv('./Outputs/LOONE_Nut_Smr_Calo.csv')
