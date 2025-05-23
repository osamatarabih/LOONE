# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:31:41 2022

@author: osamatarabih
"""

#This Script incorporates the Comprehensive LOONE Model!
import os
import pandas as pd
import numpy as np
from Model_Config import Model_Config
Working_Path = Model_Config.Working_Path
os.chdir('./Packages/Platypus-1.0.4')   
from platypus import NSGAII, Problem, Real, nondominated
os.chdir('%s'%Working_Path) 
from LOONE_Q_4All import LOONE_Q
os.chdir('%s/Code'%Working_Path) 
from LOONE_WQ import LOONE_Constituent_SimQ
from LOONE_Nut_3M import LOONE_Nut_NS
import matplotlib.pyplot as plt


############### Optimization for combined constituents
def Opt(vars):
    P_1 = vars[0]
    P_2 = vars[1] 
    NO_1 = vars[2]
    NO_2 = vars[3] 
    Chla_1 = vars[4]
    Chla_2 = vars[5] 
    S77_DV = vars[6:18]    
    S308_DV = vars[18:30]
    RegSouth_DV = vars[30:42]
    LOONE_Q_Outputs_df = pd.read_csv('./Outputs/LOONE_Q_Outputs.csv')
    NO_Chla_Out = LOONE_Constituent_SimQ(LOONE_Q_Outputs_df)

    # NOx_StL_Lds_out_df = pd.DataFrame(NO_Chla_Out[0])
    # Chla_StL_Lds_out_df = pd.DataFrame(NO_Chla_Out[1])
    NOx_Cal_Lds_out_df = pd.DataFrame(NO_Chla_Out[2])
    Chla_Cal_Lds_out_df = pd.DataFrame(NO_Chla_Out[3])
    NO_Chla_df = pd.DataFrame(NO_Chla_Out[4])
    
    TP_Out = LOONE_Nut_NS(LOONE_Q_Outputs_df)
    # P_StL_Lds_df = TP_Out[0]
    P_Cal_Lds_df = TP_Out[1]
    P_df = TP_Out[2]

    LOONE_Q_Outputs = LOONE_Q(P_1,P_2,NO_1,NO_2,Chla_1,Chla_2,S77_DV,S308_DV,RegSouth_DV,P_df['TP_Lake_S'],NO_Chla_df['NO_S'],NO_Chla_df['Sim_Chla_S'])

    # LOONE_Q_Outputs = LOONE_Q(P_1,P_2,NO_1,NO_2,Chla_1,Chla_2,S77_DV,S308_DV,P_df['TP_Lake_S'],NO_Chla_df['NO_S'],NO_Chla_df['Sim_Chla_S'])

    LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
    LOONE_Q_Outputs_df.to_csv('./Outputs/LOONE_Q_Outputs.csv')
    return (LOONE_Q_Outputs_df['S308_Q'].sum(),NOx_Cal_Lds_out_df[0].mean(),Chla_Cal_Lds_out_df[0].mean(),P_Cal_Lds_df[0].mean(),LOONE_Q_Outputs_df['Cutback'].sum())


Var_value = []
Var_value.append(Real(150,400)) 
Var_value.append(Real(50,100)) 
Var_value.append(Real(300,500)) 
Var_value.append(Real(100,300)) 
Var_value.append(Real(10,80)) 
Var_value.append(Real(10,40)) 

#S77
Var_value.append(Real(0,1200)) #1
Var_value.append(Real(0,1200)) #2
Var_value.append(Real(0,1200)) #3
Var_value.append(Real(0,1200)) #4
Var_value.append(Real(0,900)) #5
Var_value.append(Real(0,900)) #6
Var_value.append(Real(0,900)) #7
Var_value.append(Real(0,900)) #8
Var_value.append(Real(0,900)) #9
Var_value.append(Real(0,900)) #10
Var_value.append(Real(0,1200)) #11
Var_value.append(Real(0,1200)) #12
#S308
Var_value.append(Real(0,400)) #1
Var_value.append(Real(0,400)) #2
Var_value.append(Real(0,400)) #3
Var_value.append(Real(0,400)) #4
Var_value.append(Real(0,300)) #5
Var_value.append(Real(0,300)) #6
Var_value.append(Real(0,300)) #7
Var_value.append(Real(0,300)) #8
Var_value.append(Real(0,300)) #9
Var_value.append(Real(0,300)) #10
Var_value.append(Real(0,400)) #11
Var_value.append(Real(0,400)) #12

# #RegSouth
Var_value.append(Real(1000,2000)) #1
Var_value.append(Real(1000,2000)) #2
Var_value.append(Real(1600,3200)) #3
Var_value.append(Real(1400,2800)) #4
Var_value.append(Real(1200,2400)) #5
Var_value.append(Real(500,1000)) #6
Var_value.append(Real(300,600)) #7
Var_value.append(Real(500,1000)) #8
Var_value.append(Real(400,800)) #9
Var_value.append(Real(800,1600)) #10
Var_value.append(Real(800,1600)) #11
Var_value.append(Real(900,1800)) #12


problem = Problem(42,5,0)

#Set the upper amd lower limits of decision variables
problem.types[:] = Var_value
#Set the Constraint value 
# problem.constraints[0:(n_rows-2)] = ">=0" #Lake Okeechobee stage must be greater than Minimum stage(9.3)!
# problem.constraints[(n_rows-2)+1:(n_rows-2)*2] = "<=0"  #Lake Okeechobee stage must be less than Maximum stage(15.2)!
problem.function = Opt
algorithm = NSGAII(problem)
algorithm.run(1000)
results = algorithm.result
feasible_solutions = [s for s in results if s.feasible]
nondominated_solutions = nondominated(results)


np.savetxt("./Outputs/optimization_objectives_feasible_LOSOM_Scen4_StLMinQ.txt",[s.objectives[:] for s in feasible_solutions],fmt="%s")
np.savetxt("./Outputs/optimization_variables_feasible_LOSOM_Scen4_StLMinQ.txt",[s.variables[:] for s in feasible_solutions],fmt="%s")
np.savetxt("./Outputs/optimization_objectives_nondominated_LOSOM_Scen4_StLMinQ.txt",[s.objectives[:] for s in nondominated_solutions],fmt="%s")
np.savetxt("./Outputs/optimization_variables_nondominated_LOSOM_Scen4_StLMinQ.txt",[s.variables[:] for s in nondominated_solutions],fmt="%s")