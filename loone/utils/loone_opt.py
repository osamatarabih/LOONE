# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:31:41 2022

@author: osamatarabih
"""

# This Script incorporates the Comprehensive LOONE Model!
import os
import pandas as pd
import numpy as np
from loone_config.Model_Config import Model_Config

Working_Path = Model_Config.Working_Path
# os.chdir("../Packages/Platypus-1.0.4")
from platypus import NSGAII, Problem, Real, nondominated

os.chdir("%s" % Working_Path)
from loone_q.LOONE_Q import LOONE_Q
from loone.loone_nut.LOONE_NUT import LOONE_Nut
import matplotlib.pyplot as plt


def Opt(vars):
    P_1 = vars[0]
    P_2 = vars[1]
    S77_DV = vars[2:14]
    S308_DV = vars[14:26]
    LOONE_Q_Outputs_df = pd.read_csv("./Outputs/LOONE_Q_Outputs.csv")
    LOONE_Nut_out = LOONE_Nut(LOONE_Q_Outputs_df)
    StL_Lds_out_df = pd.DataFrame(LOONE_Nut_out[0])
    Calo_Lds_out_df = pd.DataFrame(LOONE_Nut_out[1])
    P_Lake_df = pd.DataFrame(LOONE_Nut_out[2])

    LOONE_Q_Outputs = LOONE_Q(P_1, P_2, S77_DV, S308_DV, P_Lake_df["TP_Lake_S"])
    LOONE_Q_Outputs_df = pd.DataFrame(LOONE_Q_Outputs[0])
    LOONE_Q_Outputs_df.to_csv("./Outputs/LOONE_Q_Outputs.csv")
    return (
        StL_Lds_out_df[0].mean(),
        Calo_Lds_out_df[0].mean(),
        LOONE_Q_Outputs_df["Cutback"].sum(),
    )


Var_value = []
Var_value.append(Real(150, 400))  # 1
Var_value.append(Real(50, 100))  # 2
# S77
Var_value.append(Real(0, 900))  # 1
Var_value.append(Real(0, 900))  # 2
Var_value.append(Real(0, 900))  # 3
Var_value.append(Real(0, 900))  # 4
Var_value.append(Real(0, 900))  # 5
Var_value.append(Real(0, 900))  # 6
Var_value.append(Real(0, 900))  # 7
Var_value.append(Real(0, 900))  # 8
Var_value.append(Real(0, 900))  # 9
Var_value.append(Real(0, 900))  # 10
Var_value.append(Real(0, 900))  # 11
Var_value.append(Real(0, 900))  # 12
# S308
Var_value.append(Real(0, 250))  # 1
Var_value.append(Real(0, 250))  # 2
Var_value.append(Real(0, 250))  # 3
Var_value.append(Real(0, 250))  # 4
Var_value.append(Real(0, 400))  # 5
Var_value.append(Real(0, 400))  # 6
Var_value.append(Real(0, 400))  # 7
Var_value.append(Real(0, 400))  # 8
Var_value.append(Real(0, 400))  # 9
Var_value.append(Real(0, 400))  # 10
Var_value.append(Real(0, 250))  # 11
Var_value.append(Real(0, 250))  # 12
# for i in range(2*(n_rows-2)):
#     Var_value.append(Real(0,1200))
# problem = Problem((n_rows-2)*2,2,(n_rows-2)*3)
problem = Problem(26, 3, 0)

# Set the upper amd lower limits of decision variables
problem.types[:] = Var_value
# Set the Constraint value
# problem.constraints[0:(n_rows-2)] = ">=0" #Lake Okeechobee stage must be greater than Minimum stage(9.3)!
# problem.constraints[(n_rows-2)+1:(n_rows-2)*2] = "<=0"  #Lake Okeechobee stage must be less than Maximum stage(15.2)!
problem.function = Opt
algorithm = NSGAII(problem)
algorithm.run(5)
results = algorithm.result
feasible_solutions = [s for s in results if s.feasible]
nondominated_solutions = nondominated(results)


np.savetxt(
    "./Outputs/optimization_objectives_feasible_LORS2008_0529.txt",
    [s.objectives[:] for s in feasible_solutions],
    fmt="%s",
)
np.savetxt(
    "./Outputs/optimization_variables_feasible_LORS2008_0529.txt",
    [s.variables[:] for s in feasible_solutions],
    fmt="%s",
)
np.savetxt(
    "./Outputs/optimization_objectives_nondominated_LORS2008_0529.txt",
    [s.objectives[:] for s in nondominated_solutions],
    fmt="%s",
)
np.savetxt(
    "./Outputs/optimization_variables_nondominated_LORS2008_0529.txt",
    [s.variables[:] for s in nondominated_solutions],
    fmt="%s",
)
###################/##############################################################################################
# Read Opt Results (Feasible and Non-dominant Solutions)
Feasible_Obj = pd.read_csv(
    "./Outputs/optimization_objectives_feasible_LORS2008_0529.csv"
)
Feasible_Var = pd.read_csv(
    "./Outputs/optimization_variables_feasible_LORS2008_0529.csv"
)
Nondominant_Obj = pd.read_csv(
    "./Outputs/optimization_objectives_nondominated_LORS2008_0529.csv"
)
Nondominant_Var = pd.read_csv(
    "./Outputs/optimization_variables_nondominated_LORS2008_0529.csv"
)

# Open operation file
totalobs = 0
objs = [0 for x in range(3)]
labels = []
for i in range(3):
    objs[i] = int(i + 1)
    if objs[i] >= 1:
        totalobs += 1
    if i == 0:
        labels.append("StL P")
    elif i == 1:
        labels.append("Cal P")
    elif i == 2:
        labels.append("Water deficit")

# Plot data 2D
if totalobs >= 2:
    for i in range(totalobs - 1):
        for j in range(i + 1, totalobs):
            x = np.array(Nondominant_Obj)[:, i]
            y = np.array(Nondominant_Obj)[:, j]
            plt.subplots(figsize=(10, 10))
            plt.grid(alpha=0.8, c="gray")
            plt.scatter(
                x, y, s=100, c="g", alpha=0.6667, edgecolors="black", linewidths=0.6667
            )
            plt.xlabel(labels[i], size=18, labelpad=5)
            plt.ylabel(labels[j], size=18, labelpad=5)
            plt.tick_params(axis="x", labelsize=18, pad=5)
            plt.tick_params(axis="y", labelsize=18, pad=5)
            plt.show()

# Plot data 3D
if totalobs >= 3:
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection="3d")
    for i in range(totalobs - 2):
        for j in range(i + 1, totalobs - 1):
            for k in range(i + 2, totalobs):
                x = np.array(Feasible_Obj)[:, i]
                y = np.array(Feasible_Obj)[:, j]
                z = np.array(Feasible_Obj)[:, k]
                ax.scatter(
                    x,
                    y,
                    z,
                    s=100,
                    c="g",
                    alpha=0.6667,
                    edgecolors="black",
                    linewidths=0.6667,
                )
                ax.view_init(20, 20)
                ax.tick_params(axis="x", labelsize=18, pad=5)
                ax.tick_params(axis="y", labelsize=18, pad=5)
                ax.tick_params(axis="z", labelsize=18, pad=5)
                ax.set_xlabel(labels[i], size=18, labelpad=20)
                ax.set_ylabel(labels[j], size=18, labelpad=20)
                ax.set_zlabel(labels[k], size=18, labelpad=20)
                plt.show()
