# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 23:44:04 2022

@author: osama
"""

#This Script calculates water kinematic viscosity function of H2O_Temperature

import pandas as pd
import numpy as np
from scipy.optimize import fsolve
from scipy import log as log
import os
# Read Mean H2O_T in LO
Working_dir = 'C:/Work/Research/LOONE/LOONE_Data_Pre' 
os.chdir('%s'%Working_dir) 
LO_Temp = pd.read_csv('./Filled_WaterT_20082023.csv')
LO_T = LO_Temp['Water_T']
nu20 = 1.0034/1E6 # m2/s (kinematic viscosity of water at T = 20 C) 

n = len(LO_T.index)

class nu_Func:
    
    def nu(T):
        nu20 = 1.0034/1E6
        def func(x):
            # return[log(x[0]/nu20)-((20-T)/(T+96))*(1.2364-1.37E-3*(20-T)+5.7E-6*(20-T)**2)]
            return[(x[0]/nu20)-10**(((20-T)/(T+96))*(1.2364-1.37E-3*(20-T)+5.7E-6*(20-T)**2))]
        sol = fsolve(func,[9.70238995692062E-07])
        nu = sol[0]
        return(nu)


nu = np.zeros(n,dtype = object)

for i in range(n):
    nu[i] = nu_Func.nu(LO_T[i])


nu_df = pd.DataFrame(LO_Temp['date'],columns=['date'])
nu_df['nu'] = nu
Working_dir = 'C:/Work/Research/LOONE/Final_Data_20082023' 
os.chdir('%s'%Working_dir) 
nu_df.to_csv('./nu_20082023.csv')
