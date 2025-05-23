# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 00:32:55 2023

@author: osamatarabih
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.stats as stats
import os
Working_dir = 'C:/LOONE_WQ_Combo' 
os.chdir('%s'%Working_dir) 
os.chdir('./Code/')
from LOONE_WQ import LOONE_Constituent_SimQ
os.chdir('%s'%Working_dir) 

# Export NO Daily Simulations
# NO_M = LOONE_NO()
# NO_M.to_csv('./NO_daily_sim.csv')

# Monthly Calibrations
NO_M = LOONE_Constituent_SimQ()
NO_M[1].to_csv('Nutrient_combined.csv')
#Read Observed NO "Monthly"
# NO_Obs = pd.read_csv('./Model_Data_Filled_20082023/LO_NO_Obs20082023.csv')
# slope, intercept, r_value, p_value, std_err = stats.linregress(NO_M['NO_M'], NO_Obs['Mean_NO'])
# # Calculate R-squared (coefficient of determination)
# r_squared = r_value**2
# print("R2_NO_M:", r_squared)

# NO_M['date'] = pd.to_datetime(NO_M['date'])
# fig, ax = plt.subplots()
# ax.plot(NO_M['date'], NO_M['NO_M'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(NO_M['date'], NO_Obs['Mean_NO'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean NO in Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('NO mg/m3')
# ax.legend()
# # plt.savefig('./OP N ts.png', dpi=600)
# plt.show()
# NO_M.to_csv('./NOx_Mar24.csv')


#########################################################################
# NH4_M = LOONE_NH4()
# #Read Observed NH4 "Monthly"
# NH4_Obs = pd.read_csv('./Model_Data_Filled_20082023/LO_NH4_M_Obs.csv')

# slope, intercept, r_value, p_value, std_err = stats.linregress(NH4_M['NH4_M'], NH4_Obs['NH4'])
# # Calculate R-squared (coefficient of determination)
# r_squared = r_value**2
# print("R2_NH4_M:", r_squared)

# NH4_M['date'] = pd.to_datetime(NH4_M['date'])
# fig, ax = plt.subplots()
# ax.plot(NH4_M['date'], NH4_M['NH4_M'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(NH4_M['date'], NH4_Obs['NH4'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean NH4 in Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('NH4 mg/m3')
# ax.legend()
# # plt.savefig('./OP N ts.png', dpi=600)
# plt.show()



#################################################################
# Chla_M = LOONE_Chla()
# #Read Observed Chla "Monthly"
# Chla_Obs = pd.read_csv('./Model_Data_Filled_20082023/LO_Chla_Monthly_M.csv')

# slope, intercept, r_value, p_value, std_err = stats.linregress(Chla_M['Sim_Chla'], Chla_Obs['Chla'])
# # Calculate R-squared (coefficient of determination)
# r_squared = r_value**2
# print("R2_Chla_M:", r_squared)

# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs['Chla'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_Mean_alt1.png', dpi=600)
# plt.show()
# Chla_M.to_csv('./Sim_Chla_Mar24.csv')



# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla_N'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs['Chla_N'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in North Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_N_alt1.png', dpi=600)
# plt.show()

# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla_S'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs['Chla_S'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in South Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_S_alt1.png', dpi=600)
# plt.show()
# Chla_M.to_csv('./Chla_Sim_Dec23.csv')
# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['External_Chla'], color='blue',linestyle = 'dashed',label = 'External_Chla')
# ax.plot(Chla_M['date'], Chla_M['Chla Load N2S'], color='red',linestyle = 'solid',label = 'Load N2S')
# ax.plot(Chla_M['date'], Chla_M['Chla Load Out'], color='green',linestyle = 'dashdot',label = 'Load Out')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Chla Loads')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla kg')
# ax.legend()
# # plt.savefig('./Chla_Sim4.png', dpi=600)
# plt.show()

# print(Chla_M['External_Chla'].mean())
# print(Chla_M['Chla Load N2S'].mean())

###############################################################
# Chla_M = LOONE_Chla_Temp()
# #Read Observed Chla "Monthly"
# Chla_Obs = pd.read_csv('./Model_Data_Filled_20082023/Obs_Chla_LO_2008-2022.csv')
# Chla_Obs_N = pd.read_csv('./Model_Data_Filled_20082023/Obs_Chla_M_N.csv')
# Chla_Obs_S = pd.read_csv('./Model_Data_Filled_20082023/Obs_Chla_M_S.csv')

# slope, intercept, r_value, p_value, std_err = stats.linregress(Chla_M['Sim_Chla'], Chla_Obs['Chla'])
# # Calculate R-squared (coefficient of determination)
# r_squared = r_value**2
# print("R2_Chla_M:", r_squared)

# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs['Chla'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_Sim4.png', dpi=600)
# plt.show()

# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla_N'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs_N['Chla_N'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in North Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_Sim4.png', dpi=600)
# plt.show()

# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla_S'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs_S['Chla_S'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in South Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_Sim4.png', dpi=600)
# plt.show()


###############################################################
# Chla_M = LOONE_Chla_TempMix()
# #Read Observed Chla "Monthly"
# Chla_Obs = pd.read_csv('./Model_Data_Filled_20082023/Obs_Chla_LO_2008-2022.csv')
# Chla_Obs_N = pd.read_csv('./Model_Data_Filled_20082023/Obs_Chla_M_N.csv')
# Chla_Obs_S = pd.read_csv('./Model_Data_Filled_20082023/Obs_Chla_M_S.csv')

# slope, intercept, r_value, p_value, std_err = stats.linregress(Chla_M['Sim_Chla'], Chla_Obs['Chla'])
# # Calculate R-squared (coefficient of determination)
# r_squared = r_value**2
# print("R2_Chla_M:", r_squared)

# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs['Chla'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_Sim3.png', dpi=600)
# plt.show()
# # Chla_M.to_csv('./Sim_Chla_Nov2023.csv')



# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla_N'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs_N['Chla_N'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in North Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_Sim4.png', dpi=600)
# plt.show()

# Chla_M['date'] = pd.to_datetime(Chla_M['date'])
# fig, ax = plt.subplots()
# ax.plot(Chla_M['date'], Chla_M['Sim_Chla_S'], color='blue',linestyle = 'dashed',label = 'Simulated')
# ax.plot(Chla_M['date'], Chla_Obs_S['Chla_S'], color='red',linestyle = 'solid',label = 'Observed')
# ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
# ax.set_title('Mean Chla in South Lake Okeechobee')
# ax.set_xlabel('date')
# ax.set_ylabel('Chla mg/m3')
# ax.legend()
# # plt.savefig('./Chla_Sim4.png', dpi=600)
# plt.show()