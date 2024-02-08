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
Working_dir = 'C:/LOONE_NO_Clean' 
os.chdir('%s'%Working_dir) 
os.chdir('./Code/')
from LOONE_WQ import LOONE_NO
os.chdir('%s'%Working_dir) 

# Run NO Module and get monthly simulated NOx
NO_M = LOONE_NO()

#Read Observed NO "Monthly"
NO_Obs = pd.read_csv('./Model_Data_Filled_20082023/LO_NO_Obs20082023.csv')

slope, intercept, r_value, p_value, std_err = stats.linregress(NO_M['NO_M'], NO_Obs['NO'])
# Calculate R-squared (coefficient of determination)
r_squared = r_value**2
print("R2_NO_M:", r_squared)

NO_M['date'] = pd.to_datetime(NO_M['date'])
fig, ax = plt.subplots()
ax.plot(NO_M['date'], NO_M['NO_M'], color='blue',linestyle = 'dashed',label = 'Simulated')
ax.plot(NO_M['date'], NO_Obs['NO'], color='red',linestyle = 'solid',label = 'Observed')
ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.set_title('Mean NO in Lake Okeechobee')
ax.set_xlabel('date')
ax.set_ylabel('NO mg/m3')
ax.legend()
# plt.savefig('./OP N ts.png', dpi=600)
plt.show()


