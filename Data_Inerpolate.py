# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:41:53 2023

@author: osamatarabih
"""
import math
from calendar import monthrange
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy import interpolate
import os
import datetime

parameter: str = 'RAD' #'RADP'
units: str = 'n' #'MICROMOLE/m^2/s' #
station: str = 'm' #'L006'#

os.chdir('C:/LOONE_WQ_Clean/Model_Data_Filled_20082023')
Data_In = pd.read_csv('./LO_RADT_20082023.csv')
Data_In = Data_In[Data_In['%s_%s_%s' % (station, parameter, units)]>=0]
Data_In = Data_In.set_index(['date'])
Data_In.index = pd.to_datetime(Data_In.index, unit='ns')
Data_df = Data_In.resample('D').mean()
Data_df = Data_df.dropna(subset=['%s_%s_%s' % (station, parameter, units)])
Data_df = Data_df.reset_index()
Data_df['Yr_M'] = pd.to_datetime(Data_df['date']).dt.to_period('M')
start_date = Data_df['date'].iloc[0]
end_date = Data_df['date'].iloc[-1]
date_rng = pd.date_range(start=start_date, end=end_date, freq='M')
date_rng = date_rng.union([date_rng[-1] + pd.DateOffset(months=1)])
Monthly_df = pd.DataFrame(date_rng, columns=['date'])
Monthly_df['Yr_M'] = pd.to_datetime(Monthly_df['date']).dt.to_period('M')
New_date = []
New_data = []
Days = []
Days_cum = []
# Set index for the two dataframes
Data_df = Data_df.set_index(['Yr_M'])
Monthly_df = Monthly_df.set_index(['Yr_M'])
for i in Monthly_df.index:
    if i in Data_df.index:
        if type(Data_df.loc[i]['date']) == pd.Timestamp:
            New_date.append(Data_df.loc[i]['date'])
            New_data.append(Data_df.loc[i]['%s_%s_%s' % (station, parameter, units)])
        else:
            for j in range(len(Data_df.loc[i]['date'])):
                New_date.append(Data_df.loc[i]['date'][j])
                New_data.append(Data_df.loc[i]['%s_%s_%s' % (station, parameter, units)][j])
    elif i not in Data_df.index:
        New_date.append(datetime.datetime(Monthly_df.loc[i]['date'].year, Monthly_df.loc[i]['date'].month, 1))
        New_data.append(np.NaN)

New_date = pd.to_datetime(New_date, format='%Y-%m-%d')
Days = New_date.strftime("%d").astype(float)
for i in range(len(Days)):
    if i == 0:
        Days_cum.append(Days[i])
    elif New_date[i].month == New_date[i-1].month:
        Days_cum.append(Days_cum[i-1]+(Days[i]-Days[i-1]))
    elif New_date[i].month != New_date[i-1].month:
        Days_cum.append(Days_cum[i-1]+Days[i]+monthrange(New_date[i-1].year, New_date[i-1].month)[1]-Days[i-1])
Final_df = pd.DataFrame()
Final_df['date'] = New_date
Final_df['Data'] = New_data
Final_df['Days'] = Days
Final_df['Days_cum'] = Days_cum
# Final_df.to_csv('C:/Work/Research/LOONE/Nitrogen Module/Interpolated_Data/In-Lake/L008_DO_No_Months_Missing_Trial.csv')  # noqa: E501
# Remove Negative Data Values
Final_df = Final_df[Final_df['Data'] >= 0]
Final_df['date'] = pd.to_datetime(Final_df['date'], format='%Y-%m-%d')
start_date = Final_df['date'].iloc[0]
end_date = Final_df['date'].iloc[-1]
date_rng_TSS_1 = pd.date_range(start=start_date, end=end_date, freq='D')
# Create a data frame with a date column
Data_df = pd.DataFrame(date_rng_TSS_1, columns=['date'])
Data_len = len(Data_df.index)
Cum_days = np.zeros(Data_len)
Data_daily = np.zeros(Data_len)
# Set initial values
Cum_days[0] = Data_df['date'].iloc[0].day
Data_daily[0] = Final_df['Data'].iloc[0]
for i in range(1, Data_len):
    Cum_days[i] = Cum_days[i-1]+1
    # Data_daily[i] = interpolate.interp1d(Final_df['Days'], Final_df['TSS'] , kind = 'linear')(Cum_days[i])
    Data_daily[i] = np.interp(Cum_days[i], Final_df['Days_cum'], Final_df['Data'])
Data_df['Data_%s'%station] = Data_daily
Data_df.to_csv('./LO_RADT_20082023.csv', index=False)