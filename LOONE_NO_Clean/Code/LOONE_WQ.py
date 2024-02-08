# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 00:32:29 2023

@author: osamatarabih
"""

#Daily Nitrogen Modeling

#import required packages
import pandas as pd
import numpy as np
from scipy.integrate import odeint
import os
Working_dir = 'C:/LOONE_NO_Clean' 
os.chdir('%s/Code/'%Working_dir) 
from LOONE_NCHLA_FNS import *
os.chdir('%s'%Working_dir) 

#### NO Simulations given Observed NH4 and observed Outflows
def LOONE_NO():
    #Read Required Data
    Q_in = pd.read_csv('./Model_Data_Filled_20082023/LO_Inflows_20082023.csv')
    Q_out = pd.read_csv('./Model_Data_Filled_20082023/LO_Outflows_20082023.csv')
    Temp_data = pd.read_csv('./Model_Data_Filled_20082023/Temp_Avg_20082023.csv')
    DO_data = pd.read_csv('./Model_Data_Filled_20082023/LO_Avg_DO_20082023.csv')
    RAD_data = pd.read_csv('./Model_Data_Filled_20082023/LO_RADT_20082023.csv')
    Storage = pd.read_csv('./Model_Data_Filled_20082023/Average_LO_Storage_20082023.csv')
    Storage_deviation = pd.read_csv('./Model_Data_Filled_20082023/Storage_Dev_20082023.csv')
    Chla_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_Merged_Chla.csv') #microgram/L
    Chla_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_Merged_Chla.csv') #microgram/L
    NOx_In = pd.read_csv('./Model_Data_Filled_20082023/Daily_NOx_External_Loads_20082023.csv') #mg
    S65E_NO_data = pd.read_csv('./Model_Data_Filled_20082023/S65E_NO_Interpolated.csv') #mg/m3
    Photoperiod = pd.read_csv('./Model_Data_Filled_20082023/PhotoPeriod_20082023.csv') 
    Photoperiod['Date'] = pd.to_datetime(Photoperiod['Date'])
    LO_DIP_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_OP.csv') #mg/m3
    LO_DIP_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_OP.csv') #mg/m3
    LO_DIN_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_DIN.csv') #mg/m3
    LO_DIN_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_DIN.csv') #mg/m3
    # RAD = RAD_data['RAD'].astype(float) * 0.484583*1000 #This one gives the best results for NO simulations
    RAD = RAD_data['RAD'].astype(float) *4.6*1000 

    Temp = Temp_data['Mean_T'].astype(float)
    DO = DO_data['DO'].astype(float)
    External_NO = NOx_In['NOx_Load_mg'].astype(float) #mg
    S65E_NO = S65E_NO_data['NO'].astype(float) #mg/m3
    
    # N-S Procedure
    N_Per = 0.43
    S_Per = 0.57
    
    volume  = Storage['Storage_m3'].astype(float) #m3
    volume_N = volume*N_Per
    volume_S = volume*S_Per
    Storage_dev = Storage_deviation['DS_dev'].astype(float) #acft
    Stage = Storage['Stage_ft'].astype(float) * 0.3048 #m
    Q_I = Q_in['Inflows_cmd'].astype(float) #m3
    Q_O = Q_out['LO_Outflows_cmd'].astype(float) #m3
    NH4_N_Obs = LO_DIN_N_data['NH4'].astype(float) #mg/m3
    NH4_S_Obs = LO_DIN_S_data['NH4'].astype(float) #mg/m3
    NH4_Obs = (NH4_N_Obs+NH4_S_Obs)/2
    NOx_N_Obs = LO_DIN_N_data['NO'].astype(float) #mg/m3
    NOx_S_Obs = LO_DIN_S_data['NO'].astype(float) #mg/m3
    NOx_Obs = (NOx_N_Obs+NOx_S_Obs)/2
    Chla_N = Chla_N_data['Chla'].astype(float) #microgram/L = mg/m3
    Chla_S = Chla_S_data['Chla'].astype(float) #microgram/L = mg/m3
    Chla = (Chla_N+Chla_S)/2 
    DIN_N_Obs = LO_DIN_N_data['DIN'].astype(float)
    DIN_S_Obs = LO_DIN_S_data['DIN'].astype(float)
    DIN_Obs = (DIN_N_Obs+DIN_S_Obs)/2
    DIP_N_Obs = LO_DIP_N_data['OP'].astype(float)
    DIP_S_Obs = LO_DIP_S_data['OP'].astype(float)
    DIP_Obs = (DIP_N_Obs+DIP_S_Obs)/2
    
    Nit_Denit_Opt = 'Opt1'
    
    #NO Atmospheric Deposition
    Atm_Deposition = (714.6*1000*1000) #mg/day
    Atm_Deposition_N = N_Per * Atm_Deposition
    Atm_Deposition_S = S_Per * Atm_Deposition
    
    #Calibration parameters 
    Cal_Par_Opt = pd.read_csv('./Model_Data_Filled_20082023/Cal_Par.csv')
    Cal_Par = Cal_Par_Opt['Par']
    
       
    #nitrification
    #Either a simple method where Nit_R_T = nitrification rate in water column (Tot_nit_R) * Temp Coeff for nitrification * (T-20)
    #Or a more detailed method R = k0 + k * fam * fox
    # two options for fam and fox 
    #Option1
    # fam = (Cam/Ks+Cam) and fox = (Cox/Ksox+Cox)
    #Option 2
    #Fam = Cam and fox = (1-foxmin) * (Cox-Coxc/Coxo-Coxc) ^10a + foxmin
    
    #The nitrification rate for the simple method
    Tot_nit_R = 0.3
    # nitrification parameters for the detailed method
    kni0 = Cal_Par[0]
    kni = Cal_Par[1]
    Theta_ni = 1.06
    
    #Denitrification
    #Either a simple method where Denit_R_T = denitrification rate in water column (Tot_denit_R) * Temp Coeff for denitrification * (T-20)
    #Or a more detailed method R = k0 + k * fni * fox
    # two options for fam and fox 
    #Option1
    # fni = (Cni/Ks+Cni) and fox = 1- (Cox/Ksox+Cox)
    #Option 2
    #Fni = Cni and fox = (Coxc-Cox/Coxc-Coxo) 
   
    #The denitrification rate for the simple method
    Tot_denit_R = 0.1
    #The denitrification rate for the simple method
    kden0 =  Cal_Par[2]
    kden =  Cal_Par[3]
    Theta_deni = 1.06
    
    #Light Attenuation due to water Kw and algae Kc
    Kw = Cal_Par[4] #1.7
    Kc = Cal_Par[5] #0.029
    #Light limitation and inhibition coefficients K1 and K2
    K1 = Cal_Par[6] #100
    K2 = Cal_Par[7] #700
    # Maximum Phytoplankton growth rate
    G_max = Cal_Par[8] #1.5
    
    #Half Saturation Coefficients
    K_NH = Cal_Par[9] #0.1
    K_Nitr = Cal_Par[10]
    K_DIN = K_NH + K_Nitr
    KP = Cal_Par[11]
    K_TN = Cal_Par[12] #0.1
    YNOChla = Cal_Par[13] #0.1
    
    #Temperatures
    T_opt = Cal_Par[14] #15 #C
    T_min = Cal_Par[15] #10 #C
    T_max = Cal_Par[16] #30 #C
    # Dissolved Oxygen    
    KDO = Cal_Par[17]
    #Sediment Release
    S_NO = Cal_Par[18] #mg/m2/d
    
    Theta_NO = 1.06
    NO_Temp_Adj = Theta_NO**(Temp-20)  
    Area = 1730 *1E6 #m2
    Area_N = Area * N_Per
    Area_S = Area * S_Per
    
    X = len(Q_in.index)
    
    #Nitrogen and Phosphorus Limiting
    f_P = DIP_Obs/(KP+DIP_Obs)
    f_N = DIN_Obs/(K_DIN+DIN_Obs)
    
    f_P_N = DIP_N_Obs/(KP+DIP_N_Obs)
    f_P_S = DIP_S_Obs/(KP+DIP_S_Obs)
    
    f_N_N = DIN_N_Obs/(K_DIN+DIN_N_Obs)
    f_N_S = DIN_S_Obs/(K_DIN+DIN_S_Obs)
    
    ### In case f_N is calculated
    # f_N_N = np.zeros(X,dtype = object)
    # f_N_S = np.zeros(X,dtype = object)
    # f_N = np.zeros(X,dtype = object)

    Q_I_M = np.zeros(X,dtype = object)
    Q_O_M = np.zeros(X,dtype = object)
    Q_N2S = np.zeros(X,dtype = object)
    External_NO_M = np.zeros(X,dtype = object)

    Nit_R_N = np.zeros(X,dtype = object)
    Nit_R_N_T = np.zeros(X,dtype = object)
    Nit_R_S = np.zeros(X,dtype = object)
    Nit_R_S_T = np.zeros(X,dtype = object)
    Denit_R_N = np.zeros(X,dtype = object)
    Denit_R_N_T = np.zeros(X,dtype = object)
    Denit_R_S = np.zeros(X,dtype = object)
    Denit_R_S_T = np.zeros(X,dtype = object)
    
    fT = np.zeros(X,dtype = object)
    fL = np.zeros(X,dtype = object)
    NO_N = np.zeros(X,dtype = object)
    NO_S = np.zeros(X,dtype = object)
    NO_MEAN = np.zeros(X,dtype = object)
            
    #Water depth could be stage - B.L. (variable) or a constant average depth (e.g., 2.5 m)

    z = np.zeros(X,dtype = object)
    # B_L = 1 #m
    # z = Stage - B_L #m

    for i in range(X):
        z[i] = 2.5
    
    fox_min = Cal_Par[19]
    DO_Cr_n = Cal_Par[20]
    DO_Opt_n = Cal_Par[21]
    a = Cal_Par[22]
    DO_Cr_d = Cal_Par[23]
    DO_Opt_d = Cal_Par[24]
    
    Grazing = 0
    Theta_G = 1.06
    Grazing_T = Grazing*Theta_G**(Temp-20)
    
    ### Initial Values
    NO_N[0] = NOx_N_Obs[0]
    NO_S[0] = NOx_S_Obs[0]
    
    NO_MEAN[0] = (NO_N[0]+NO_S[0])*0.5
    
    # Map the monthly values to the DataFrame
    Nitro_Model_Output = pd.DataFrame(Q_in['date'],columns=['date'])
    Nitro_Model_Output['date'] = pd.to_datetime(Nitro_Model_Output['date'])
    print("LOONE Nitrogen Module is Running!")
    for i in range(X-1):
        print(Nitro_Model_Output['date'].iloc[i])
        
        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48 #m3/d
            Q_O_M[i] = Q_O[i]
            External_NO_M[i] = External_NO[i] + Q_I_M[i] * S65E_NO[i]            
    
        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48 #m3/d
            Q_I_M[i] = Q_I[i]
            External_NO_M[i] = External_NO[i]         
    
        Q_N2S[i] = (Q_I_M[i]*1 + Q_O_M[i]*0)   
    
        t = np.linspace(1,2,num = 2)
        Nit_R_N[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0,kni,NH4_N_Obs[i],K_NH,DO[i],KDO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_N_T[i] = Nit_R_N[i]*Theta_ni**(Temp[i]-20) 
        Nit_R_S[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0,kni,NH4_S_Obs[i],K_NH,DO[i],KDO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_S_T[i] = Nit_R_S[i]*Theta_ni**(Temp[i]-20) 
    
        Denit_R_N[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0,kden,NO_N[i],K_Nitr,DO[i],KDO,DO_Cr_d,DO_Opt_d)
        Denit_R_N_T[i] = Denit_R_N[i]*Theta_deni**(Temp[i]-20) 
        Denit_R_S[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0,kden,NO_S[i],K_Nitr,DO[i],KDO,DO_Cr_d,DO_Opt_d)
        Denit_R_S_T[i] = Denit_R_S[i]*Theta_deni**(Temp[i]-20) 
        
        # f_N_N[i] = (NO_N[i]+ NH4_N_Obs[i])/(NO_N[i]+ NH4_N_Obs[i]+K_DIN)
        # f_N_S[i] = (NO_S[i]+ NH4_S_Obs[i])/(NO_S[i]+ NH4_S_Obs[i]+K_DIN)
        # f_N[i] = (NO_MEAN[i]+ NH4_Obs[i])/(NO_MEAN[i]+ NH4_Obs[i]+K_DIN)
        
        # f_N_N[i] = (NO_N[i]+ NH4_N[i])/(NO_N[i]+ NH4_N[i]+K_DIN)
        # f_N_S[i] = (NO_S[i]+ NH4_S[i])/(NO_S[i]+ NH4_S[i]+K_DIN)
        # f_N[i] = (NO_MEAN[i]+ NH4_Obs[i])/(NO_MEAN[i]+ NH4_Obs[i]+K_DIN)
    
    
        fT[i] = f_T_NO_alt1(Temp[i],T_opt,T_min,T_max)
        # fL[i] = f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw,Kc,Sim_Chla[i],z[i],K1,K2)*1.5 if Nitro_Model_Output['date'].iloc[i].month in (6,7,8,9,10) else f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw,Kc,Sim_Chla[i],z[i],K1,K2)*0.75
        fL[i] = f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw,Kc,Chla[i],z[i],K1,K2)

        DfEq_Res_N = odeint(NOx_N_DiffEq,NO_N[i],t,args=(External_NO[i],Atm_Deposition_N,Q_N2S[i],Nit_R_N_T[i],Denit_R_N_T[i],NH4_N_Obs[i],volume_N[i],G_max,fT[i],fL[i],f_P_N[i],f_N_N[i],K_NH,K_TN,YNOChla,Chla_N[i],S_NO,Area_N,NO_Temp_Adj[i],))
        NO_N[i+1] = DfEq_Res_N[:,0][1]
        DfEq_Res_S = odeint(NOx_S_DiffEq,NO_S[i],t,args=(Atm_Deposition_S,Q_N2S[i],Q_O_M[i],NO_N[i],Nit_R_S_T[i],Denit_R_S_T[i],NH4_S_Obs[i],volume_S[i],G_max,fT[i],fL[i],f_P_S[i],f_N_S[i],K_NH,K_TN,YNOChla,Chla_S[i],S_NO,Area_S,NO_Temp_Adj[i],))
        NO_S[i+1] = DfEq_Res_S[:,0][1]
        NO_MEAN[i+1] = (NO_N[i+1] + NO_S[i+1])*0.5
   
    print("Exporting Module Outputs!")

    Nitro_Model_Output['NO_N'] = pd.to_numeric(NO_N)
    Nitro_Model_Output['NO_S'] = pd.to_numeric(NO_S)
    Nitro_Model_Output['NO_M'] = pd.to_numeric(NO_MEAN)    
    
    Nitro_Model_Output = Nitro_Model_Output.set_index(['date'])
    Nitro_Model_Output.index = pd.to_datetime(Nitro_Model_Output.index, unit = 'ns')
    Nitro_Mod_Out_M =  Nitro_Model_Output.resample('M').mean()
    Nitro_Model_Output = Nitro_Model_Output.reset_index()
    Nitro_Mod_Out_M = Nitro_Mod_Out_M.reset_index()

    return(Nitro_Mod_Out_M)

#### NO model given simulated Outflows from LOONE Q Module
def LOONE_NO_SimQ(LOONE_Q_Outputs):
    #Read Required Data
    Q_in = pd.read_csv('./Model_Data_Filled_20082023/LO_Inflows_20082023.csv')
    Temp_data = pd.read_csv('./Model_Data_Filled_20082023/Temp_Avg_20082023.csv')
    DO_data = pd.read_csv('./Model_Data_Filled_20082023/LO_Avg_DO_20082023.csv')
    RAD_data = pd.read_csv('./Model_Data_Filled_20082023/LO_RADT_20082023.csv')
    Storage = pd.read_csv('./Model_Data_Filled_20082023/Average_LO_Storage_20082023.csv')
    Storage_deviation = pd.read_csv('./Model_Data_Filled_20082023/Storage_Dev_20082023.csv')
    Chla_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_Merged_Chla.csv') #microgram/L
    Chla_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_Merged_Chla.csv') #microgram/L
    NOx_In = pd.read_csv('./Model_Data_Filled_20082023/Daily_NOx_External_Loads_20082023.csv') #mg
    S65E_NO_data = pd.read_csv('./Model_Data_Filled_20082023/S65E_NO_Interpolated.csv') #mg/m3
    Photoperiod = pd.read_csv('./Model_Data_Filled_20082023/PhotoPeriod_20082023.csv') 
    Photoperiod['Date'] = pd.to_datetime(Photoperiod['Date'])
    LO_DIP_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_OP.csv') #mg/m3
    LO_DIP_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_OP.csv') #mg/m3
    LO_DIN_N_data = pd.read_csv('./Model_Data_Filled_20082023/N_DIN.csv') #mg/m3
    LO_DIN_S_data = pd.read_csv('./Model_Data_Filled_20082023/S_DIN.csv') #mg/m3
    RAD = RAD_data['RAD'].astype(float) * 4.6*1000 

    Temp = Temp_data['Mean_T'].astype(float)
    DO = DO_data['DO'].astype(float)
    External_NO = NOx_In['NOx_Load_mg'].astype(float) #mg
    S65E_NO = S65E_NO_data['NO'].astype(float) #mg/m3
    # N-S Procedure
    N_Per = 0.43
    S_Per = 0.57
    
    volume  = Storage['Storage_m3'].astype(float) #m3
    volume_N = volume*N_Per
    volume_S = volume*S_Per
    Storage_dev = Storage_deviation['DS_dev'].astype(float) #acft
    Stage = Storage['Stage_ft'].astype(float) * 0.3048 #m
    Q_I = Q_in['Inflows_cmd'].astype(float) #m3
    Q_O = (LOONE_Q_Outputs['S77EW'] *0.028316847 + ((LOONE_Q_Outputs['TotRegEW'] + LOONE_Q_Outputs['TotRegSo'])/70.0456)) * 3600 * 24
    S77_Q = LOONE_Q_Outputs['S77_Q']
    S308_Q = LOONE_Q_Outputs['S308_Q']
    TotRegSo = LOONE_Q_Outputs['TotRegSo'] #acft/day

    NH4_N_Obs = LO_DIN_N_data['NH4'].astype(float) #mg/m3
    NH4_S_Obs = LO_DIN_S_data['NH4'].astype(float) #mg/m3
    NH4_Obs = (NH4_N_Obs+NH4_S_Obs)/2
    NOx_N_Obs = LO_DIN_N_data['NO'].astype(float) #mg/m3
    NOx_S_Obs = LO_DIN_S_data['NO'].astype(float) #mg/m3
    NOx_Obs = (NOx_N_Obs+NOx_S_Obs)/2
    Chla_N = Chla_N_data['Chla'].astype(float) #microgram/L = mg/m3
    Chla_S = Chla_S_data['Chla'].astype(float) #microgram/L = mg/m3
    Chla = (Chla_N+Chla_S)/2 
    DIN_N_Obs = LO_DIN_N_data['DIN'].astype(float)
    DIN_S_Obs = LO_DIN_S_data['DIN'].astype(float)
    DIN_Obs = (DIN_N_Obs+DIN_S_Obs)/2
    DIP_N_Obs = LO_DIP_N_data['OP'].astype(float)
    DIP_S_Obs = LO_DIP_S_data['OP'].astype(float)
    DIP_Obs = (DIP_N_Obs+DIP_S_Obs)/2
    
    Nit_Denit_Opt = 'Opt1'
    
    #NO Atmospheric Deposition
    Atm_Deposition = (714.6*1000*1000) #mg/day
    Atm_Deposition_N = N_Per * Atm_Deposition
    Atm_Deposition_S = S_Per * Atm_Deposition
    
    #Calibration parameters 
    Cal_Par_Opt = pd.read_csv('./Model_Data_Filled_20082023/Cal_Par.csv')
    Cal_Par = Cal_Par_Opt['Par']
    
       
    #nitrification
    #The nitrification rate for the simple method
    Tot_nit_R = 0.3
    # nitrification parameters for the detailed method
    kni0 = Cal_Par[0]
    kni = Cal_Par[1]
    Theta_ni = 1.06
    
    #Denitrification
    #The denitrification rate for the simple method
    Tot_denit_R = 0.1
    #The denitrification rate for the simple method
    kden0 =  Cal_Par[2]
    kden =  Cal_Par[3]
    Theta_deni = 1.06
    
    #Light Attenuation due to water Kw and algae Kc
    Kw = Cal_Par[4] #1.7
    Kc = Cal_Par[5] #0.029
    #Light limitation and inhibition coefficients K1 and K2
    K1 = Cal_Par[6] #100
    K2 = Cal_Par[7] #700
    # Maximum Phytoplankton growth rate
    G_max = Cal_Par[8] #1.5
    
    #Half Saturation Coefficients
    K_NH = Cal_Par[9] #0.1
    K_Nitr = Cal_Par[10]
    K_DIN = K_NH + K_Nitr
    KP = Cal_Par[11]
    K_TN = Cal_Par[12] #0.1
    YNOChla = Cal_Par[13] #0.1
    
    #Temperatures
    T_opt = Cal_Par[14] #15 #C
    T_min = Cal_Par[15] #10 #C
    T_max = Cal_Par[16] #30 #C
    # Dissolved Oxygen    
    KDO = Cal_Par[17]
    #Sediment Release
    S_NO = Cal_Par[18] #mg/m2/d
    
    Theta_NO = 1.06
    NO_Temp_Adj = Theta_NO**(Temp-20)  
    Area = 1730 *1E6 #m2
    Area_N = Area * N_Per
    Area_S = Area * S_Per
    
    X = len(Q_in.index)
    
    #Nitrogen and Phosphorus Limiting
    f_P = DIP_Obs/(KP+DIP_Obs)
    f_N = DIN_Obs/(K_DIN+DIN_Obs)
    
    f_P_N = DIP_N_Obs/(KP+DIP_N_Obs)
    f_P_S = DIP_S_Obs/(KP+DIP_S_Obs)
    
    f_N_N = DIN_N_Obs/(K_DIN+DIN_N_Obs)
    f_N_S = DIN_S_Obs/(K_DIN+DIN_S_Obs)
    
    Q_I_M = np.zeros(X,dtype = object)
    Q_O_M = np.zeros(X,dtype = object)
    Q_N2S = np.zeros(X,dtype = object)
    External_NO_M = np.zeros(X,dtype = object)

    Nit_R_N = np.zeros(X,dtype = object)
    Nit_R_N_T = np.zeros(X,dtype = object)
    Nit_R_S = np.zeros(X,dtype = object)
    Nit_R_S_T = np.zeros(X,dtype = object)
    Denit_R_N = np.zeros(X,dtype = object)
    Denit_R_N_T = np.zeros(X,dtype = object)
    Denit_R_S = np.zeros(X,dtype = object)
    Denit_R_S_T = np.zeros(X,dtype = object)
    
    fT = np.zeros(X,dtype = object)
    fL = np.zeros(X,dtype = object)
    NO_N = np.zeros(X,dtype = object)
    NO_S = np.zeros(X,dtype = object)
    NO_MEAN = np.zeros(X,dtype = object)
    
    NO_Load_Cal = np.zeros(X,dtype = object)
    NO_Load_StL = np.zeros(X,dtype = object)
    NO_Load_South = np.zeros(X,dtype = object)

    
    z = np.zeros(X,dtype = object)
    for i in range(X):
        z[i] = 2.5
    
    fox_min = Cal_Par[19]
    DO_Cr_n = Cal_Par[20]
    DO_Opt_n = Cal_Par[21]
    a = Cal_Par[22]
    DO_Cr_d = Cal_Par[23]
    DO_Opt_d = Cal_Par[24]
    
    
    Grazing = 0
    Theta_G = 1.06
    Grazing_T = Grazing*Theta_G**(Temp-20)
    NO_N[0] = NOx_N_Obs[0]
    NO_S[0] = NOx_S_Obs[0]
    
    NO_MEAN[0] = (NO_N[0]+NO_S[0])*0.5
    
    # Map the monthly values to the DataFrame
    Nitro_Model_Output = pd.DataFrame(Q_in['date'],columns=['date'])
    Nitro_Model_Output['date'] = pd.to_datetime(Nitro_Model_Output['date'])
    print("LOONE Nitrogen Module is Running!")
    for i in range(X-1):
        print(Nitro_Model_Output['date'].iloc[i])
        
        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48 #m3/d
            Q_O_M[i] = Q_O[i]
            External_NO_M[i] = External_NO[i] + Q_I_M[i] * S65E_NO[i]            
    
        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48 #m3/d
            Q_I_M[i] = Q_I[i]
            External_NO_M[i] = External_NO[i]         
    
        Q_N2S[i] = (Q_I_M[i]*1 + Q_O_M[i]*0)   
    
        t = np.linspace(1,2,num = 2)
        Nit_R_N[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0,kni,NH4_N_Obs[i],K_NH,DO[i],KDO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_N_T[i] = Nit_R_N[i]*Theta_ni**(Temp[i]-20) 
        Nit_R_S[i] = Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0,kni,NH4_S_Obs[i],K_NH,DO[i],KDO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_S_T[i] = Nit_R_S[i]*Theta_ni**(Temp[i]-20) 
    
        Denit_R_N[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0,kden,NO_N[i],K_Nitr,DO[i],KDO,DO_Cr_d,DO_Opt_d)
        Denit_R_N_T[i] = Denit_R_N[i]*Theta_deni**(Temp[i]-20) 
        Denit_R_S[i] = Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0,kden,NO_S[i],K_Nitr,DO[i],KDO,DO_Cr_d,DO_Opt_d)
        Denit_R_S_T[i] = Denit_R_S[i]*Theta_deni**(Temp[i]-20) 
    
        fT[i] = f_T_NO_alt1(Temp[i],T_opt,T_min,T_max)
        fL[i] = f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw,Kc,Chla[i],z[i],K1,K2)

        DfEq_Res_N = odeint(NOx_N_DiffEq,NO_N[i],t,args=(External_NO[i],Atm_Deposition_N,Q_N2S[i],Nit_R_N_T[i],Denit_R_N_T[i],NH4_N_Obs[i],volume_N[i],G_max,fT[i],fL[i],f_P_N[i],f_N_N[i],K_NH,K_TN,YNOChla,Chla_N[i],S_NO,Area_N,NO_Temp_Adj[i],))
        NO_N[i+1] = DfEq_Res_N[:,0][1]
        DfEq_Res_S = odeint(NOx_S_DiffEq,NO_S[i],t,args=(Atm_Deposition_S,Q_N2S[i],Q_O_M[i],NO_N[i],Nit_R_S_T[i],Denit_R_S_T[i],NH4_S_Obs[i],volume_S[i],G_max,fT[i],fL[i],f_P_S[i],f_N_S[i],K_NH,K_TN,YNOChla,Chla_S[i],S_NO,Area_S,NO_Temp_Adj[i],))
        NO_S[i+1] = DfEq_Res_S[:,0][1]
        NO_MEAN[i+1] = (NO_N[i+1] + NO_S[i+1])*0.5
          
        
        NO_Load_Cal[i] = S77_Q[i]*0.028316847*3600*24*NO_S[i] #mg/d P
        NO_Load_StL[i] = S308_Q[i]*0.028316847*3600*24*NO_S[i] #mg/d P
        NO_Load_South[i] = TotRegSo[i]*1233.48 *NO_S[i] #mg/d P

    print("Exporting Module Outputs!")

    Nitro_Model_Output['NO_N'] = pd.to_numeric(NO_N)
    Nitro_Model_Output['NO_S'] = pd.to_numeric(NO_S)
    Nitro_Model_Output['NO_M'] = pd.to_numeric(NO_MEAN)    
    
    Nitro_Model_Output = Nitro_Model_Output.set_index(['date'])
    Nitro_Model_Output.index = pd.to_datetime(Nitro_Model_Output.index, unit = 'ns')
    Nitro_Mod_Out_M =  Nitro_Model_Output.resample('M').mean()
    Nitro_Model_Output = Nitro_Model_Output.reset_index()
    Nitro_Mod_Out_M = Nitro_Mod_Out_M.reset_index()
    return(Nitro_Mod_Out_M)