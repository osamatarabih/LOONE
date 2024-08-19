# Daily Nitrate-Nitrite NOx and Chlorophyll-a Modeling

# import required packages
import os
import argparse
import pandas as pd
import numpy as np
from scipy.integrate import odeint
from loone.utils import loone_nchla_fns, load_config
from loone.data import Data as DClass


def LOONE_Constituent_SimQ(workspace: str):
    """Daily Nitrate-Nitrite NOx and Chlorophyll-a Modeling
    
    Args:
        workspace (str): The working directory path.
    
    Returns:
        list: The results of the simulation.
    """
    # Read in the config file
    config = load_config(workspace)
    
    # Initialize the Data object
    data = DClass(workspace)
    
    # Read Required Data
    Q_in = pd.read_csv(os.path.join(workspace, 'LO_Inflows_BK.csv'))
    Temp_data = pd.read_csv(os.path.join(workspace, 'Filled_WaterT.csv'))
    DO_data = pd.read_csv(os.path.join(workspace, 'LO_DO_Clean_daily.csv'))
    RAD_data = pd.read_csv(os.path.join(workspace, 'LO_RADT_data.csv'))
    Storage = pd.read_csv(os.path.join(workspace, 'Average_LO_Storage_3MLag.csv'))
    Chla_N_data = pd.read_csv(os.path.join(workspace, 'N_Merged_Chla.csv')) # microgram/L
    Chla_S_data = pd.read_csv(os.path.join(workspace, 'S_Merged_Chla.csv')) # microgram/L
    NOx_In = pd.read_csv(os.path.join(workspace, 'LO_External_Loadings_NO.csv')) # mg
    S65E_NO_data = pd.read_csv(os.path.join(workspace, 'water_quality_S65E_NITRATE+NITRITE-N_Interpolated.csv')) # mg/m3
    Chla_In = pd.read_csv(os.path.join(workspace, 'Chla_Loads_In.csv')) # mg
    S65E_Chla_data = pd.read_csv(os.path.join(workspace, 'S65E_Chla_Merged.csv')) # mg/m3

    Photoperiod = pd.read_csv(os.path.join(workspace, 'PhotoPeriod.csv'))
    Photoperiod['Date'] = pd.to_datetime(Photoperiod['Date'])
    LO_DIP_N_data = pd.read_csv(os.path.join(workspace, 'N_OP.csv')) # mg/m3
    LO_DIP_S_data = pd.read_csv(os.path.join(workspace, 'S_OP.csv')) # mg/m3
    LO_DIN_N_data = pd.read_csv(os.path.join(workspace, 'N_DIN.csv')) # mg/m3
    LO_DIN_S_data = pd.read_csv(os.path.join(workspace, 'S_DIN.csv')) # mg/m3

    date_start = Q_in['date'].iloc[0]
    
    RAD = RAD_data['Mean_RADT'].astype(float) * 4.6*1000 

    Temp = Temp_data['Water_T'].astype(float)
    DO = DO_data['Mean_DO'].astype(float)
    External_NO = NOx_In['External_NO_Ld_mg'].astype(float) # mg
    S65E_NO = (S65E_NO_data[S65E_NO_data['date'] >= date_start]['Data'] * 1000).astype(float) # mg/m3
    
    External_Chla = Chla_In['Chla_Loads'].astype(float)*3 # mg
    S65E_Chla = S65E_Chla_data['Data'].astype(float)

    # N-S Procedure
    N_Per = 0.43
    S_Per = 0.57
    
    Storage_dev = data.Stroage_dev_df['DS_dev'].astype(float) # acft
    Q_I = Q_in['Inflows_cmd'].astype(float) # m3

   
    # Simulated Q
    # S77_Q = LOONE_Q_Outputs['S77_Q'] * 0.0283168 * 86400 # cfs to cubic meters per day
    # S308_Q = LOONE_Q_Outputs['S308_Q'] * 0.0283168 * 86400 # cfs to cubic meters per day
    # TotRegSo = LOONE_Q_Outputs['TotRegSo'] # acft/day
    # TotRegEW = LOONE_Q_Outputs['TotRegEW'] # acft/day

    # Simulated Stage and Storage
    # Stage = LOONE_Q_Outputs['Stage'] * 0.3048 # ft to m
    # volume = LOONE_Q_Outputs['Storage'] * 1233.48 # acft to m3
    
    # Observed S77 S308 South
    Outflows_Obs = pd.read_csv(os.path.join(workspace, f'Flow_df_3MLag.csv'))
    S77_Q = Outflows_Obs['S77_Out']
    S308_Q = Outflows_Obs['S308_Out']
    TotRegSo = Outflows_Obs[['S351_Out','S354_Out','S352_Out','L8_Out']].sum(axis=1)/1233.48    # m3/day to acft
    
    # Observed Stage and Storage 
    Stage = Storage['Stage_ft'].astype(float) * 0.3048 # m
    volume  = Storage['Storage_cmd'].astype(float) # m3
    
    volume_N = volume*N_Per
    volume_S = volume*S_Per

    Q_O = S77_Q + S308_Q + TotRegSo * 1233.48  # cmd

    NH4_N_Obs = LO_DIN_N_data['NH4'].astype(float) # mg/m3
    NH4_S_Obs = LO_DIN_S_data['NH4'].astype(float) # mg/m3
    NH4_Obs = (NH4_N_Obs+NH4_S_Obs)/2
    NOx_N_Obs = LO_DIN_N_data['NO'].astype(float) # mg/m3
    NOx_S_Obs = LO_DIN_S_data['NO'].astype(float) # mg/m3
    NOx_Obs = (NOx_N_Obs+NOx_S_Obs)/2
    Chla_N = Chla_N_data['Chla'].astype(float) # microgram/L = mg/m3
    Chla_S = Chla_S_data['Chla'].astype(float) # microgram/L = mg/m3
    Chla = (Chla_N+Chla_S)/2 
    DIN_N_Obs = LO_DIN_N_data['DIN'].astype(float)
    DIN_S_Obs = LO_DIN_S_data['DIN'].astype(float)
    DIN_Obs = (DIN_N_Obs+DIN_S_Obs)/2
    DIP_N_Obs = LO_DIP_N_data['OP'].astype(float)
    DIP_S_Obs = LO_DIP_S_data['OP'].astype(float)
    DIP_Obs = (DIP_N_Obs+DIP_S_Obs)/2
    
    Nit_Denit_Opt = 'Opt1'
    
    # NO Atmospheric Deposition
    Atm_Deposition = (714.6*1000*1000) # mg/day
    Atm_Deposition_N = N_Per * Atm_Deposition
    Atm_Deposition_S = S_Per * Atm_Deposition
    
    # Calibration parameters 
    Cal_Par_NO = data.Cal_Par['Par_NO']
    Cal_Par_Chla = data.Cal_Par['Par_Chla']
 
       
    # nitrification
    # Either a simple method where Nit_R_T = nitrification rate in water column (Tot_nit_R) * Temp Coeff for nitrification * (T-20)
    # Or a more detailed method R = k0 + k * fam * fox
    # two options for fam and fox 
    # Option1
    # fam = (Cam/Ks+Cam) and fox = (Cox/Ksox+Cox)
    # Option 2
    # Fam = Cam and fox = (1-foxmin) * (Cox-Coxc/Coxo-Coxc) ^10a + foxmin
    # The nitrification rate for the simple method
    Tot_nit_R = 0.3
    # nitrification parameters for the detailed method
    kni0_NO = Cal_Par_NO[0]
    kni0_Chla = Cal_Par_Chla[0]
    kni_NO = Cal_Par_NO[1]
    kni_Chla = Cal_Par_Chla[1]

    Theta_ni = 1.06
    
    # Denitrification
    # Either a simple method where Denit_R_T = denitrification rate in water column (Tot_denit_R) * Temp Coeff for denitrification * (T-20)
    # Or a more detailed method R = k0 + k * fni * fox
    # two options for fam and fox 
    # Option1
    # fni = (Cni/Ks+Cni) and fox = 1- (Cox/Ksox+Cox)
    # Option 2
    # Fni = Cni and fox = (Coxc-Cox/Coxc-Coxo) 
    # The denitrification rate for the simple method
    Tot_denit_R = 0.1
    # The denitrification rate for the simple method
    kden0_NO =  Cal_Par_NO[2]
    kden0_Chla =  Cal_Par_Chla[2]

    kden_NO =  Cal_Par_NO[3]
    kden_Chla =  Cal_Par_Chla[3]

    Theta_deni = 1.06
    
    # Light Attenuation due to water Kw and algae Kc
    Kw_NO = Cal_Par_NO[4] # 1.7
    Kw_Chla = Cal_Par_Chla[4] # 1.7
    Kc_NO = Cal_Par_NO[5] # 0.029
    Kc_Chla = Cal_Par_Chla[5] # 0.029
    # Light limitation and inhibition coefficients K1 and K2
    K1_NO = Cal_Par_NO[6] # 100
    K1_Chla = Cal_Par_Chla[6] # 100
    K2_NO = Cal_Par_NO[7] # 700
    K2_Chla = Cal_Par_Chla[7] # 700
    # Maximum Phytoplankton growth rate
    G_max_NO = Cal_Par_NO[8] # 1.5
    G_max_Chla = Cal_Par_Chla[8] # 1.5

    # Half Saturation Coefficients
    K_NH_NO = Cal_Par_NO[9] # 0.1
    K_NH_Chla = Cal_Par_Chla[9] # 0.1
    K_Nitr_NO = Cal_Par_NO[10]
    K_Nitr_Chla = Cal_Par_Chla[10]
    K_DIN_NO = K_NH_NO + K_Nitr_NO
    K_DIN_Chla = K_NH_Chla + K_Nitr_Chla
    KP_NO = Cal_Par_NO[11]
    KP_Chla = Cal_Par_Chla[11]
    K_TN_NO = Cal_Par_NO[12] # 0.1
    K_TN_Chla = Cal_Par_Chla[12] # 0.1
    YNOChla_NO = Cal_Par_NO[13] # 0.1
    YNOChla_Chla = Cal_Par_Chla[13] # 0.1

    
    # Temperatures
    T_opt_NO = Cal_Par_NO[14] # 15 # C
    T_min_NO = Cal_Par_NO[15] # 10 # C
    T_max_NO = Cal_Par_NO[16] # 30 # C
    
    T_opt_Chla = Cal_Par_Chla[14] # 15 # C
    T_min_Chla = Cal_Par_Chla[15] # 10 # C
    T_max_Chla = Cal_Par_Chla[16] # 30 # C

    # Dissolved Oxygen    
    KDO_NO = Cal_Par_NO[17]
    KDO_Chla = Cal_Par_Chla[17]

    # Sediment Release
    S_NO_NO = Cal_Par_NO[18] # mg/m2/d
    S_NO_Chla = Cal_Par_Chla[18] # mg/m2/d

    Theta_NO = 1.06
    NO_Temp_Adj = Theta_NO**(Temp-20)  
    Area = 1730 *1E6 # m2
    Area_N = Area * N_Per
    Area_S = Area * S_Per
    
    X = len(Q_in.index)
    
    # Nitrogen and Phosphorus Limiting
    f_P_NO = DIP_Obs/(KP_NO+DIP_Obs)
    f_P_Chla = DIP_Obs/(KP_Chla+DIP_Obs)

    f_N_NO = DIN_Obs/(K_DIN_NO+DIN_Obs)
    f_N_Chla = DIN_Obs/(K_DIN_Chla+DIN_Obs)

    f_P_N_NO = DIP_N_Obs/(KP_NO+DIP_N_Obs)
    f_P_N_Chla = DIP_N_Obs/(KP_Chla+DIP_N_Obs)

    f_P_S_NO = DIP_S_Obs/(KP_NO+DIP_S_Obs)
    f_P_S_Chla = DIP_S_Obs/(KP_Chla+DIP_S_Obs)

    f_N_N_NO = DIN_N_Obs/(K_DIN_NO+DIN_N_Obs)
    f_N_N_Chla = DIN_N_Obs/(K_DIN_Chla+DIN_N_Obs)

    f_N_S_NO = DIN_S_Obs/(K_DIN_NO+DIN_S_Obs)
    f_N_S_Chla = DIN_S_Obs/(K_DIN_Chla+DIN_S_Obs)

    Q_I_M = np.zeros(X,dtype = object)
    Q_O_M = np.zeros(X,dtype = object)
    Q_N2S = np.zeros(X,dtype = object)
    External_NO_M = np.zeros(X,dtype = object)
    External_Chla_M = np.zeros(X,dtype = object)

    Nit_R_N_NO = np.zeros(X,dtype = object)
    Nit_R_N_T_NO = np.zeros(X,dtype = object)
    Nit_R_S_NO = np.zeros(X,dtype = object)
    Nit_R_S_T_NO = np.zeros(X,dtype = object)
    Denit_R_N_NO = np.zeros(X,dtype = object)
    Denit_R_N_T_NO = np.zeros(X,dtype = object)
    Denit_R_S_NO = np.zeros(X,dtype = object)
    Denit_R_S_T_NO = np.zeros(X,dtype = object)
    

    fT_NO = np.zeros(X,dtype = object)
    fL_NO = np.zeros(X,dtype = object)
    
    fT_Chla = np.zeros(X,dtype = object)
    fL_Chla = np.zeros(X,dtype = object)

    NO_N = np.zeros(X,dtype = object)
    NO_S = np.zeros(X,dtype = object)
    NO_MEAN = np.zeros(X,dtype = object)
    
    NO_Load_Cal = np.zeros(X,dtype = object)
    NO_Load_StL = np.zeros(X,dtype = object)
    NO_Load_South = np.zeros(X,dtype = object)

    Sim_Chla_N = np.zeros(X,dtype = object)
    Sim_Chla_S = np.zeros(X,dtype = object)
    Sim_Chla = np.zeros(X,dtype = object)
    
    Chla_Load_Cal = np.zeros(X,dtype = object)
    Chla_Load_StL = np.zeros(X,dtype = object)
    Chla_Load_South = np.zeros(X,dtype = object)

    
    # Water depth could be stage - B.L. (variable) or a constant average depth (e.g., 2.5 m)

    z = np.zeros(X,dtype = object)
    # B_L = 1 # m
    # z = Stage - B_L # m

    for i in range(X):
        z[i] = 2.5
    
    fox_min = Cal_Par_NO[19]
    DO_Cr_n = Cal_Par_NO[20]
    DO_Opt_n = Cal_Par_NO[21]
    a = Cal_Par_NO[22]
    DO_Cr_d = Cal_Par_NO[23]
    DO_Opt_d = Cal_Par_NO[24]
   
    vc = Cal_Par_Chla[28] # phytoplankton settling velocity (m/d).

    K_r = Cal_Par_Chla[25]
    Theta_r = 1.06
    K_r_T = K_r*Theta_r**(Temp-20) 

    YNHChla = Cal_Par_Chla[26] #0.05
    # Chla mortality
    k_m = Cal_Par_Chla[27]
    Theta_m = 1.06
    K_m_T = k_m*Theta_m**(Temp-20) 
    
    Im = Cal_Par_Chla[31]

    Grazing = 0
    Theta_G = 1.06
    Grazing_T = Grazing*Theta_G**(Temp-20)
    # Initial values
    NO_N[0] = NOx_N_Obs[0]
    NO_S[0] = NOx_S_Obs[0]
    
    NO_MEAN[0] = (NO_N[0]+NO_S[0])*0.5
    
    Sim_Chla_N[0] = Chla_N[0]
    Sim_Chla_S[0] = Chla_S[0]
    Sim_Chla[0] = Chla[0]

    Nitro_Model_Output = pd.DataFrame(Q_in['date'],columns=['date'])
    Nitro_Model_Output['date'] = pd.to_datetime(Nitro_Model_Output['date'])
    print("LOONE Nitrogen Module is Running!")
    for i in range(X-1):
        # print(Nitro_Model_Output['date'].iloc[i])
        
        if Storage_dev[i] >= 0:
            Q_I_M[i] = Q_I[i] + Storage_dev[i] * 1233.48 # m3/d
            Q_O_M[i] = Q_O[i]
            External_NO_M[i] = External_NO[i] + Q_I_M[i] * S65E_NO[i]            
            External_Chla_M[i] = External_Chla[i] + Q_I_M[i] * S65E_Chla[i]    

        else:
            Q_O_M[i] = Q_O[i] - Storage_dev[i] * 1233.48 # m3/d
            Q_I_M[i] = Q_I[i]
            External_NO_M[i] = External_NO[i]         
            External_Chla_M[i] = External_Chla[i]

        Q_N2S[i] = (Q_I_M[i]*1 + Q_O_M[i]*0)
    
        t = np.linspace(1,2,num = 2)
        Nit_R_N_NO[i] = loone_nchla_fns.Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0_NO,kni_NO,NH4_N_Obs[i],K_NH_NO,DO[i],KDO_NO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_N_T_NO[i] = Nit_R_N_NO[i]*Theta_ni**(Temp[i]-20) 
        Nit_R_S_NO[i] = loone_nchla_fns.Nit_Rate('%s'%Nit_Denit_Opt,Tot_nit_R,kni0_NO,kni_NO,NH4_S_Obs[i],K_NH_NO,DO[i],KDO_NO,fox_min,DO_Cr_n,DO_Opt_n,a)  
        Nit_R_S_T_NO[i] = Nit_R_S_NO[i]*Theta_ni**(Temp[i]-20) 
    
        Denit_R_N_NO[i] = loone_nchla_fns.Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0_NO,kden_NO,NO_N[i],K_Nitr_NO,DO[i],KDO_NO,DO_Cr_d,DO_Opt_d)
        Denit_R_N_T_NO[i] = Denit_R_N_NO[i]*Theta_deni**(Temp[i]-20) 
        Denit_R_S_NO[i] = loone_nchla_fns.Denit_Rate('%s'%Nit_Denit_Opt,Tot_denit_R,kden0_NO,kden_NO,NO_S[i],K_Nitr_NO,DO[i],KDO_NO,DO_Cr_d,DO_Opt_d)
        Denit_R_S_T_NO[i] = Denit_R_S_NO[i]*Theta_deni**(Temp[i]-20) 
    
    
        fT_NO[i] = loone_nchla_fns.f_T_alt1(Temp[i],T_opt_NO,T_min_NO,T_max_NO)
        fL_NO[i] = loone_nchla_fns.f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw_NO,Kc_NO,Sim_Chla[i],z[i],K1_NO,K2_NO)

        DfEq_Res_N = odeint(loone_nchla_fns.NOx_N_DiffEq,NO_N[i],t,args=(External_NO[i],Atm_Deposition_N,Q_N2S[i],Nit_R_N_T_NO[i],Denit_R_N_T_NO[i],NH4_N_Obs[i],volume_N[i],G_max_NO,fT_NO[i],fL_NO[i],f_P_N_NO[i],f_N_N_NO[i],K_NH_NO,K_TN_NO,YNOChla_NO,Sim_Chla_N[i],S_NO_NO,Area_N,NO_Temp_Adj[i],))
        NO_N[i+1] = DfEq_Res_N[:,0][1]
        DfEq_Res_S = odeint(loone_nchla_fns.NOx_S_DiffEq,NO_S[i],t,args=(Atm_Deposition_S,Q_N2S[i],Q_O_M[i],NO_N[i],Nit_R_S_T_NO[i],Denit_R_S_T_NO[i],NH4_S_Obs[i],volume_S[i],G_max_NO,fT_NO[i],fL_NO[i],f_P_S_NO[i],f_N_S_NO[i],K_NH_NO,K_TN_NO,YNOChla_NO,Sim_Chla_S[i],S_NO_NO,Area_S,NO_Temp_Adj[i],))
        NO_S[i+1] = DfEq_Res_S[:,0][1]
        NO_MEAN[i+1] = (NO_N[i+1] + NO_S[i+1])*0.5
          
        
        NO_Load_Cal[i] = S77_Q[i]*NO_S[i] # mg/d P
        NO_Load_StL[i] = S308_Q[i]*NO_S[i] # mg/d P
        NO_Load_South[i] = TotRegSo[i]*1233.48 *NO_S[i] # mg/d P
        
        
        ##### Chla 
        fT_Chla[i] = loone_nchla_fns.f_T__Chla_alt1(Temp[i],T_opt_Chla,T_min_Chla,T_max_Chla,Nitro_Model_Output['date'].iloc[i].month)
        
        fL_Chla[i] = loone_nchla_fns.f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw_Chla,Kc_Chla,Sim_Chla[i],z[i],K1_Chla,K2_Chla)*1.2 if Nitro_Model_Output['date'].iloc[i].month in (6,7,8,9,10) else loone_nchla_fns.f_L_alt1(Photoperiod['Data'].iloc[i],RAD[i],Kw_Chla,Kc_Chla,Sim_Chla[i],z[i],K1_Chla,K2_Chla)*1

        Sim_Chla_N[i+1] = loone_nchla_fns.Chla_N_alt1(External_Chla_M[i], Q_N2S[i], Sim_Chla_N[i], vc, z[i], K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla[i], fL_Chla[i], f_P_N_Chla[i], f_N_N_Chla[i], volume_N[i])
        Sim_Chla_S[i+1] = loone_nchla_fns.Chla_S_alt1(Q_N2S[i], Q_O_M[i], Sim_Chla_N[i], Sim_Chla_S[i], vc, z[i], K_m_T[i], K_r_T[i], 0, G_max_Chla, fT_Chla[i], fL_Chla[i], f_P_S_Chla[i], f_N_S_Chla[i], volume_S[i])
        Sim_Chla[i+1] = (Sim_Chla_N[i+1] + Sim_Chla_S[i+1])/2

        Chla_Load_Cal[i] = S77_Q[i]*Sim_Chla_S[i] # mg/d P
        Chla_Load_StL[i] = S308_Q[i]*Sim_Chla_S[i] # mg/d P
        Chla_Load_South[i] = TotRegSo[i]*1233.48 *Sim_Chla_S[i] # mg/d P

    print("Exporting Module Outputs!")

    Nitro_Model_Output['NO_N'] = pd.to_numeric(NO_N)
    Nitro_Model_Output['NO_S'] = pd.to_numeric(NO_S)
    Nitro_Model_Output['NO_M'] = pd.to_numeric(NO_MEAN)    
    
    Nitro_Model_Output['Sim_Chla'] = pd.to_numeric(Sim_Chla)
    Nitro_Model_Output['Sim_Chla_N'] = pd.to_numeric(Sim_Chla_N)
    Nitro_Model_Output['Sim_Chla_S'] = pd.to_numeric(Sim_Chla_S)
    Nitro_Model_Output['Water Temp'] = pd.to_numeric(Temp) # C

    Nitro_Model_Output = Nitro_Model_Output.set_index(['date'])
    Nitro_Model_Output.index = pd.to_datetime(Nitro_Model_Output.index, unit = 'ns')
    Nitro_Mod_Out_M =  Nitro_Model_Output.resample('ME').mean()
    Nitro_Model_Output = Nitro_Model_Output.reset_index()
    Nitro_Mod_Out_M = Nitro_Mod_Out_M.reset_index()

    Constit_Loads_df = pd.DataFrame(Q_in['date'], columns=['date'])
    Constit_Loads_df['date'] = pd.to_datetime(Constit_Loads_df['date'])

    Constit_Loads_df['NO_Load_Cal'] = pd.to_numeric(NO_Load_Cal)/1E9 # tons
    Constit_Loads_df['NO_Load_StL'] = pd.to_numeric(NO_Load_StL)/1E9 # tons
    Constit_Loads_df['NO_Load_South'] = pd.to_numeric(NO_Load_South)/1E9 # tons
    Constit_Loads_df['Chla_Load_Cal'] = pd.to_numeric(Chla_Load_Cal)/1E6 # Kgs
    Constit_Loads_df['Chla_Load_StL'] = pd.to_numeric(Chla_Load_StL)/1E6 # Kgs
    Constit_Loads_df['Chla_Load_South'] = pd.to_numeric(Chla_Load_South)/1E6 # Kgs

    
    Constit_Loads_df = Constit_Loads_df.set_index('date')
    Constit_Loads_df.index = pd.to_datetime(Constit_Loads_df.index, unit = 'ns')
    Constit_Loads_M = Constit_Loads_df.resample('ME').sum()
    Constit_Loads_M = Constit_Loads_M.reset_index()
    
    
    Smr_Mnth_NOx_StL = []
    Smr_Mnth_NOx_Cal = []
    Smr_Mnth_Chla_StL = []
    Smr_Mnth_Chla_Cal = []
    for i in range(len(Constit_Loads_M.index)):
        if Constit_Loads_M['date'].iloc[i].month in [5,6,7,8,9,10]:
            Smr_Mnth_NOx_StL.append(Constit_Loads_M['NO_Load_StL'].iloc[i])
            Smr_Mnth_Chla_StL.append(Constit_Loads_M['Chla_Load_StL'].iloc[i])

            Smr_Mnth_NOx_Cal.append(Constit_Loads_M['NO_Load_Cal'].iloc[i])
            Smr_Mnth_Chla_Cal.append(Constit_Loads_M['Chla_Load_Cal'].iloc[i])

    Smr_Mnth_NOx_StL_arr = np.asarray(Smr_Mnth_NOx_StL)
    Smr_Mnth_Chla_StL_arr = np.asarray(Smr_Mnth_Chla_StL)

    Smr_Mnth_NOx_Cal_arr = np.asarray(Smr_Mnth_NOx_Cal)
    Smr_Mnth_Chla_Cal_arr = np.asarray(Smr_Mnth_Chla_Cal)

    if config['sim_type'] in [0, 1]:
        return [Constit_Loads_M,Nitro_Mod_Out_M,Smr_Mnth_NOx_StL_arr,Smr_Mnth_Chla_StL_arr,Smr_Mnth_NOx_Cal_arr,Smr_Mnth_Chla_Cal_arr]
    else:
        return [Smr_Mnth_NOx_StL_arr,Smr_Mnth_Chla_StL_arr,Smr_Mnth_NOx_Cal_arr,Smr_Mnth_Chla_Cal_arr,Nitro_Model_Output]
    
    

    # Algae_Opt_Mnth_NOx_StL = []
    # Algae_Opt_Mnth_NOx_Cal = []
    # Algae_Opt_Mnth_Chla_StL = []
    # Algae_Opt_Mnth_Chla_Cal = []

    # for i in range(len(Constit_Loads_M.index)):
    #     if Nitro_Model_Output['Water Temp'].iloc[i] >= 25: 
    #         Algae_Opt_Mnth_NOx_StL.append(Constit_Loads_M['NO_Load_StL'].iloc[i])
    #         Algae_Opt_Mnth_Chla_StL.append(Constit_Loads_M['Chla_Load_StL'].iloc[i])
    #         Algae_Opt_Mnth_NOx_Cal.append(Constit_Loads_M['NO_Load_Cal'].iloc[i])
    #         Algae_Opt_Mnth_Chla_Cal.append(Constit_Loads_M['Chla_Load_Cal'].iloc[i])

    # Algae_Opt_Mnth_NOx_StL_arr = np.asarray(Algae_Opt_Mnth_NOx_StL)
    # Algae_Opt_Mnth_Chla_StL_arr = np.asarray(Algae_Opt_Mnth_Chla_StL)

    # Algae_Opt_Mnth_NOx_Cal_arr = np.asarray(Algae_Opt_Mnth_NOx_Cal)
    # Algae_Opt_Mnth_Chla_Cal_arr = np.asarray(Algae_Opt_Mnth_Chla_Cal)
    # if config['sim_type'] in [0, 1]:
    #     return[Constit_Loads_M,Nitro_Mod_Out_M,Algae_Opt_Mnth_NOx_StL_arr,Algae_Opt_Mnth_Chla_StL_arr,Algae_Opt_Mnth_NOx_Cal_arr,Algae_Opt_Mnth_Chla_Cal_arr]
    # else:
    #     return[Algae_Opt_Mnth_NOx_StL_arr,Algae_Opt_Mnth_Chla_StL_arr,Algae_Opt_Mnth_NOx_Cal_arr,Algae_Opt_Mnth_Chla_Cal_arr,Nitro_Model_Output]
   

if __name__ == "__main__":
    # Parse Arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "workspace",
        nargs=1,
        help="The path to the working directory.",
    )
    
    args = argparser.parse_args()
    workspace = args.workspace[0]
    
    # Run LOONE_Constituent_SimQ
    LOONE_Constituent_SimQ(workspace)
    