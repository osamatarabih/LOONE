# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 16:51:22 2020

@author: osama
"""
class LONINO_FNs:
    def RF_Cls(Wkly_Trib_Cond_NetRF):
        if Wkly_Trib_Cond_NetRF >= 8:
            Cls = 6
        elif 4 <= Wkly_Trib_Cond_NetRF < 8:
            Cls = 5
        elif 2 <= Wkly_Trib_Cond_NetRF < 4:
            Cls = 4
        elif -1 <= Wkly_Trib_Cond_NetRF < 2:
            Cls = 3
        elif -3 <= Wkly_Trib_Cond_NetRF < -1:
            Cls = 2
        elif -10 <= Wkly_Trib_Cond_NetRF < -3:
            Cls = 1
        elif Wkly_Trib_Cond_NetRF == -999999:
            Cls = 3
        else:
            Cls = 0
        return(Cls)
    def MainTrib_Cls(Wkly_Trib_Cond_S65E):
        if Wkly_Trib_Cond_S65E >= 9000:
            Cls = 6
        elif 6000 <= Wkly_Trib_Cond_S65E < 9000:
            Cls = 5
        elif 3500 <= Wkly_Trib_Cond_S65E < 6000:
            Cls = 4
        elif 500 <= Wkly_Trib_Cond_S65E < 3500:
            Cls = 3
        elif 200 <= Wkly_Trib_Cond_S65E < 500:
            Cls = 2
        elif 0 <= Wkly_Trib_Cond_S65E < 200:
            Cls = 1
        elif Wkly_Trib_Cond_S65E == -999999:
            Cls = 1
        else:
            Cls = 0
        return(Cls)
    def Palmer_Cls(Wkly_Trib_Cond_Palmer):
        if Wkly_Trib_Cond_Palmer >= 4:
            Cls = 6
        elif 2.9999 <= Wkly_Trib_Cond_Palmer < 4:
            Cls = 5
        elif 1.5 <= Wkly_Trib_Cond_Palmer < 2.9999:
            Cls = 4
        elif -1.4999 <= Wkly_Trib_Cond_Palmer < 1.5:
            Cls = 3
        elif -3 <= Wkly_Trib_Cond_Palmer < -1.4999:
            Cls = 2
        elif -5 <= Wkly_Trib_Cond_Palmer < -3:
            Cls = 1
        elif Wkly_Trib_Cond_Palmer == -999999:
            Cls = 1
        else:
            Cls = 0
        return(Cls)
    def NetInflow_Cls(Wkly_Trib_Cond_NetInf):
        if Wkly_Trib_Cond_NetInf >= 15000:
            Cls = 6
        elif 6000 <= Wkly_Trib_Cond_NetInf < 15000:
            Cls = 5
        elif 2500 <= Wkly_Trib_Cond_NetInf < 6000:
            Cls = 4
        elif 500 <= Wkly_Trib_Cond_NetInf < 2500:
            Cls = 3
        elif -5000 <= Wkly_Trib_Cond_NetInf < 500:
            Cls = 2
        elif -10000 <= Wkly_Trib_Cond_NetInf < -5000:
            Cls = 1
        elif Wkly_Trib_Cond_NetInf == -999999:
            Cls = 1
        else:
            Cls = 0
        return(Cls)
    def LONINO_Seas_cls(LONINO_df_LONINO_Seas):
        if LONINO_df_LONINO_Seas >= 2.0001:
            Cls = 4
        elif 1.8 <= LONINO_df_LONINO_Seas < 2.0001:
            Cls = 3
        elif 1.39 <= LONINO_df_LONINO_Seas < 1.8:
            Cls = 2
        elif -10 <= LONINO_df_LONINO_Seas < 1.39:
            Cls = 1
        elif LONINO_df_LONINO_Seas == -999999:
            Cls = 1
        else:
            Cls = 0
        return(Cls)
    def LONINO_M_Seas_cls(LONINO_df_LONINO_Mult_Seas):
        if LONINO_df_LONINO_Mult_Seas >= 2.0001:
            Cls = 4
        elif 1.18 <= LONINO_df_LONINO_Mult_Seas < 2.0001:
            Cls = 3
        elif 0.5001 <= LONINO_df_LONINO_Mult_Seas < 1.18:
            Cls = 2
        elif -10 <= LONINO_df_LONINO_Mult_Seas < 0.5001:
            Cls = 1
        elif LONINO_df_LONINO_Mult_Seas == -999999:
            Cls = 1
        else:
            Cls = 0
        return(Cls)