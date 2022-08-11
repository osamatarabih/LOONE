import numpy as np
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 12:24:42 2020

@author: osamatarabih
"""

#This script runs the main functions of the LOONE Model
class LO_FNs:
    def WSM_Zone(Stage_LO,WSM4,WSM3,WSM2,WSM1): #Note that here WSM_Zone(i) return value at time step (i) in LO_Model dataFrame. 
        WSM = 0
        if Stage_LO < WSM4:
            WSM = 4
        elif Stage_LO < WSM3:
            WSM = 3
        elif Stage_LO < WSM2:
            WSM = 2
        elif Stage_LO < WSM1:
            WSM = 1
        else:
            WSM = 0
        return(WSM)
    def Max_Supply(WSM_Zone,LOSA_dmd_Dailydmd,Pre_defined_Variables_Z1_cutback,Pre_defined_Variables_Z2_cutback,Pre_defined_Variables_Z3_cutback,Pre_defined_Variables_Z4_cutback):
        M_S = 0
        if WSM_Zone == 1:
            M_S = LOSA_dmd_Dailydmd * (1-Pre_defined_Variables_Z1_cutback)
        elif WSM_Zone == 2:
            M_S = LOSA_dmd_Dailydmd * (1-Pre_defined_Variables_Z2_cutback)
        elif WSM_Zone == 3:
            M_S = LOSA_dmd_Dailydmd * (1-Pre_defined_Variables_Z3_cutback)
        elif WSM_Zone == 4:
            M_S = LOSA_dmd_Dailydmd * (1-Pre_defined_Variables_Z4_cutback)
        else:
            M_S = LOSA_dmd_Dailydmd
        return(M_S)
    def LOSA_Supply(WSM_Zone,LOSA_dmd,Max_Supply,Pre_defined_Variables_Opt_LOSAws):
        Sup = 0
        if Pre_defined_Variables_Opt_LOSAws == 1:
            if WSM_Zone == 0:
                Sup = LOSA_dmd
            else:
                Sup = min(LOSA_dmd, Max_Supply)
        elif Pre_defined_Variables_Opt_LOSAws == 2:
            Sup = 0
        else:
            pass
        return(Sup)
    #
    def Zone_Code(Stage_LO,WSMs_A,WSMs_B,WSMs_C,WSMs_D3,WSMs_D2,WSMs_D1,WSMs_D0,WSMs_WSM1):
        Zcd = 2
        if Stage_LO > WSMs_A:
            Zcd = 9
        elif Stage_LO > WSMs_B:
            Zcd = 8
        elif Stage_LO > WSMs_C:
            Zcd = 7
        elif Stage_LO > WSMs_D3:
            Zcd = 6
        elif Stage_LO > WSMs_D2:
            Zcd = 5
        elif Stage_LO > WSMs_D1:
            Zcd = 4
        elif Stage_LO > WSMs_D0:
            Zcd = 3
        elif Stage_LO < WSMs_WSM1:
            Zcd = 1
        else:
            Zcd = 2
        return(Zcd)
    def LO_Zone(Zone_Code):
        Z = 0
        if Zone_Code == 1:
            Z = 'SSM'
        elif Zone_Code == 2:
            Z = 'E'
        elif Zone_Code == 3:
            Z = 'D0'
        elif Zone_Code == 4:
            Z = 'D1'
        elif Zone_Code == 5:
            Z = 'D2'
        elif Zone_Code == 6:
            Z = 'D3'
        elif Zone_Code == 7:
            Z = 'C'
        elif Zone_Code == 8:
            Z = 'B'
        elif Zone_Code == 9:
            Z = 'A' 
        return(Z)
    def DecTree_Relslevel(Zone_Code,Zone_D_Rel_Code,Zone_C_Rel_Code,Zone_B_Rel_Code): #Note that DecTree_Relslevel(i) returns value at time step (i+2) in LO_Model dataFrame.
        DTR = 0
        if Zone_Code == 1:
            DTR = 0
        elif Zone_Code == 2:
            DTR = 0
        elif Zone_Code == 3:
            DTR = -1
        elif Zone_Code == 4:
            DTR = Zone_D_Rel_Code
        elif Zone_Code == 5:
            if Zone_D_Rel_Code == -1:
                DTR = -1
            else: 
                DTR = min(4, Zone_D_Rel_Code*2)
        elif Zone_Code == 6:
            if Zone_D_Rel_Code == -1:
                DTR = -1
            else: 
                DTR = min(4, Zone_D_Rel_Code*3)
        elif Zone_Code == 7:
            DTR = Zone_C_Rel_Code
        elif Zone_Code == 8:
            DTR = Zone_B_Rel_Code
        else:
            DTR = 6
        return(DTR)
    def PlsDay(DayFlags,DecTree_Relslevel,Pre_defined_Variables_PlsDay_Switch):
        PlsD = 0
        if DayFlags == 'P.A.Day1':
            if DecTree_Relslevel > 0 and DecTree_Relslevel < 4:
                PlsD = 1
            else:
                PlsD = 0
        else:
            if Pre_defined_Variables_PlsDay_Switch == 1:
                if 'PlsD' in locals():
                    if PlsD > 0 and PlsD < 10 and DecTree_Relslevel < 4:
                        PlsD = PlsD + 1
                    elif DecTree_Relslevel > 0 and DecTree_Relslevel < 4:
                        PlsD = 1
                    else:
                        PlsD = 0
                else:
                    PlsD = 0
            elif Pre_defined_Variables_PlsDay_Switch == 0:
                if 'PlsD' in locals():
                    if PlsD > 0 and PlsD < 10:
                        PlsD = PlsD + 1
                    elif DecTree_Relslevel > 0 and DecTree_Relslevel < 4:
                        PlsD = 1
                    else:
                        PlsD = 0
                else:
                    PlsD = 0
        return(PlsD)
    def Release_Level(Prev_RL,Stage_LO,Tributary_Condition,PlsDay,Zone_Code,DecTree_Relslevel,Pre_defined_Variables_MaxQstgTrigger):
        RL = Prev_RL
        if Stage_LO > Pre_defined_Variables_MaxQstgTrigger and Tributary_Condition == 6:#Note that at i =2 then, (i-2) = 0 which indicates 1/1/1965 for TC_LONINO because it starts at 1/1/1965! 
            RL = 6
        elif PlsDay == 0 or Zone_Code > 7:
            RL = DecTree_Relslevel
        elif PlsDay == 1:
            RL = DecTree_Relslevel
        else:
            RL = RL
        return(RL)
    def ZoneCodeminus1Code(Zone_Code,WSMs_WSM1,WSMs_D0,WSMs_D1,WSMs_D2,WSMs_D3,WSMs_C,WSMs_B,WSMs_A): 
        Sam = 0
        if (Zone_Code - 1) == 1: 
            Sam = WSMs_WSM1 
        elif (Zone_Code - 1) == 2: 
            Sam = WSMs_D0
        elif (Zone_Code - 1) == 3:
            Sam = WSMs_D1 
        elif (Zone_Code - 1) == 4: 
            Sam = WSMs_D2
        elif (Zone_Code - 1) == 5:
            Sam = WSMs_D3 
        elif (Zone_Code - 1) == 6:
            Sam = WSMs_C
        elif (Zone_Code - 1) == 7: 
            Sam = WSMs_B 
        elif (Zone_Code - 1) == 8:
            Sam = WSMs_A
        return(Sam)
    def ZoneCodeCode(Zone_Code,WSMs_WSM1,WSMs_D0,WSMs_D1,WSMs_D2,WSMs_D3,WSMs_C,WSMs_B,WSMs_A):
        Roty = 0
        if Zone_Code == 1: 
            Roty = WSMs_WSM1 
        elif Zone_Code == 2: 
            Roty = WSMs_D0
        elif Zone_Code == 3:
            Roty = WSMs_D1 
        elif Zone_Code == 4: 
            Roty = WSMs_D2
        elif Zone_Code == 5:
            Roty = WSMs_D3 
        elif Zone_Code == 6:
            Roty = WSMs_C
        elif Zone_Code == 7: 
            Roty = WSMs_B 
        elif Zone_Code == 8:
            Roty = WSMs_A
        return(Roty)
    def Fraction_of_Zone_height(Zone_Code,Stage_LO,ZoneCodeminus1Code,ZoneCodeCode):
        FZH = 0
        if Zone_Code == 1 or Zone_Code == 9:
            FZH = 1
        else:
            FZH = ((Stage_LO - ZoneCodeminus1Code) / (ZoneCodeCode - ZoneCodeminus1Code + 0.01))
        return(FZH)
    def ReLevelCode_1(Release_Level,Pre_defined_Variables_dstar_D1,Pre_defined_Variables_dstar_D2,Pre_defined_Variables_dstar_D3,Pre_defined_Variables_dstar_C,Pre_defined_Variables_dstar_B):
        Hnur = 0
        if (Release_Level + 2) == 1:
            Hnur = -99
        elif (Release_Level + 2) == 2:
            Hnur = -99
        elif (Release_Level + 2) == 3:
            Hnur = Pre_defined_Variables_dstar_D1
        elif (Release_Level + 2) == 4:
            Hnur = Pre_defined_Variables_dstar_D2
        elif (Release_Level + 2) == 5:
            Hnur = Pre_defined_Variables_dstar_D3
        elif (Release_Level + 2) == 6:
            Hnur = Pre_defined_Variables_dstar_C
        elif (Release_Level + 2) == 7:
            Hnur = Pre_defined_Variables_dstar_B
        elif (Release_Level + 2) == 8:
            Hnur = -99
        return(Hnur)
    def ReLevelCode_2(Release_Level,Pre_defined_Variables_astar_D1,Pre_defined_Variables_astar_D2,Pre_defined_Variables_astar_D3,Pre_defined_Variables_astar_C,Pre_defined_Variables_astar_B):
        Hnur = 0
        if (Release_Level + 2) == 1:
            Hnur = 0
        elif (Release_Level + 2) == 2:
            Hnur = 0
        elif (Release_Level + 2) == 3:
            Hnur = Pre_defined_Variables_astar_D1
        elif (Release_Level + 2) == 4:
            Hnur = Pre_defined_Variables_astar_D2
        elif (Release_Level + 2) == 5:
            Hnur = Pre_defined_Variables_astar_D3
        elif (Release_Level + 2) == 6:
            Hnur = Pre_defined_Variables_astar_C
        elif (Release_Level + 2) == 7:
            Hnur = Pre_defined_Variables_astar_B
        elif (Release_Level + 2) == 8:
            Hnur = 0
        return(Hnur)
    def ReLevelCode_3_S80(Release_Level,Pre_defined_Variables_bstar_S80_D1,Pre_defined_Variables_bstar_S80_D2,Pre_defined_Variables_bstar_S80_D3,Pre_defined_Variables_bstar_S80_C,Pre_defined_Variables_bstar_S80_B):
        Hnur = 0
        if (Release_Level + 2) == 1:
            Hnur = 1
        elif (Release_Level + 2) == 2:
            Hnur = 2
        elif (Release_Level + 2) == 3:
            Hnur = Pre_defined_Variables_bstar_S80_D1
        elif (Release_Level + 2) == 4:
            Hnur = Pre_defined_Variables_bstar_S80_D2
        elif (Release_Level + 2) == 5:
            Hnur = Pre_defined_Variables_bstar_S80_D3
        elif (Release_Level + 2) == 6:
            Hnur = Pre_defined_Variables_bstar_S80_C
        elif (Release_Level + 2) == 7:
            Hnur = Pre_defined_Variables_bstar_S80_B
        elif (Release_Level + 2) == 8:
            Hnur = 1
        return(Hnur)
    def Outlet2DS_Mult(Seasons_Season,Seasons_Month,dh_7days,ReLevelCode_1,Fraction_of_Zone_height,ReLevelCode_2,ReLevelCode_3_S80,Pre_defined_Variables_Opt_QregMult):
        Mult = 0
        if Pre_defined_Variables_Opt_QregMult == 0 or (Pre_defined_Variables_Opt_QregMult == 1 and Seasons_Season > 2) or\
        (Pre_defined_Variables_Opt_QregMult == 2 and Seasons_Month > 4 and Seasons_Month < 9):
            Mult = 1
        elif dh_7days < ReLevelCode_1 and Fraction_of_Zone_height < ReLevelCode_2:
            Mult = ReLevelCode_3_S80
        return(Mult)
    def Outlet2DS_Mult_2(month,day,PlsDay,S80_Mult_i_PlsDay,S80_Mult_i,Pre_defined_Variables_Opt_QregMult):
        S80_M2 = 0
        if Pre_defined_Variables_Opt_QregMult == 1 and (month == 6 or month == 11)\
        and day < 10 and PlsDay > 1 and PlsDay < 11:
            S80_M2 = S80_Mult_i_PlsDay
        else:
            S80_M2 = S80_Mult_i
        return(S80_M2)
    
    def Outlet2DSRS(Release_Level,S80_RegRelRates_Zone_D1,S80avgL1,Pulses_S_80_L1,S80_Mult_2,CE_SLE_turns_SLEturn,S80_RegRelRates_Zone_D2,S80avgL2,Pulses_S_80_L2,S80_RegRelRates_Zone_D3,S80avgL3,Pulses_S_80_L3,S80_RegRelRates_Zone_C,S80_RegRelRates_Zone_B,S80_RegRelRates_Zone_A):
        S = 0
        if (Release_Level + 2) ==1:
            S = 0
        elif (Release_Level + 2) ==2:
            S = 0
        elif (Release_Level + 2) ==3:
            S = S80_RegRelRates_Zone_D1/S80avgL1*Pulses_S_80_L1*S80_Mult_2*CE_SLE_turns_SLEturn
        elif (Release_Level + 2) ==4:
            S = S80_RegRelRates_Zone_D2/S80avgL2*Pulses_S_80_L2*S80_Mult_2*CE_SLE_turns_SLEturn
        elif (Release_Level + 2) ==5:
            S = S80_RegRelRates_Zone_D3/S80avgL3*Pulses_S_80_L3*S80_Mult_2*CE_SLE_turns_SLEturn
        elif (Release_Level + 2) ==6:
            S = S80_RegRelRates_Zone_C*S80_Mult_2*CE_SLE_turns_SLEturn
        elif (Release_Level + 2) ==7:
            S = S80_RegRelRates_Zone_B*S80_Mult_2*CE_SLE_turns_SLEturn
        elif (Release_Level + 2) ==8:
            S = S80_RegRelRates_Zone_A
        return(S)
    
    def Sum_Outlet2USRG1(day,S308RG1):
        Sigma = 0
        if day ==1:
            Sigma = S308RG1
        else:
            Sigma = S308RG1 + Sigma
        return(Sigma)
    
    def Outlet2DSBS(Release_Level,Sum_S308RG1,VLOOKUP1_c,StlEst_baseflow,Pre_defined_Variables_Option_S80Baseflow):
        BS = 0
        if Release_Level == -1:
            if Pre_defined_Variables_Option_S80Baseflow ==0:
                BS = StlEst_baseflow
            elif Sum_S308RG1/31 > VLOOKUP1_c:
                BS = 0
            else:
                BS = VLOOKUP1_c
        else:
            BS = 0
        return(BS)
    
    def Outlet2USBK(Stage_LO,WSMs_D1,S308RG,C44RO,S308BK_data,Pre_defined_Variables_Opt_S308,Pre_defined_Variables_S308BK_Const,Pre_defined_Variables_S308_BK_Thr):
        S308 = 0
        if Pre_defined_Variables_Opt_S308 ==1:
            if Pre_defined_Variables_S308BK_Const == 1:
                if Stage_LO < min(Pre_defined_Variables_S308_BK_Thr, WSMs_D1-0.25) and S308RG == 0:
                    S308 = C44RO
                else:
                    S308 = 0
            else:
                S308 = S308BK_data
        else:
            S308 = 0
        return(S308)
    
    def Outlet2USBS(S80BS,S308RG1,ROeast,Pre_defined_Variables_Option_S80Baseflow):
        S308 = 0
        if Pre_defined_Variables_Option_S80Baseflow == 0:
            S308 = max(0, S80BS-(S308RG1+ROeast))
        else:
            S308 = max(0, S80BS-S308RG1)
        return(S308)
    
    def Sum_Outlet2USBK(day,S308BK):
        SBK = 0
        if day ==1:
            SBK = S308BK
        else:
            SBK = S308BK + SBK
        return(SBK)
    def Outlet2USRG_Code(S308RG1,S308BS,S308RG_data,STEST_data,Pre_defined_Variables_Option_RegS77S308):
        Code_1 = 0
        if Pre_defined_Variables_Option_RegS77S308 + 1 == 1:
            Code_1 = S308RG1 + S308BS
        elif Pre_defined_Variables_Option_RegS77S308 + 1 == 2:
            Code_1 = 0
        elif Pre_defined_Variables_Option_RegS77S308 + 1 == 3:
            Code_1 = S308RG_data + STEST_data
        elif Pre_defined_Variables_Option_RegS77S308 + 1 == 4:
            Code_1 = 0
        return Code_1
    def Outlet2USRG(S308RG_Code,S308RG_data,STEST_data,Pre_defined_Variables_Opt_S308,Pre_defined_Variables_S308RG_Const):
        S_1 = 0
        if Pre_defined_Variables_Opt_S308 ==1:
            if Pre_defined_Variables_S308RG_Const == 1:
                S_1 = S308RG_Code
            else:
                S_1 = S308RG_data + STEST_data
        else:
            S_1 = 0
        return(S_1)
    def S80(ROeast,S308RG,S80_data,Pre_defined_Variables_S80_Const):
        S_80 = 0
        if Pre_defined_Variables_S80_Const == 1:
            S_80 = ROeast + S308RG
        else: 
            S_80 = S80_data
        return(S_80)
    def ReLevelCode_3_S77(Release_Level,Pre_defined_Variables_bstar_S77_D1,Pre_defined_Variables_bstar_S77_D2,Pre_defined_Variables_bstar_S77_D3,Pre_defined_Variables_bstar_S77_C,Pre_defined_Variables_bstar_S77_B):
        Hnur = 0
        if (Release_Level + 2) == 1:
            Hnur = 1
        elif (Release_Level + 2) == 2:
            Hnur = 1
        elif (Release_Level + 2) == 3:
            Hnur = Pre_defined_Variables_bstar_S77_D1
        elif (Release_Level + 2) == 4:
            Hnur = Pre_defined_Variables_bstar_S77_D2
        elif (Release_Level + 2) == 5:
            Hnur = Pre_defined_Variables_bstar_S77_D3
        elif (Release_Level + 2) == 6:
            Hnur = Pre_defined_Variables_bstar_S77_C
        elif (Release_Level + 2) == 7:
            Hnur = Pre_defined_Variables_bstar_S77_B
        elif (Release_Level + 2) == 8:
            Hnur = 1
        return(Hnur)
    def Outlet1US_Mult(Seasons_Season,Seasons_Month,dh_7days,ReLevelCode_1,Fraction_of_Zone_height,ReLevelCode_2,ReLevelCode_3_S77,Pre_defined_Variables_Opt_QregMult):
        Mult = 0
        if Pre_defined_Variables_Opt_QregMult == 0 or (Pre_defined_Variables_Opt_QregMult == 1 and Seasons_Season > 2) or\
        (Pre_defined_Variables_Opt_QregMult == 2 and Seasons_Month > 4 and Seasons_Month < 9):
            Mult = 1
        elif dh_7days < ReLevelCode_1 and Fraction_of_Zone_height < ReLevelCode_2:
            Mult = ReLevelCode_3_S77
        return(Mult)
    def Outlet1US_Mult_2(month,day,PlsDay,S77_Mult_i_PlsDay,S77_Mult_i,Pre_defined_Variables_Opt_QregMult):
        S77_M2 = 0
        if Pre_defined_Variables_Opt_QregMult == 1 and (month == 6 or month == 11)\
        and day < 10 and PlsDay > 1 and PlsDay < 11:
            S77_M2 = S77_Mult_i_PlsDay
        else:
            S77_M2 = S77_Mult_i
        return(S77_M2)
    def Outlet1USRS(Release_Level,S77_RegRelRates_Zone_D1,S77avgL1,Pulses_S_77_L1,S77_Mult_2,C43RO,CE_SLE_turns_CEturn,S77_RegRelRates_Zone_D2,S77avgL2,Pulses_S_77_L2,Zone_Code,S77_RegRelRates_Zone_D3,S77avgL3,Pulses_S_77_L3,S77_RegRelRates_Zone_C,S77_RegRelRates_Zone_B,S77_RegRelRates_Zone_A,Pre_defined_Variables_Opt_Outlet1DSRG):
        S = 0
        if (Release_Level + 2) ==1:
            S = 0
        elif (Release_Level + 2) ==2:
            S = 0
        elif (Release_Level + 2) ==3:
            S = max(0, S77_RegRelRates_Zone_D1/S77avgL1*Pulses_S_77_L1*S77_Mult_2 - Pre_defined_Variables_Opt_Outlet1DSRG * C43RO) * CE_SLE_turns_CEturn
        elif (Release_Level + 2) ==4:
            S = max(0, S77_RegRelRates_Zone_D2/S77avgL2*Pulses_S_77_L2*S77_Mult_2 - Pre_defined_Variables_Opt_Outlet1DSRG * C43RO) * CE_SLE_turns_CEturn
        elif (Release_Level + 2) ==5 and Zone_Code < 7:
            S = max(0, S77_RegRelRates_Zone_D3/S77avgL3*Pulses_S_77_L3*S77_Mult_2 - Pre_defined_Variables_Opt_Outlet1DSRG * C43RO) * CE_SLE_turns_CEturn
        elif (Release_Level + 2) ==5 and Zone_Code >= 7:
            S = max(0, S77_RegRelRates_Zone_D3/S77avgL3*Pulses_S_77_L3*S77_Mult_2) * CE_SLE_turns_CEturn
        elif (Release_Level + 2) ==6:
            S = S77_RegRelRates_Zone_C*S77_Mult_2*CE_SLE_turns_CEturn
        elif (Release_Level + 2) ==7:
            S = S77_RegRelRates_Zone_B*S77_Mult_2*CE_SLE_turns_CEturn
        elif (Release_Level + 2) ==8:
            S = S77_RegRelRates_Zone_A
        return(S)
    def Sum_Outlet1USRS(day,Outlet1USRS):
        Sigma = 0
        if day ==1:
            Sigma = Outlet1USRS
        else:
            Sigma = Outlet1USRS + Sigma
        return(Sigma)
    def Outlet1USBK(Stage_LO,Outlet1USRS,Outlet1USBSAP,Outlet1USEWS,C43RO,S77BK_data,Pre_defined_Variables_Outlet1USBK_Switch,Pre_defined_Variables_Outlet1USBK_Threshold):
        Y = 0
        if Pre_defined_Variables_Outlet1USBK_Switch == 1:
            if Stage_LO < Pre_defined_Variables_Outlet1USBK_Threshold and Outlet1USRS == 0 and\
            Outlet1USBSAP == 0 and Outlet1USEWS == 0:
                Y = max(0,0.3*C43RO-5)
            else:
                Y = 0
        else:
            Y = S77BK_data
        return(Y)
    # def Outlet1DSBS(Release_Level,Sum_Outlet1USRS,VLOOKUP2_c,CalEst_baseflow,Pre_defined_Variables_Option_S77Baseflow):
    #     BS = 0
    #     if Release_Level == -1:
    #         if Pre_defined_Variables_Option_S77Baseflow ==0:
    #             BS = CalEst_baseflow
    #         elif Sum_Outlet1USRS/31 > VLOOKUP2_c: 
    #             BS = 0
    #         else:
    #             BS = VLOOKUP2_c 
    #     else:
    #         BS = 0
    #     return(BS)
    def Outlet1DSBS(Release_Level,Sum_Outlet1USRS,VLOOKUP2_c,CalEst_baseflow,Pre_defined_Variables_Option_S77Baseflow):
        BS = 0
        if Release_Level >= -1:
            if Pre_defined_Variables_Option_S77Baseflow ==0:
                BS = CalEst_baseflow
            elif Sum_Outlet1USRS/31 > VLOOKUP2_c: 
                BS = 0
            else:
                BS = VLOOKUP2_c 
        else:
            BS = 0
        return(BS)
    # def Outlet1DSBS(Release_Level,Sum_Outlet1USRS,VLOOKUP2_c,CalEst_baseflow,Pre_defined_Variables_Option_S77Baseflow):
    #     BS = 0
    #     if Release_Level == -1 or Release_Level == 1 or Release_Level == 4: 
    #         if Pre_defined_Variables_Option_S77Baseflow ==0:
    #             BS = CalEst_baseflow
    #         elif Sum_Outlet1USRS/31 > VLOOKUP2_c: 
    #             BS = 0
    #         else:
    #             BS = VLOOKUP2_c 
    #     else:
    #         BS = 0
    #     return(BS)
    def Outlet1USBS(Outlet1DSBS,Outlet1USRS,ROwest,Pre_defined_Variables_Option_S77Baseflow):
        Z = 0
        if Pre_defined_Variables_Option_S77Baseflow == 0:
            Z = max(0,Outlet1DSBS-(Outlet1USRS+ROwest))
        else:
            Z = max(0,Outlet1DSBS-Outlet1USRS)
        return(Z)
    def Outlet1USBSAP(Outlet1USBS,Post_Ap_Baseflow,Pre_defined_Variables_Opt_AdapProt):
        S77 = 0
        if Pre_defined_Variables_Opt_AdapProt == 0:
            S77 = Outlet1USBS
        else:
            S77 = Post_Ap_Baseflow
        return(S77)
    def Outlet1USEWS(Post_AP_EWS,CAEST_data,Pre_defined_Variables_Outlet1USEWS_Switch,Pre_defined_Variables_Opt_AdapProt):
        X = 0
        if Pre_defined_Variables_Outlet1USEWS_Switch ==1:
            if Pre_defined_Variables_Opt_AdapProt == 0:
                X = 0
            else:
                X = Post_AP_EWS
        else:
            X = CAEST_data
        return(X)
    def Outlet1USREG(Outlet1USRS,Outlet1USBSAP,S77RG_data,Pre_defined_Variables_Outlet1USREG_Switch,Pre_defined_Variables_Option_RegS77S308):
        W = 0
        if Pre_defined_Variables_Outlet1USREG_Switch == 1:
            if Pre_defined_Variables_Option_RegS77S308 + 1 == 1:
                W = Outlet1USRS + Outlet1USBSAP
            elif Pre_defined_Variables_Option_RegS77S308 + 1 == 2:
                W = 0
            elif Pre_defined_Variables_Option_RegS77S308 + 1 == 3:
                W = S77RG_data
            elif Pre_defined_Variables_Option_RegS77S308 + 1 == 4:
                W = 0
        else:
            W = S77RG_data
        return(W)
    def Outlet1DS(Outlet1USREG,Outlet1USEWS,ROwest,Outlet1DS_data,Pre_defined_Variables_Outlet1DS_Switch):
        S = 0
        if Pre_defined_Variables_Outlet1DS_Switch == 1:
            S = Outlet1USREG + Outlet1USEWS + ROwest
        else:
            S = Outlet1DS_data
        return(S)
    def Choose_WCA(RegWCA_data,Pre_defined_Variables_Option_RegWCA,Pre_defined_Variables_Constant_RegWCA):
        Ch = 0
        if Pre_defined_Variables_Option_RegWCA == 1:
            Ch = 0
        elif Pre_defined_Variables_Option_RegWCA == 2:
            Ch = RegWCA_data
        elif Pre_defined_Variables_Option_RegWCA == 3:
            Ch = 0
        elif Pre_defined_Variables_Option_RegWCA == 4:
            Ch = Pre_defined_Variables_Constant_RegWCA
        return(Ch)
    def Choose_L8C51(RegL8C51_data,Pre_defined_Variables_Option_RegL8C51,Pre_defined_Variables_Constant_RegL8C51):
        Ch = 0
        if Pre_defined_Variables_Option_RegL8C51 == 1:
            Ch = 0
        elif Pre_defined_Variables_Option_RegL8C51 == 2:
            Ch = RegL8C51_data
        elif Pre_defined_Variables_Option_RegL8C51 == 3:
            Ch = 0
        elif Pre_defined_Variables_Option_RegL8C51 == 4:
            Ch = Pre_defined_Variables_Constant_RegL8C51
        return(Ch)
    def ET(et_dry_data,Stage2ar,et_litoral_data,Stage2marsh,et_open_data,ETVOL_data,Pre_defined_Variables_ET_Switch):
        E = 0
        if Pre_defined_Variables_ET_Switch == 1:
            E = ((et_dry_data/12)*(466000-Stage2ar)) + ((et_litoral_data/12)*(Stage2marsh))+\
            ((et_open_data/12)*(Stage2ar-Stage2marsh))
        else:
            E = ETVOL_data
        return(E)
    def Choose_WSA_1(WSMs_WSM1,Pre_defined_Variables_Opt_WSA,Pre_defined_Variables_WSAtrig2,Pre_defined_Variables_WSAoff2):
        Ch1 = 0
        if Pre_defined_Variables_Opt_WSA == 1:
            Ch1 = Pre_defined_Variables_WSAtrig2
        elif Pre_defined_Variables_Opt_WSA == 2:
            Ch1 = WSMs_WSM1 + Pre_defined_Variables_WSAoff2
        else:
            Ch1 = np.nan
        return(Ch1)
    def Choose_WSA_2(WSMs_WSM1,Pre_defined_Variables_Opt_WSA,Pre_defined_Variables_WSAtrig1,Pre_defined_Variables_WSAoff1):
        Ch2 = 0
        if Pre_defined_Variables_Opt_WSA == 1:
            Ch2 = Pre_defined_Variables_WSAtrig1
        elif Pre_defined_Variables_Opt_WSA == 2:
            Ch2 = WSMs_WSM1 + Pre_defined_Variables_WSAoff1
        else:
            Ch2 = np.nan
        return(Ch2)
    def WSA_MIA(Are_WCA_stages_too_low,LONINO_Seasonal_Classes,Stage_LO,Choose_WSA_1,MIA_cfs,S3PMP,Choose_WSA_2,Pre_defined_Variables_Opt_WSA,Pre_defined_Variables_WSA_THC,Pre_defined_Variables_MIAcap2,Pre_defined_Variables_MIAcap1):
        WSA = 0
        if Pre_defined_Variables_Opt_WSA == 0 or Are_WCA_stages_too_low == True or LONINO_Seasonal_Classes > Pre_defined_Variables_WSA_THC:
            WSA = 0
        elif Stage_LO < Choose_WSA_1:
            WSA = max(0,min(Pre_defined_Variables_MIAcap2, MIA_cfs)-S3PMP)
        elif Stage_LO < Choose_WSA_2:
            WSA = max(0,min(Pre_defined_Variables_MIAcap1, MIA_cfs)-S3PMP)
        else:
            WSA = 0
        return(WSA)
    def WSA_NNR(Are_WCA_stages_too_low,LONINO_Seasonal_Classes,Stage_LO,Choose_WSA_1,NNR_cfs,S2PMP,Choose_WSA_2,Pre_defined_Variables_Opt_WSA,Pre_defined_Variables_WSA_THC,Pre_defined_Variables_NNRcap2,Pre_defined_Variables_NNRcap1):
        WSA = 0
        if Pre_defined_Variables_Opt_WSA == 0 or Are_WCA_stages_too_low == True or LONINO_Seasonal_Classes > Pre_defined_Variables_WSA_THC:
            WSA = 0
        elif Stage_LO < Choose_WSA_1:
            WSA = max(0,min(Pre_defined_Variables_NNRcap2, NNR_cfs)-S2PMP)
        elif Stage_LO < Choose_WSA_2:
            WSA = max(0,min(Pre_defined_Variables_NNRcap1, NNR_cfs)-S2PMP)
        else:
            WSA = 0
        return(WSA)
    def Storage(DayFlags,Storage_i,StartStorage,Storage_iplus1,DSto_iplus2):
        Sto = 0
        if DayFlags == 'CS start date':
            Sto = Storage_i
        elif DayFlags == 'SimDay1':
            Sto = StartStorage
        else:
            Sto = Storage_iplus1 + DSto_iplus2
        return(Sto)
    def Lake_Stage(stg2sto_Storage_iplus2,EOD_Stg_data,Pre_defined_Variables_Option_Stage):
        Sta = 0
        if Pre_defined_Variables_Option_Stage+1 == 1:
            Sta = stg2sto_Storage_iplus2
        elif Pre_defined_Variables_Option_Stage+1 == 2:
            Sta = np.nan
        elif Pre_defined_Variables_Option_Stage+1 == 3:
            Sta = EOD_Stg_data    
        elif Pre_defined_Variables_Option_Stage+1 == 4:
            Sta = np.nan
        elif Pre_defined_Variables_Option_Stage+1 == 5:
            Sta = np.nan
        return(Sta)
                            