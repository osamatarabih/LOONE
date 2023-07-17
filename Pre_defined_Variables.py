# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 12:29:27 2020

@author: osamatarabih
"""
#A Class contains all the pre-defined variables
#I will include different values for variables in this class.
class Pre_defined_Variables:
    Schedule = 'LORS20082023'
    startyear = 2008
    endyear = 2023
    startdate_entry = startyear,1,1 
    begdateCS_entry = startyear,1,1 
    enddate_entry = endyear,3,31
    # enddate_TC = endyear+1,1,1
    enddate_TC = endyear,4,1
    # Month_N = (endyear-startyear+1)*12
    Month_N = (endyear-startyear)*12 + 3

    Opt_NewTree = 1 #if New Tree Decision is used enter 1 else enter 0.
    Code = 6
    Multiplier = 100
    TCI = 1
    Opt_NetInflow = 2
    NetInf_const = 0
    startstage = 10.268 #ft
    begstageCS = startstage
    Opt_LOSAdmd = 1
    Mult_LOSA = 100
    Opt_LOSAws = 1 #the option for LOSA Daily Supply where 1: Calculated Function of WSM, 2: Set values to zeros.
    Opt_DecTree = 1 #if Tree Decision is used enter 1 else enter 0.
    Zone_C_MetFcast_Indicator = 1 #0 to use the same tree classifications as SLONINO or 1 to use the same tree classifications as LORS2008 and SFWMM.
    WCA3a_REG_Zone = 'ERTP:TopE'
    WCA3Aoffset = 0
    WCA328min = 7.5
    WCA3NWmin = 11
    WCA217min = 11.1
    Opt_WCAlimitWSA = 2
    CSflag = 1 
    PlsDay_Switch = 0 #0: pulse day counter continues to 10 even if release level increases, 1: pulse day counter is set to zero if release level increases during the 10-day pulse.
    MaxQstgTrigger = 20 # the maximum stage trigger for maximum discharge if Trib_cond. = XWet.
    Opt_QregMult = 0 #option for using multipliers 0: don't use, 1: apply only during dry season (Nov-May), 2: apply year-round.
    AlternateHighQyrs = 0
    Option_S80Baseflow = 0 #0: baseflow supplements daily C44RO, 1: baseflow supplements monthly C44RO.
    S308BK_Const = 1 #method to calculate S308BK 0: Use Input data or 1: Simulate S308BK.
    S308_BK_Thr = 14.5 #threshold for S308 Backflow to occur.
    Opt_S308 = 1 #option for using S308 for LOREG and Backflow (1: Y and 0: N).
    S308RG_Const = 1 #0: use Input data (overrides Option_RegS77S308) and 1: simulate S308RG.
    Option_RegS77S308 = 0 #0: simulate discharge using LOONE, 1: same values as SFWMM LEC2005 Base, 2: Based on SFWMM 46yr run (see SFWMMdata4LOOPS).
    S80_Const = 1 #0: to use Input data or 1: to simulate S80.
    Opt_Outlet1DSRG = 0 #(Option for Regulatory Releases to CE) = where 0:  All releases >ZoneD0 measured at S77 and 1: All Pulse Releases measured at S79.
    THC_threshold = 2
    LowChance = 50
    Opt_LChance_line = 1
    Opt_Date_Targ_Stg = 1
    OptSalFcast = 3
    CE_SalThreshold = 5
    Late_Dry_Season_Option = 0
    Opt_NoAP_above_BF_SB = 1
    Opt_AdapProt = 1
    Opt_CEews_LOWSM = 0
    Opt_THCbypLateDS = 1
    APCB1 = 100
    APCB2 = 100
    APCB3 = 100
    APCB4 = 100
    CalEst_ews = 300
    Outlet1USEWS_Switch = 1
    Outlet1USBK_Switch = 1 #option for S77BK simulation 0: Use input data or 1: Simulate with LOONE.
    Outlet1USBK_Threshold = 11.1
    Option_S77Baseflow = 0 #0: baseflow supplements daily C43RO, 1: baseflow supplements monthly C43RO.
    Outlet1USREG_Switch = 1
    Outlet1DS_Switch = 1
    MaxCap_RegWCA = 4000
    Multiplier_RegWCA = 1
    Option_RegWCA = 2
    Constant_RegWCA = 400
    MaxCap_RegL8C51 = 500
    Multiplier_RegL8C51 = 1
    Option_RegL8C51 = 2
    Constant_RegL8C51 = 200
    ET_Switch = 0
    Opt_WSA = 0 #Options for Lake O water supply augmentation (WSA) operation (0 = no WSA,1 = use flat trigger stages to activate WSA operation,2 = trigger stages defined using offsets from LOWSM WST line to activate WSA operation)
    WSA_THC = 2
    WSAtrig1 = 12.5
    WSAtrig2 = 11.5
    WSAoff1 = 0.5
    WSAoff2 = 0
    MIAcap1 = 860
    MIAcap2 = 1720
    NNRcap1 = 900
    NNRcap2 = 1800
    Option_Stage = 0
    #Water demand cutback for each WSM Zone 
    Z1_cutback = 0.15
    Z2_cutback = 0.3
    Z3_cutback = 0.45
    Z4_cutback = 0.6
##################################################################################    
    dstar_B = 99        
    dstar_C = 99
    dstar_D3 = 99
    dstar_D2 = 99
    dstar_D1 = 99
    astar_B = 1
    astar_C = 1
    astar_D3 = 1
    astar_D2 = 1
    astar_D1 = 1
    bstar_S77_B = 1 
    bstar_S77_C = 1
    bstar_S77_D3 = 0.5
    bstar_S77_D2 = 0.5
    bstar_S77_D1 = 0.5
    bstar_S80_B = 1
    bstar_S80_C = 1
    bstar_S80_D3 = 0.5
    bstar_S80_D2 = 0.5
    bstar_S80_D1 = 0.5

    
    