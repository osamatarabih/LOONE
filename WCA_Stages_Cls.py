# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 18:30:33 2022

@author: osama
"""
#This following section calculates the parameter states (trib conds, stage tests, seasonal & multi-seasonal LONINO) and sets a 4-digit code for the branch.  The branch code is used with the Routing sheet to determine release rates.
def WCA_Stages_Cls(TC_LONINO_df):
    import pandas as pd
    from datetime import datetime
    import numpy as np
    from scipy import interpolate
    from calendar import monthrange  
    import os
    from Model_Config import Model_Config 
    Working_Path = Model_Config.Working_Path
    os.chdir('%s'%Working_Path) 
    from Pre_defined_Variables import Pre_defined_Variables 
    from Additional_Fncs import Add_Fn
    
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    startdate = datetime(year, month, day).date() 
    year, month, day = map(int, Pre_defined_Variables.startdate_entry)
    begdateCS = datetime(year, month, day).date()
    year, month, day = map(int, Pre_defined_Variables.enddate_entry)
    enddate = datetime(year, month, day).date()
    
    date_rng_5 = pd.date_range(start = startdate, end = enddate, freq ='D')
    
    #Read the WCA Stage data
    WCA_Stages = pd.read_csv('./Data/%s/%s/WCA_Stages_Inputs_%s.csv'% (Pre_defined_Variables.Schedule,'ts_data',Pre_defined_Variables.Schedule))
    #Read WCA3A_REG inputs 
    #Note that I added a date column in the first column and I copied values of Feb28 to the Feb29!
    WCA3A_REG = pd.read_csv('./Data/%s/%s/WCA3A_REG_Inputs_%s.csv'% (Pre_defined_Variables.Schedule,'fixed_data',Pre_defined_Variables.Schedule))

    #generate WCA Stage dataframe 
    WCA_Stages_df = pd.DataFrame(date_rng_5, columns = ['Date'])
    WCA_Stages_df['3A-NW'] = WCA_Stages['3A-NW']
    WCA_Stages_df['2A-17'] = WCA_Stages['2A-17']
    WCA_Stages_df['3A-3'] = WCA_Stages['3A-3']
    WCA_Stages_df['3A-4'] = WCA_Stages['3A-4']
    WCA_Stages_df['3A-28'] = WCA_Stages['3A-28']
    WCA_Stages_df['3gage Avg'] =(WCA_Stages_df['3A-3'] + WCA_Stages_df['3A-4'] + WCA_Stages_df['3A-28'])/3
    #Generate a daily date range for one year (2020)
    date_rng_10 = pd.date_range(start = '1/1/2020', end = '12/31/2020', freq ='D') 
    WCA3A_REG['Date'] = date_rng_10
    Simperioddays = len(WCA_Stages_df.index)
    Oneyeardays = len(WCA3A_REG.index)
    #This Following Loop assigns the 366 values of the predifined WCA3a_REG_Zone (One random year) in the WCA3A_REG dataframe 
    #to the entire study period.
    def Replicate_WCA3 (year, day_num): #Where x is the percentage value (i.e., 10,20,30,40,50,60)
        leap_day_val = WCA3A_REG['%s'%Pre_defined_Variables.WCA3a_REG_Zone].iloc[59]
        if Add_Fn.leap_year(year) == True:
            day_num_adj = day_num
        else:
            day_num_adj = day_num + (1 if day_num >= 60 else 0)
        day_value = leap_day_val if day_num_adj == 60 and Add_Fn.leap_year(year) == True else WCA3A_REG['%s'%Pre_defined_Variables.WCA3a_REG_Zone].iloc[day_num_adj-1]
        return(day_value) 
    Iden_WCA3Areg = np.zeros(Simperioddays)
    for i in range(Simperioddays):
        Iden_WCA3Areg[i] = Replicate_WCA3(WCA_Stages_df['Date'].iloc[i].year, WCA_Stages_df['Date'].iloc[i].timetuple().tm_yday)  
    Iden_WCA3Areg_c = [x for x in Iden_WCA3Areg if ~np.isnan(x)]
    WCA_Stages_df['%s'%Pre_defined_Variables.WCA3a_REG_Zone] = Iden_WCA3Areg_c
    #Define WCA-2A Schedule
    def Replicate_WCA2A (year, day_num): #Where x is the percentage value (i.e., 10,20,30,40,50,60)
        leap_day_val = WCA3A_REG['WCA2Areg'].iloc[59]
        if Add_Fn.leap_year(year) == True:
            day_num_adj = day_num
        else:
            day_num_adj = day_num + (1 if day_num >= 60 else 0)
        day_value = leap_day_val if day_num_adj == 60 and Add_Fn.leap_year(year) == True else WCA3A_REG['WCA2Areg'].iloc[day_num_adj-1]
        return(day_value) 
    WCA_2A_Sched = np.zeros(Simperioddays)
    for i in range(Simperioddays):
        WCA_2A_Sched[i] = Replicate_WCA2A(WCA_Stages_df['Date'].iloc[i].year, WCA_Stages_df['Date'].iloc[i].timetuple().tm_yday)  
    WCA_2A_Sched_c = [x for x in WCA_2A_Sched if ~np.isnan(x)]
    WCA_Stages_df['WCA-2A Schedule'] = WCA_2A_Sched_c
    #The Following Lines Determine the predifined WCA3A_REG Zone + Offset, and determine some TRUE or
    #FALSE Logical Tests in #SFWMMdata4LOOPS Tab in the Spreadsheet LOOPS Model.
    Iden_WCA3Areg_plus_offset = np.zeros(Simperioddays)
    thr_gag_les_sch_pl_off = np.zeros(Simperioddays)
    thrA28_less_min = np.zeros(Simperioddays)
    thrANW_less_min = np.zeros(Simperioddays)
    twoA17_less_min = np.zeros(Simperioddays)
    WCA_3A_Stg_gr_Iden_WCA3Areg = np.zeros(Simperioddays)
    WCA_2A_Stg_gr_Sched = np.zeros(Simperioddays)
    E3Aor2A = np.zeros(Simperioddays)
    BWCAs_le_reg = np.zeros(Simperioddays)
    All_Cells = np.zeros(Simperioddays)
    for i in range(Simperioddays):
        Iden_WCA3Areg_plus_offset[i] = WCA_Stages_df['%s'%Pre_defined_Variables.WCA3a_REG_Zone].iloc[i] + Pre_defined_Variables.WCA3Aoffset
    
    #3gage < Sched+offset?
        if WCA_Stages_df['3gage Avg'].iloc[i] < Iden_WCA3Areg_plus_offset[i]:
            TF = True
        else:
            TF = False
        thr_gag_les_sch_pl_off[i] = TF
    #FIXME (WCA-3A floor elevation = 7.5')
    #3A-28 < min?
        if WCA_Stages_df['3A-28'].iloc[i] < Pre_defined_Variables.WCA328min:
            TF = True
        else:
            TF = False
        thrA28_less_min[i] = TF
    #FIXME (LSEL=10.5'!!!)
    #3A-NW < min?
        if WCA_Stages_df['3A-NW'].iloc[i] < Pre_defined_Variables.WCA3NWmin:
            TF = True
        else:
            TF = False
        thrANW_less_min[i] = TF
    #FIXME (LSEL=11.1')
    #2A-17 < min?
        if WCA_Stages_df['2A-17'].iloc[i] < Pre_defined_Variables.WCA217min:
            TF = True
        else:
            TF = False
        twoA17_less_min[i] = TF
    #WCA-3A Stg>iden
        if WCA_Stages_df['3gage Avg'].iloc[i] >= Iden_WCA3Areg_plus_offset[i]:
            TF = True
        else:
            TF = False
        WCA_3A_Stg_gr_Iden_WCA3Areg[i] = TF
    #WCA-2A Stg>Sched
        if WCA_Stages_df['2A-17'].iloc[i] >= WCA_Stages_df['WCA-2A Schedule'].iloc[i]:
            TF = True
        else:
            TF = False
        WCA_2A_Stg_gr_Sched[i] = TF
    #Either 3A or 2A are above
        if WCA_3A_Stg_gr_Iden_WCA3Areg[i] == True or WCA_2A_Stg_gr_Sched[i] == True:
            TF = True
        else:
            TF = False
        E3Aor2A[i] = TF   
    #(both WCAs<reg)?
        if E3Aor2A[i] ==True:
            TF = False
        elif E3Aor2A[i] == False:
            TF = True
        BWCAs_le_reg[i] = TF
    #Are WCA stages too low? (not OK to do WSA)
    #Options for constraining/limiting WSA based on WCA stages
    #0 = no constraint on WSA from WCA stages
    #1 = Offset named "WCA3Aoffset" is added to the 9.5-10.5' WCA3A regulation schedule to establish a lower limit to allow WSA to Lake O.
    #WSA is allowed only if WCA-3A stage is above this lower limit. 
    #Currently LOOPS uses the WCA-3A 3-gage avg stage from the 41-yr SFWMM simulation of the LORS-08.
    #2= This option allows multiple WCA sites to limit WSA to Lake O.  Minimum stages named "WCA2A17min", "WCA3ANmin", "WCA3A28min", etc.  Minimum stage values compared to the SFWMM-simulated stages at these sites.  WSA is allowed only if WCA stages at these sites are all above their respective minimums.  
    #See worksheet "SFWMMdata4LOOPS" beginning in column BF
    WCA_Stages_df['%s+offset'%Pre_defined_Variables.WCA3a_REG_Zone] = Iden_WCA3Areg_plus_offset
    WCA_Stages_df['3gage < Sched+offset?'] = thr_gag_les_sch_pl_off
    WCA_Stages_df['3A-28 < min?'] = thrA28_less_min
    WCA_Stages_df['3A-NW < min?'] = thrANW_less_min
    WCA_Stages_df['2A-17 < min?'] = twoA17_less_min
    WCA_Stages_df['WCA-3A Stg>Iden_WCA3Areg'] = WCA_3A_Stg_gr_Iden_WCA3Areg
    WCA_Stages_df['WCA-2A Stg>Sched'] = WCA_2A_Stg_gr_Sched
    WCA_Stages_df['Either 3A or 2A are above'] = E3Aor2A
    WCA_Stages_df['(both WCAs<reg)?'] = BWCAs_le_reg
    #First Cell
    if Pre_defined_Variables.Opt_WCAlimitWSA + 1 == 1:
        TF = False
    elif Pre_defined_Variables.Opt_WCAlimitWSA + 1 == 2:
        TF = WCA_Stages_df['3gage < Sched+offset?'].iloc[0]
    else:
        if WCA_Stages_df['3A-28 < min?'].iloc[0] == True or WCA_Stages_df['3A-NW < min?'].iloc[0] == True or WCA_Stages_df['2A-17 < min?'].iloc[0] == True:
            TF = True
        else:
            TF = False
        All_Cells[0] = TF
    if WCA_Stages_df['(both WCAs<reg)?'].iloc[1] == True or (WCA_Stages_df['3A-28 < min?'].iloc[1] == True or\
                     WCA_Stages_df['3A-NW < min?'].iloc[1] == True or WCA_Stages_df['2A-17 < min?'].iloc[1] == True):
        TF = True
    else:
        TF = False
        All_Cells[1] = TF
    for i in range(2,Simperioddays):
        if WCA_Stages_df['Either 3A or 2A are above'].iloc[i] == True and WCA_Stages_df['3A-28'].iloc[i] >= Pre_defined_Variables.WCA328min and\
                     WCA_Stages_df['3A-NW'].iloc[i] >= Pre_defined_Variables.WCA3NWmin and WCA_Stages_df['2A-17'].iloc[i] >= Pre_defined_Variables.WCA217min:
                         TF = False
        else:
                         TF = True
        All_Cells[i] = TF
    WCA_Stages_df['Are WCA stages too low?'] = All_Cells
    return(WCA_Stages_df)