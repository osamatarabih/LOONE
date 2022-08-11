# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 18:47:55 2022

@author: osama
"""

class Add_Fn:
    def leap_year(y):
        if y % 400 == 0:
            return True
        if y % 100 == 0:
            return False
        if y % 4 == 0:
            return True
        else:
            return False
        
    def Replicate (year, day_num,x,Targ_Stg): #Where x is the percentage value (i.e., 10,20,30,40,50,60)
        leap_day_val = Targ_Stg['%d%s' % (x,'%')].iloc[59]
        if Add_Fn.leap_year(year) == True:
            day_num_adj = day_num
        else:
            day_num_adj = day_num + (1 if day_num >= 60 else 0)
        day_value = leap_day_val if day_num_adj == 60 and Add_Fn.leap_year(year) == True else Targ_Stg['%d%s' % (x,'%')].iloc[day_num_adj-1]
        return(day_value) 