# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 17:04:54 2021

@author: osama
"""

class TP_Variables:
    import pandas as pd
    Z_sed = 0.05 #m
    Per_H2O_M = 0#85 #%
    Per_H2O_S = 0#20 #%
    Per_H2O_R = 0#20  #%
    Per_H2O_P = 0#85 #%
    N_Per = 0.43
    S_Per = 0.57
    ####
    Bulk_density_M = 0.15 #g/cm3
    Bulk_density_S = 1.213 #g/cm3
    Bulk_density_R = 1.213 #g/cm3
    Bulk_density_P = 0.14 #g/cm3
    ####
    Particle_density_M = 1.2 #g/cm3
    Particle_density_S = 2.56 #g/cm3
    Particle_density_R = 2.56 #g/cm3
    Particle_density_P = 1.2 #g/cm3
    ####
    A_Mud_N = 377415128 #m2 in 1988!
    A_Mud_S = 394290227 #m2 in 1988!
    A_Sand_N = 237504380 #m2 in 1988!
    A_Sand_S = 117504905 #m2 in 1988!
    A_Rock_N = 17760274 #m2 in 1988!
    A_Rock_S = 141327951 #m2 in 1988!
    A_Peat_N = 97497728 #m2 in 1988!
    A_Peat_S = 301740272 #m2 in 1988!
    A_N = A_Mud_N + A_Sand_N + A_Rock_N + A_Peat_N
    A_S = A_Mud_S + A_Sand_S + A_Rock_S + A_Peat_S
    A_tot = A_N + A_S
    Γ_inf = 91 #(mg/kg)
    #####Monthly
    v_burial_M = 0.0000003#1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)
    v_burial_S = 0.0000003#1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)
    v_burial_R = 0.0000003#1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)
    v_burial_P = 0.0000003#1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)
    ######
    # v_diff_M = 0.00319465967020851 #(m/day)#0.0906586660119641#(m/month)
    # v_diff_S = 0.00928768403057645 #(m/day)#0.0906586660119641#(m/month)
    # v_diff_R = 0.0436856478923439 #(m/day)#0.0906586660119641#(m/month)
    # v_diff_P = 0.0091905820248605 #(m/day)#0.0906586660119641#(m/month)
    # ####
    # K_decomp_M = 0.0000128902720168035 #(1/day)#0.0030004261875638304#(1/month)
    # K_decomp_S = 0.0000654610832125583 #(1/day)#0.0030004261875638304#(1/month)
    # K_decomp_R = 0.000404228978540542 #(1/day)#0.0030004261875638304#(1/month)
    # K_decomp_P = 0.000144961790358067 #(1/day)#0.0030004261875638304#(1/month)
    # ###
    # K_des_M = 0.00106466721892536 #(1/day)#0.01903931431934496#(1/month)
    # K_des_S = 0.000106506328394249 #(1/day)#0.01903931431934496#(1/month)
    # K_des_R = 0.000140539038484399  #(1/day)#0.01903931431934496#(1/month)
    # K_des_P = 0.00281682073810924 #(1/day)#0.01903931431934496#(1/month)
    # ####
    # K_ads_M = 0.000024387119180877 #(m3/mg.day)#6.000418535050426e-05#(m3/mg.month)
    # K_ads_S = 0.000015524224161757 #(m3/mg.day)#6.000418535050426e-05#(m3/mg.month)
    # K_ads_R = 0.0000351158788363223 #(m3/mg.day)#6.000418535050426e-05#(m3/mg.month)
    # K_ads_P = 0.0000151884647904376  #(m3/mg.day)#6.000418535050426e-05#(m3/mg.month)
    ###
    # v_settle = 0.03024 #0.004 #(m/day)#0.12732776837140575 #(m/month) 
    
    # Read Calibration Outputs
    Cal_Res = pd.read_csv('C:/Work/Research/LOONE/Model To be Published/LOONE_Model/Data/nondominated_Sol_var.csv')
    Par = Cal_Res['Par']
    v_diff_M = Par[0]
    v_diff_S = Par[1]
    v_diff_R = Par[2]
    v_diff_P = Par[3]
    ####
    K_decomp_M = Par[4]
    K_decomp_S = Par[5]
    K_decomp_R = Par[6]
    K_decomp_P = Par[7]
    ###
    K_des_M = Par[8]
    K_des_S = Par[9]
    K_des_R = Par[10]
    K_des_P = Par[11]
    ####
    K_ads_M = Par[12]
    K_ads_S = Par[13]
    K_ads_R = Par[14]
    K_ads_P = Par[15]
    # v_settle = Par[16]

