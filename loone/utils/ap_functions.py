def thc_class_normal_or_above(tributary_hydrologic_condition, thc_threshold):
    """
    Determines if the tributary hydrologic condition (THC) is normal or above.

    Args:
        tributary_hydrologic_condition (float): The tributary hydrologic condition to check.
        thc_threshold (float): The threshold for determining if the condition is normal or above.

    Returns:
        bool: True if the tributary hydrologic condition is normal or above, False otherwise.
    """
    return tributary_hydrologic_condition > thc_threshold


def LStgCorres(
    month,
    WSMs_WSM1,
    Targ_Stg_df_Pre_defined_Variables_LowChance,
    Pre_defined_Variables_Opt_LChance_line,
    Pre_defined_Variables_LowChance,
):
    """
    Calculates the lake stage that is to be used with LowChance_Check.

    Args:
        month (int): The current month of the year (1-12).
        WSMs_WSM1 (float): A parameter from a water surface model.
        Targ_Stg_df_Pre_defined_Variables_LowChance (float): The target stage value under low chance conditions.
        Pre_defined_Variables_Opt_LChance_line (int): The option value being used for Opt_LChance_line (in pre_defined_variables.py).
        Pre_defined_Variables_LowChance (int): A predefined variable related to the low chance conditions.

    Returns:
        float: The calculated corresponding lake stage.
    """
    if Pre_defined_Variables_Opt_LChance_line == 1:
        if month > 5 and month < 10:
            St = WSMs_WSM1
        else:
            St = WSMs_WSM1 + 0.5
    elif Pre_defined_Variables_LowChance == 1 or (month > 5 and month < 10):
        St = -999
    else:
        St = Targ_Stg_df_Pre_defined_Variables_LowChance
    return St


def LowChance_Check(lake_o_stage_ap, lstg_corres):
    """
    Checks if the stage (water level) of Lake Okeechobee is greater than a given corresponding stage.
    Used to determine if there is a low chance of the stage falling below a certain level on june 1st.

    Args:
        Lake_O_Stage_AP (float): The stage (water level) of Lake Okeechobee.
        LStgCorres (float): The corresponding stage (water level) to compare with.

    Returns:
        bool: True if the Lake Okeechobee stage (water level) is greater than the corresponding stage, False otherwise.
    """
    return lake_o_stage_ap > lstg_corres


def Forecast_D_Sal(Pre_defined_Variables_OptSalFcast):
    """
    Calculates the forecasted daily salinity at a specific location based on the provided option?

    Args:
        OptSalFcast (int): The option being used for salinity forecast.

    Returns:
        float: The forecasted daily salinity at the specified location.
    """
    import numpy as np

    if Pre_defined_Variables_OptSalFcast == 1:
        pass
    elif Pre_defined_Variables_OptSalFcast == 2:
        pass
    elif Pre_defined_Variables_OptSalFcast == 3:
        F = np.nan
    return F


def n30d_mavg(Pre_defined_Variables_OptSalFcast):
    """
    Calculates the 30-day moving average salinity (psu).

    Args:
        Pre_defined_Variables_OptSalFcast (int): Option for salinity forecast.

    Returns:
        float: The 30-day moving average salinity.
    """
    import numpy as np

    if Pre_defined_Variables_OptSalFcast == 1:
        pass
    elif Pre_defined_Variables_OptSalFcast == 2:
        pass
    elif Pre_defined_Variables_OptSalFcast == 3:
        F = np.nan
    return F


def n30davgForecast(
    Estuary_needs_water,
    n30d_mavg_i2iplus13,
    Pre_defined_Variables_OptSalFcast,
    Pre_defined_Variables_CE_SalThreshold,
):
    """
    Calculates whether the 30 day average forecast salinity is greater than the defined CE_SalThreshold psu within 2wks.

    Args:
        Estuary_needs_water (float): The current water needs of the Estuary.
        n30d_mavg_i2iplus13 (list): A list of two weeks worth of data from a list of 30-day moving average values.
        Pre_defined_Variables_OptSalFcast (int): The salinity forecast option being used.
        Pre_defined_Variables_CE_SalThreshold (float): The predefined salinity threshold value for Caloosahatchee Estuary.

    Returns:
        bool: True if the 30-day average forecast salinity is greater than the defined CE_SalThreshold psu within 2wks, False otherwise.
    """
    if Pre_defined_Variables_OptSalFcast == 3:
        Est = Estuary_needs_water
    elif max(n30d_mavg_i2iplus13) > Pre_defined_Variables_CE_SalThreshold:
        Est = True
    else:
        Est = False
    return Est


def LORS08_bf_rel(S77BS_AP, n30davgForecast):
    """
    Returns whether or not the LORS08 schedule suggests baseflow release.

    Args:
        S77BS_AP (float): The baseflow release from Lake Okeechobee. ?
        n30davgForecast (bool): whether the 30 day average forecast salinity is greater than the defined CE_SalThreshold psu within 2wks.

    Returns:
        bool: True if the LORS08 schedule suggests baseflow release, False otherwise.
    """
    if S77BS_AP > 0:
        return n30davgForecast

    return False


def LDS_LC6_1(
    Late_Dry_Season,
    LowChance_Check,
    Pre_defined_Variables_Late_Dry_Season_Option,
):
    """
    Calculates the value of LDS_LC6_1 based on the provided parameters. LDS_LC6_1 stands for Late Dry Season Low Chance 6/1.

    Args:
        Late_Dry_Season (bool): Whether or not it is a late dry season.
        LowChance_Check (bool): The result of the LowChance_Check function (see above).
        Pre_defined_Variables_Late_Dry_Season_Option (int): The value being used for the Late_Dry_Season_Option.

    Returns:
        bool: PLACEHOLDER.
    """
    if Pre_defined_Variables_Late_Dry_Season_Option == 1:
        if Late_Dry_Season == True:
            LD = LowChance_Check
        else:
            LD = False
    else:
        LD = LowChance_Check
    return LD


def S_O(LDS_LC6_1, S77BS_AP):
    if LDS_LC6_1 == False:
        if S77BS_AP > 0:
            S = 1
        else:
            S = 0
    else:
        S = 0
    return S


def All_4(LORS08_bf_rel, LDS_LC6_1):
    if LORS08_bf_rel == True:
        All = LDS_LC6_1
    else:
        All = False
    return All


def Sabf(Lake_O_Schedule_Zone):
    if Lake_O_Schedule_Zone > 3:
        Sa = True
    else:
        Sa = False
    return Sa


def Swbf(Lake_O_Schedule_Zone):
    if Lake_O_Schedule_Zone == 3:
        Sw = True
    else:
        Sw = False
    return Sw


def Swbu(Lake_O_Schedule_Zone):
    if Lake_O_Schedule_Zone == 2:
        Swb = True
    else:
        Swb = False
    return Swb


def All_4andStage(All_4, Sabf):
    if All_4 == True and Sabf == True:
        AAA = True
    else:
        AAA = False
    return AAA


def All_4andStagein(All_4, Sabf):
    if All_4 == True and Sabf == False:
        AAA = True
    else:
        AAA = False
    return AAA


def P_AP_BF_Stg(Lake_O_Schedule_Zone, S77BS_AP):
    if (
        Lake_O_Schedule_Zone > 3 and Lake_O_Schedule_Zone < 7
    ) and S77BS_AP > 0:
        P = True
    else:
        P = False
    return P


def Logic_test_1(
    All_4andStage, P_AP_BF_Stg, Pre_defined_Variables_Opt_NoAP_above_BF_SB
):
    if Pre_defined_Variables_Opt_NoAP_above_BF_SB == 0:
        Logic = All_4andStage
    else:
        Logic = P_AP_BF_Stg
    return Logic


def Post_Ap_Baseflow(Logic_test_1, S77BS_AP, All_4andStagein, Choose_1):
    if Logic_test_1 == True:
        PB = S77BS_AP
    elif All_4andStagein == True:
        PB = min(S77BS_AP, Choose_1)
    else:
        PB = 0
    return PB


def S77RSplusPreAPS77bsf(
    S77RS_Pre_AP_S77_Baseflow,
    Lake_O_Schedule_Zone,
    Pre_defined_Variables_Opt_CEews_LOWSM,
):
    if Pre_defined_Variables_Opt_CEews_LOWSM == 0:
        if S77RS_Pre_AP_S77_Baseflow == 0 and Lake_O_Schedule_Zone > 1:
            S77RSPAPS77bsf = True
        else:
            S77RSPAPS77bsf = False
    elif S77RS_Pre_AP_S77_Baseflow == 0:
        S77RSPAPS77bsf = True
    else:
        S77RSPAPS77bsf = False
    return S77RSPAPS77bsf


def AndEstNeedsLakeWater(n30davgForecast, S77RSplusPreAPS77bsf):
    if n30davgForecast == True and S77RSplusPreAPS77bsf == True:
        AENLW = True
    else:
        AENLW = False
    return AENLW


def AndLowChance61stagelessth11(LDS_LC6_1, AndEstNeedsLakeWater):
    if LDS_LC6_1 == True and AndEstNeedsLakeWater == True:
        ALC = True
    else:
        ALC = False
    return ALC


def ATHCnora(
    AndLowChance61stagelessth11,
    THC_Class_normal_or_above,
    Late_Dry_Season,
    Pre_defined_Variables_Opt_THCbypLateDS,
):
    if Pre_defined_Variables_Opt_THCbypLateDS == 1:
        if AndLowChance61stagelessth11 == True and (
            THC_Class_normal_or_above == True or Late_Dry_Season == True
        ):
            AT = True
        else:
            AT = False
    else:
        if (
            AndLowChance61stagelessth11 == True
            and THC_Class_normal_or_above == True
        ):
            AT = True
        else:
            AT = False
    return AT


def Choose_PAPEWS_1(
    WSM_Zone,
    Pre_defined_Variables_APCB1,
    Pre_defined_Variables_APCB2,
    Pre_defined_Variables_APCB3,
    Pre_defined_Variables_APCB4,
):
    if WSM_Zone == 1:
        C1 = Pre_defined_Variables_APCB1 / 100
    elif WSM_Zone == 2:
        C1 = Pre_defined_Variables_APCB2 / 100
    elif WSM_Zone == 3:
        C1 = Pre_defined_Variables_APCB3 / 100
    elif WSM_Zone == 4:
        C1 = Pre_defined_Variables_APCB4 / 100
    else:
        C1 = -9999
    return C1


def Choose_PAPEWS_2(
    Pre_defined_Variables_Opt_CEews_LOWSM, Pre_defined_Variables_CalEst_ews
):
    if Pre_defined_Variables_Opt_CEews_LOWSM + 1 == 1:
        C2 = 0
    elif Pre_defined_Variables_Opt_CEews_LOWSM + 1 == 2:
        C2 = Pre_defined_Variables_CalEst_ews
    elif Pre_defined_Variables_Opt_CEews_LOWSM + 1 == 3:
        C2 = 300
    return C2


def Post_AP_EWS(
    ATHCnora,
    WSM_Zone,
    Choose_PAPEWS_1,
    Choose_PAPEWS_2,
    Pre_defined_Variables_CalEst_ews,
):
    if ATHCnora == True:
        if WSM_Zone == 0:
            PAEW = Pre_defined_Variables_CalEst_ews
        else:
            PAEW = (1 - Choose_PAPEWS_1) * Choose_PAPEWS_2
    else:
        PAEW = 0
    return PAEW
