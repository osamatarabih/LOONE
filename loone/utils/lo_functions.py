# This script holds the main functions of the LOONE Model

import numpy as np


def WSM_Zone(stage_lo: float, wsm4: float, wsm3: float, wsm2: float, wsm1: float) -> int:
    """
    Determine the Water Shortage Management (WSM) zone based on the current lake stage.

    Args:
        stage_lo(float): The current lake stage.
        wsm4(float): The stage at which WSM zone 4 is located.
        wsm3(float): The stage at which WSM zone 3 is located.
        wsm2(float): The stage at which WSM zone 2 is located.
        wsm1(float): The stage at which WSM zone 1 is located.

    Returns:
        int: The WSM zone code (0, 1, 2, 3, or 4).
    """
    if stage_lo < wsm4:
        return 4
    elif stage_lo < wsm3:
        return 3
    elif stage_lo < wsm2:
        return 2
    elif stage_lo < wsm1:
        return 1
    else:
        return 0


def Max_Supply(
    wsm_zone: int,
    losa_demand_daily: float,
    cutback_z1: float,
    cutback_z2: float,
    cutback_z3: float,
    cutback_z4: float
) -> float:
    """
    Calculate the maximum supply based on the WSM zone and cutbacks.

    Args:
        wsm_zone(int): The current WSM zone (0, 1, 2, 3, or 4).
        losa_demand_daily(float): The daily demand for the Lake Okeechobee Service Area (LOSA).
        cutback_z1(float): The cutback percentage for WSM zone 1.
        cutback_z2(float): The cutback percentage for WSM zone 2.
        cutback_z3(float): The cutback percentage for WSM zone 3.
        cutback_z4(float): The cutback percentage for WSM zone 4.

    Returns:
        float: The maximum supply based on the WSM zone and cutbacks.
    """
    if wsm_zone == 1:
        max_supply = losa_demand_daily * (1 - cutback_z1)
    elif wsm_zone == 2:
        max_supply = losa_demand_daily * (1 - cutback_z2)
    elif wsm_zone == 3:
        max_supply = losa_demand_daily * (1 - cutback_z3)
    elif wsm_zone == 4:
        max_supply = losa_demand_daily * (1 - cutback_z4)
    else:
        max_supply = losa_demand_daily
    return max_supply


def LOSA_Supply(
    wsm_zone: int, losa_demand: float, max_supply: float, losa_supply_option: int
) -> float:
    """
    Calculate the LOSA supply based on the WSM zone and cutbacks.

    Args:
        wsm_zone(int): The current WSM zone (0, 1, 2, 3, or 4).
        losa_demand(float): The daily demand for the Lake Okeechobee Service Area (LOSA).
        max_supply(float): The maximum supply based on the WSM zone and cutbacks.
        losa_supply_option(int): The option for calculating the LOSA supply (1 or 2).

    Returns:
        float: The LOSA supply.
    """
    # Option 1: Supply is limited to the max supply during water shortage
    if losa_supply_option == 1:
        if wsm_zone == 0:
            supply = losa_demand
        else:
            supply = min(losa_demand, max_supply)
    # Option 2: No supply during water shortage
    elif losa_supply_option == 2:
        supply = 0
    else:
        raise ValueError("Invalid LOSA supply option")
    return supply


#
def Zone_Code(
    stage_lo: float,
    wsms_a: float,
    wsms_b: float,
    wsms_c: float,
    wsms_d3: float,
    wsms_d2: float,
    wsms_d1: float,
    wsms_d0: float,
    wsms_wsm1: float,
) -> int:
    """
    Determine the zone code based on the current lake stage and the WSMs.

    Args:
        stage_lo (float): The current lake stage.
        wsms_a (float): The stage at which WSM zone A is located.
        wsms_b (float): The stage at which WSM zone B is located.
        wsms_c (float): The stage at which WSM zone C is located.
        wsms_d3 (float): The stage at which WSM zone D3 is located.
        wsms_d2 (float): The stage at which WSM zone D2 is located.
        wsms_d1 (float): The stage at which WSM zone D1 is located.
        wsms_d0 (float): The stage at which WSM zone D0 is located.
        wsms_wsm1 (float): The stage at which WSM zone 1 is located.

    Returns:
        int: The zone code (1, 2, 3, 4, 5, 6, 7, 8, or 9).
    """
    if stage_lo > wsms_a:
        return 9
    elif stage_lo > wsms_b:
        return 8
    elif stage_lo > wsms_c:
        return 7
    elif stage_lo > wsms_d3:
        return 6
    elif stage_lo > wsms_d2:
        return 5
    elif stage_lo > wsms_d1:
        return 4
    elif stage_lo > wsms_d0:
        return 3
    elif stage_lo < wsms_wsm1:
        return 1
    else:
        return 2


def LO_Zone(zone_code: int) -> str:
    """
    Return the LO zone as a string given the zone code.

    Args:
        zone_code (int): The zone code (1, 2, 3, 4, 5, 6, 7, 8, or 9).

    Returns:
        str: The LO zone as a string.
    """
    zone_codes = {
        1: "SSM",
        2: "E",
        3: "D0",
        4: "D1",
        5: "D2",
        6: "D3",
        7: "C",
        8: "B",
        9: "A",
    }

    return zone_codes.get(zone_code, 0)


def DecTree_Relslevel(
    zone_code: int, zone_d_rel_code: int, zone_c_rel_code: int, zone_b_rel_code: int
) -> int:
    """
    Determine the decision tree release level based on the zone code and release codes.

    Args:
        zone_code (int): The zone code (1, 2, 3, 4, 5, 6, 7, 8, or 9).
        zone_d_rel_code (int): The zone D release code (-1, 1, 4, or 6).
        zone_c_rel_code (int): The zone C release code (-1, 1, 4, or 6).
        zone_b_rel_code (int): The zone B release code (3, 5, or 6).

    Returns:
        int: The decision tree release level (0, -1, 1, 2, 3, 4, 5, or 6).
    """
    if zone_code == 1:
        rel_level = 0
    elif zone_code == 2:
        rel_level = 0
    elif zone_code == 3:
        rel_level = -1
    elif zone_code == 4:
        rel_level = zone_d_rel_code
    elif zone_code == 5:
        if zone_d_rel_code == -1:
            rel_level = -1
        else:
            rel_level = min(4, zone_d_rel_code * 2)
    elif zone_code == 6:
        if zone_d_rel_code == -1:
            rel_level = -1
        else:
            rel_level = min(4, zone_d_rel_code * 3)
    elif zone_code == 7:
        rel_level = zone_c_rel_code
    elif zone_code == 8:
        rel_level = zone_b_rel_code
    else:
        rel_level = 6
    return rel_level


def PlsDay(day_flags: str, dec_tree_relslevel: int, pls_day_switch: int) -> int:
    """Determine the number of days since the last pulse event.

    Args:
        day_flags (str): The day flags string (e.g. "P.A.Day1", "P.A.Day2", etc.).
        dec_tree_relslevel (int): The decision tree release level.
        pls_day_switch (int): The switch to turn on or off the pulse day counter.

    Returns:
        int: The number of days since the last pulse event.
    """
    pls_d = 0
    if day_flags == "P.A.Day1":
        if dec_tree_relslevel > 0 and dec_tree_relslevel < 4:
            pls_d = 1
        else:
            pls_d = 0
    else:
        if pls_day_switch == 1:
            if pls_d > 0 and pls_d < 10 and dec_tree_relslevel < 4:
                pls_d += 1
            elif dec_tree_relslevel > 0 and dec_tree_relslevel < 4:
                pls_d = 1
            else:
                pls_d = 0
        elif pls_day_switch == 0:
            if pls_d > 0 and pls_d < 10:
                pls_d += 1
            elif dec_tree_relslevel > 0 and dec_tree_relslevel < 4:
                pls_d = 1
            else:
                pls_d = 0
    return pls_d


def release_level(
    prev_rl: int,
    stage_lo: float,
    tributary_condition: int,
    pls_day: int,
    zone_code: int,
    dec_tree_relslevel: int,
    max_qstg_trigger: float,
) -> int:
    """Determine the release level based on the stage of Lake Okeechobee, the tributary condition, the pulse day,
    the zone code, and the decision tree release level.

    Args:
        prev_rl (int): The previous release level.
        stage_lo (float): The stage of Lake Okeechobee.
        tributary_condition (int): The tributary condition.
        pls_day (int): The number of days since the last pulse event.
        zone_code (int): The zone code.
        dec_tree_relslevel (int): The decision tree release level.
        max_qstg_trigger (float): The maximum discharge trigger.

    Returns:
        int: The new release level.
    """
    release_level = prev_rl
    if stage_lo > max_qstg_trigger and tributary_condition == 6:
        release_level = 6
    elif pls_day in (0, 1) or zone_code > 7:
        release_level = dec_tree_relslevel
    return release_level


def ZoneCodeminus1Code(
    Zone_Code: int,
    WSMs_WSM1: float,
    WSMs_D0: float,
    WSMs_D1: float,
    WSMs_D2: float,
    WSMs_D3: float,
    WSMs_C: float,
    WSMs_B: float,
    WSMs_A: float,
) -> float:
    """Determine the stage corresponding to the given zone code minus one.

    Args:
        Zone_Code (int): The zone code.
        WSMs_WSM1 (float): The stage at which WSM zone 1 is located.
        WSMs_D0 (float): The stage at which WSM zone D0 is located.
        WSMs_D1 (float): The stage at which WSM zone D1 is located.
        WSMs_D2 (float): The stage at which WSM zone D2 is located.
        WSMs_D3 (float): The stage at which WSM zone D3 is located.
        WSMs_C (float): The stage at which WSM zone C is located.
        WSMs_B (float): The stage at which WSM zone B is located.
        WSMs_A (float): The stage at which WSM zone A is located.

    Returns:
        float: The stage corresponding to the given zone code minus one.
    """
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
    return Sam


def ZoneCodeCode(
    Zone_Code: int,
    WSMs_WSM1: float,
    WSMs_D0: float,
    WSMs_D1: float,
    WSMs_D2: float,
    WSMs_D3: float,
    WSMs_C: float,
    WSMs_B: float,
    WSMs_A: float,
) -> float:
    """
    Determine the stage corresponding to a given zone code.

    Args:
        Zone_Code (int): The zone code (1 to 8).
        WSMs_WSM1 (float): The stage at which WSM zone 1 is located.
        WSMs_D0 (float): The stage at which WSM zone D0 is located.
        WSMs_D1 (float): The stage at which WSM zone D1 is located.
        WSMs_D2 (float): The stage at which WSM zone D2 is located.
        WSMs_D3 (float): The stage at which WSM zone D3 is located.
        WSMs_C (float): The stage at which WSM zone C is located.
        WSMs_B (float): The stage at which WSM zone B is located.
        WSMs_A (float): The stage at which WSM zone A is located.

    Returns:
        float: The stage corresponding to the given zone code.
    """
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
    return Roty


def Fraction_of_Zone_height(
    zone_code: int, stage_lo: float, prev_zone_stage: float, curr_zone_stage: float
) -> float:
    """
    Calculate the fraction of the current zone's height based on the current stage.

    Args:
        zone_code (int): The zone code (1 to 9).
        stage_lo (float): The current lake stage.
        prev_zone_stage (float): The stage at which the previous zone is located.
        curr_zone_stage (float): The stage at which the current zone is located.

    Returns:
        float: The fraction of the current zone's height (0 to 1).
    """
    if zone_code == 1 or zone_code == 9:
        fraction = 1
    else:
        fraction = (stage_lo - prev_zone_stage) / (
            curr_zone_stage - prev_zone_stage + 0.01
        )
    return fraction


def ReLevelCode_1(
    release_level: int,
    dstar_D1: float,
    dstar_D2: float,
    dstar_D3: float,
    dstar_C: float,
    dstar_B: float,
) -> float:
    """
    Determine the release level code for a given release level.

    Args:
        release_level (int): The release level (-2 to 5).
        dstar_D1 (float): The pre defined variable dstar D1.
        dstar_D2 (float): The pre defined variable dstar D2.
        dstar_D3 (float): The pre defined variable dstar D3.
        dstar_C (float): The pre defined variable dstar C.
        dstar_B (float): Thepre defined variable dstar B.

    Returns:
        int: The Hnur value for the given release level.
    """
    Hnur = 0
    if (release_level + 2) == 1:
        Hnur = -99
    elif (release_level + 2) == 2:
        Hnur = -99
    elif (release_level + 2) == 3:
        Hnur = dstar_D1
    elif (release_level + 2) == 4:
        Hnur = dstar_D2
    elif (release_level + 2) == 5:
        Hnur = dstar_D3
    elif (release_level + 2) == 6:
        Hnur = dstar_C
    elif (release_level + 2) == 7:
        Hnur = dstar_B
    elif (release_level + 2) == 8:
        Hnur = -99
    return Hnur


def ReLevelCode_2(
    release_level: int,
    astar_D1: float,
    astar_D2: float,
    astar_D3: float,
    astar_C: float,
    astar_B: float,
) -> float:
    """
    Determine the release level code for a given release level.

    Args:
        release_level (int): The release level (-2 to 5).
        astar_D1 (float): The pre defined variable astar D1.
        astar_D2 (float): The pre defined variable astar D2.
        astar_D3 (float): The pre defined variable astar D3.
        astar_C (float): The pre defined variable astar C.
        astar_B (float): Thepre defined variable astar B.

    Returns:
        float: The Hnur value for the given release level.
    """
    Hnur = 0
    if (release_level + 2) == 1:
        Hnur = 0
    elif (release_level + 2) == 2:
        Hnur = 0
    elif (release_level + 2) == 3:
        Hnur = astar_D1
    elif (release_level + 2) == 4:
        Hnur = astar_D2
    elif (release_level + 2) == 5:
        Hnur = astar_D3
    elif (release_level + 2) == 6:
        Hnur = astar_C
    elif (release_level + 2) == 7:
        Hnur = astar_B
    elif (release_level + 2) == 8:
        Hnur = 0
    return Hnur


def ReLevelCode_3_S80(
    release_level: int,
    bstar_S80_D1: float,
    bstar_S80_D2: float,
    bstar_S80_D3: float,
    bstar_S80_C: float,
    bstar_S80_B: float,
) -> float:
    """
    Determine the release level code for a given release level.

    Args:
        release_level (int): The release level (-2 to 5).
        bstar_S80_D1 (float): The pre defined variable bstar S80 D1.
        bstar_S80_D2 (float): The pre defined variable bstar S80 D2.
        bstar_S80_D3 (float): The pre defined variable bstar S80 D3.
        bstar_S80_C (float): The pre defined variable bstar S80 C.
        bstar_S80_B (float): The pre defined variable bstar S80 B.

    Returns:
        float: The Hnur value for the given release level.
    """
    Hnur = 0
    if (release_level + 2) == 1:
        Hnur = 1
    elif (release_level + 2) == 2:
        Hnur = 2
    elif (release_level + 2) == 3:
        Hnur = bstar_S80_D1
    elif (release_level + 2) == 4:
        Hnur = bstar_S80_D2
    elif (release_level + 2) == 5:
        Hnur = bstar_S80_D3
    elif (release_level + 2) == 6:
        Hnur = bstar_S80_C
    elif (release_level + 2) == 7:
        Hnur = bstar_S80_B
    elif (release_level + 2) == 8:
        Hnur = 1
    return Hnur


def Outlet2DS_Mult(
    season: int,
    month: int,
    seven_day_dh: float,
    level_code_1: float,
    zone_height_fraction: float,
    level_code_2: float,
    level_code_3_s80: float,
    qreg_mult_option: int,
) -> float:
    """
    Calculate the multiplier for outlet 2DS based on the given parameters.

    Args:
        season (int): The season.
        month (int): The month.
        seven_day_dh (float): The seven day discharge height.
        level_code_1 (float): The release level code 1.
        zone_height_fraction (float): The fraction of zone height.
        level_code_2 (float): The release level code 1.
        level_code_3_s80 (float): The release level code 1 S80.
        qreg_mult_option (int): The predefined variable for Qreg multiplier option.

    Returns:
        float: The calculated multiplier based on the input conditions.
    """
    multiplier = 0
    if (
        qreg_mult_option == 0
        or (qreg_mult_option == 1 and season > 2)
        or (qreg_mult_option == 2 and 4 < month < 9)
    ):
        multiplier = 1
    elif seven_day_dh < level_code_1 and zone_height_fraction < level_code_2:
        multiplier = level_code_3_s80
    return multiplier


def Outlet2DS_Mult_2(
    month: int,
    day: int,
    pulse_day: int,
    s80_mult_i_pulse_day: float,
    s80_mult_i: float,
    pre_defined_variables_opt_qreg_mult: int,
) -> float:
    """
    Calculate the multiplier for outlet 2DS based on the given parameters.

    Args:
        month (int): The month.
        day (int): The day of the month.
        pulse_day (int): The pulse day.
        s80_mult_i_pulse_day (float): The pre-defined variable S80 mult i pulse day.
        s80_mult_i (float): The pre-defined variable S80 mult i.
        pre_defined_variables_opt_qreg_mult (int): The pre-defined variable option
            for Qreg multiplier.

    Returns:
        float: The calculated multiplier based on the input conditions.
    """
    s80_mult_2 = 0
    if (
        pre_defined_variables_opt_qreg_mult == 1
        and (month == 6 or month == 11)
        and day < 10
        and 11 > pulse_day > 1
    ):
        s80_mult_2 = s80_mult_i_pulse_day
    else:
        s80_mult_2 = s80_mult_i
    return s80_mult_2


def Outlet2DSRS(
    release_level: int,
    s80_reg_rel_rates_zone_d1: float,
    s80_avg_l1: float,
    pulses_s_80_l1: float,
    s80_mult_2: float,
    ce_sle_turns_sle_turn: float,
    s80_reg_rel_rates_zone_d2: float,
    s80_avg_l2: float,
    pulses_s_80_l2: float,
    s80_reg_rel_rates_zone_d3: float,
    s80_avg_l3: float,
    pulses_s_80_l3: float,
    s80_reg_rel_rates_zone_c: float,
    s80_reg_rel_rates_zone_b: float,
    s80_reg_rel_rates_zone_a: float,
) -> float:
    """
    Calculate the outlet 2d srs.

    Args:
        release_level (int): The release level.
        s80_reg_rel_rates_zone_d1 (float): The release rate for zone D1.
        s80_avg_l1 (float): The average lake level for zone D1.
        pulses_s_80_l1 (float): The pulses for zone D1.
        s80_mult_2 (float): The multiplier for zone D2.
        ce_sle_turns_sle_turn (float): The CE SLE turns SLE turn value.
        s80_reg_rel_rates_zone_d2 (float): The release rate for zone D2.
        s80_avg_l2 (float): The average lake level for zone D2.
        pulses_s_80_l2 (float): The pulses for zone D2.
        s80_reg_rel_rates_zone_d3 (float): The release rate for zone D3.
        s80_avg_l3 (float): The average lake level for zone D3.
        pulses_s_80_l3 (float): The pulses for zone D3.
        s80_reg_rel_rates_zone_c (float): The release rate for zone C.
        s80_reg_rel_rates_zone_b (float): The release rate for zone B.
        s80_reg_rel_rates_zone_a (float): The release rate for zone A.

    Returns:
        float: The calculated release rate.
    """
    S = 0
    if (release_level + 2) == 1:
        S = 0
    elif (release_level + 2) == 2:
        S = 0
    elif (release_level + 2) == 3:
        S = (
            s80_reg_rel_rates_zone_d1
            / s80_avg_l1
            * pulses_s_80_l1
            * s80_mult_2
            * ce_sle_turns_sle_turn
        )
    elif (release_level + 2) == 4:
        S = (
            s80_reg_rel_rates_zone_d2
            / s80_avg_l2
            * pulses_s_80_l2
            * s80_mult_2
            * ce_sle_turns_sle_turn
        )
    elif (release_level + 2) == 5:
        S = (
            s80_reg_rel_rates_zone_d3
            / s80_avg_l3
            * pulses_s_80_l3
            * s80_mult_2
            * ce_sle_turns_sle_turn
        )
    elif (release_level + 2) == 6:
        S = s80_reg_rel_rates_zone_c * s80_mult_2 * ce_sle_turns_sle_turn
    elif (release_level + 2) == 7:
        S = s80_reg_rel_rates_zone_b * s80_mult_2 * ce_sle_turns_sle_turn
    elif (release_level + 2) == 8:
        S = s80_reg_rel_rates_zone_a
    return S


def Outlet_Rel_Sim(release_level, release_rate):
    """
    Simulate the outlet release based on the release level and rates.

    Args:
        release_level (int): The release level.
        release_rate (float): The release rate chosen by the user.

    Returns:
        float: The simulated outlet release.
    """
    s = 0
    if (release_level + 2) == 1:
        s = 0
    elif (release_level + 2) == 2:
        s = 0
    else:
        s = release_rate
    return s


def Sum_Outlet2USRG1(day: int, S308RG1: float):
    """
    Calculate the cumulative sum of the S308 RG1 releases.

    Args:
        day (int): The day number.
        S308RG1 (float): The S308 RG1 release rate.

    Returns:
        float: The cumulative sum of the S308 RG1 releases.
    """
    sigma = 0
    if day == 1:
        sigma = S308RG1
    else:
        sigma = S308RG1 + sigma
    return sigma


def Outlet2DSBS(
    release_level: int,
    sum_s308rg1: float,
    vlookup1_c: float,
    stl_est_baseflow: float,
    option_s80_baseflow: int,
) -> float:
    """
    Calculate the Outlet2DSBS value from the given arguments.

    Args:
        release_level (int): The release level.
        sum_s308rg1 (float): The cumulative sum of S308 RG1 releases.
        vlookup1_c (float): The value from the VLOOKUP table.
        stl_est_baseflow (float): The estimated baseflow from the STL model.
        option_s80_baseflow (int): The option for the S80 baseflow.

    Returns:
        float: The Outlet2DSBS value.
    """
    baseflow = 0
    if release_level == -1:
        if option_s80_baseflow == 0:
            baseflow = stl_est_baseflow
        elif sum_s308rg1 / 31 <= vlookup1_c:
            baseflow = vlookup1_c
    return baseflow


def Outlet2USBK(
    lake_o_stage: float,
    wsm_d1: float,
    s308_release: float,
    c44ro: float,
    s308_backflow_data: float,
    s308_option: int,
    s308_constant: int,
    s308_threshold: float,
) -> float:
    """
    Calculate the S308 backflow rate.

    Args:
        lake_o_stage (float): The current stage of Lake Okeechobee.
        wsm_d1 (float): The water surface elevation of the Water Supply Management
            Zone 1.
        s308_release (float): The release rate from S308.
        c44ro (float): The release rate from C-44.
        s308_backflow_data (float): The backflow rate from S308.
        s308_option (int): The option for calculating the backflow rate from S308.
        s308_constant (int): The constant for calculating the backflow rate from S308.
        s308_threshold (float): The threshold stage for calculating the backflow rate
            from S308.

    Returns:
        float: The S308 backflow rate.
    """
    s308 = 0.0
    if s308_option == 1:
        if s308_constant == 1:
            if (
                lake_o_stage < min(s308_threshold, wsm_d1 - 0.25)
                and s308_release == 0
            ):
                s308 = c44ro
        else:
            s308 = s308_backflow_data
    return s308


def Outlet2USBS(
    s80bs: float, s308rg1: float, ro_east: float, option_s80_baseflow: int
) -> float:
    """
    Calculate the S308 baseflow rate.

    Args:
        s80bs (float): The S80 baseflow rate.
        s308rg1 (float): The release rate from S308.
        ro_east (float): The release rate from RO East.
        option_s80_baseflow (int): The option for calculating the S308 baseflow rate.

    Returns:
        float: The S308 baseflow rate.
    """
    s308 = 0
    if option_s80_baseflow == 0:
        s308 = max(0, s80bs - (s308rg1 + ro_east))
    else:
        s308 = max(0, s80bs - s308rg1)
    return s308


def Sum_Outlet2USBK(day: int, S308BK: float) -> float:
    """
    Calculate the sum of the S308 BK releases.

    Args:
        day (int): The day number.
        S308BK (float): The S308 BK release rate.

    Returns:
        float: The cumulative sum of the S308 BK releases.
    """
    SBK = 0
    if day == 1:
        SBK = S308BK
    else:
        SBK = S308BK + SBK
    return SBK


def Outlet2USRG_Code(
    S308RG1: float,
    S308BS: float,
    S308RG_data: float,
    STEST_data: float,
    Pre_defined_Variables_Option_RegS77S308: int,
) -> float:
    """
    Calculate the Outlet2USRG code based on the predefined option.

    Args:
        S308RG1 (float): The first release group value for S308.
        S308BS (float): The baseflow value for S308.
        S308RG_data (float): The release group data for S308.
        STEST_data (float): The test data value.
        Pre_defined_Variables_Option_RegS77S308 (int): The predefined option for regulating S77 and S308.

    Returns:
        float: The calculated Outlet2USRG code based on the option selected.
    """
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


def Outlet2USRG(
    s308rg_code: float,
    s308rg_data: float,
    stest_data: float,
    opt_s308: int,
    s308rg_const: int,
) -> float:
    """
    Calculate the Outlet2USRG value based on the predefined option.

    Args:
        s308rg_code (float): The S308 release group code.
        s308rg_data (float): The S308 release group data.
        stest_data (float): The test data value.
        opt_s308 (int): The predefined option for S308.
        s308rg_const (int): The predefined constant for S308 release group.

    Returns:
        float: The calculated Outlet2USRG value based on the option selected.
    """
    s_1 = 0
    if opt_s308 == 1:
        if s308rg_const == 1:
            s_1 = s308rg_code
        else:
            s_1 = s308rg_data + stest_data
    return s_1


def S80(
    ROeast: float,
    S308RG: float,
    S80_data: float,
    Pre_defined_Variables_S80_Const: int
) -> float:
    """
    Calculate the S80 value based on predefined constants and input parameters.

    Args:
        ROeast (float): The release rate from RO East.
        S308RG (float): The release rate from S308.
        S80_data (float): The base S80 data value.
        Pre_defined_Variables_S80_Const (int): The predefined constant option.

    Returns:
        float: The calculated S80 value.
    """
    S_80 = 0.0
    if Pre_defined_Variables_S80_Const == 1:
        S_80 = ROeast + S308RG
    else:
        S_80 = S80_data
    return S_80


def ReLevelCode_3_S77(
    Release_Level: int,
    Pre_defined_Variables_bstar_S77_D1: float,
    Pre_defined_Variables_bstar_S77_D2: float,
    Pre_defined_Variables_bstar_S77_D3: float,
    Pre_defined_Variables_bstar_S77_C: float,
    Pre_defined_Variables_bstar_S77_B: float,
) -> float:
    """
    Determine the release level code for a given release level.

    Args:
        Release_Level (int): The release level (-2 to 5).
        Pre_defined_Variables_bstar_S77_D1 (float): The pre defined variable bstar S77 D1.
        Pre_defined_Variables_bstar_S77_D2 (float): The pre defined variable bstar S77 D2.
        Pre_defined_Variables_bstar_S77_D3 (float): The pre defined variable bstar S77 D3.
        Pre_defined_Variables_bstar_S77_C (float): The pre defined variable bstar S77 C.
        Pre_defined_Variables_bstar_S77_B (float): The pre defined variable bstar S77 B.

    Returns:
        float: The Hnur value for the given release level.
    """
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
    return Hnur


def Outlet1US_Mult(
    Seasons_Season: int,
    Seasons_Month: int,
    dh_7days: float,
    ReLevelCode_1: float,
    Fraction_of_Zone_height: float,
    ReLevelCode_2: float,
    ReLevelCode_3_S77: float,
    Pre_defined_Variables_Opt_QregMult: int,
) -> float:
    """
    Calculate the multiplier for Outlet1US based on the given parameters.

    Args:
        Seasons_Season (int): The season.
        Seasons_Month (int): The month.
        dh_7days (float): The 7 day discharge height.
        ReLevelCode_1 (float): The release level code 1.
        Fraction_of_Zone_height (float): The fraction of zone height.
        ReLevelCode_2 (float): The release level code 2.
        ReLevelCode_3_S77 (float): The release level code 3 S77.
        Pre_defined_Variables_Opt_QregMult (int): The predefined variable option for Qreg multiplier.

    Returns:
        float: The calculated multiplier based on the input conditions.
    """
    Mult = 0
    if (
        Pre_defined_Variables_Opt_QregMult == 0
        or (Pre_defined_Variables_Opt_QregMult == 1 and Seasons_Season > 2)
        or (
            Pre_defined_Variables_Opt_QregMult == 2
            and Seasons_Month > 4
            and Seasons_Month < 9
        )
    ):
        Mult = 1
    elif dh_7days < ReLevelCode_1 and Fraction_of_Zone_height < ReLevelCode_2:
        Mult = ReLevelCode_3_S77
    return Mult


def Outlet1US_Mult_2(
    month: int,
    day: int,
    PlsDay: int,
    S77_Mult_i_PlsDay: float,
    S77_Mult_i: float,
    Pre_defined_Variables_Opt_QregMult: int,
) -> float:
    """
    Calculate the multiplier for Outlet1US based on the given parameters.

    Args:
        month (int): The month.
        day (int): The day of the month.
        PlsDay (int): The pulse day.
        S77_Mult_i_PlsDay (float): The S77 multiplier for the pulse day.
        S77_Mult_i (float): The S77 multiplier for the day.
        Pre_defined_Variables_Opt_QregMult (int): The predefined variable option for Qreg multiplier.

    Returns:
        float: The calculated multiplier based on the input conditions.
    """
    S77_M2 = 0
    if (
        Pre_defined_Variables_Opt_QregMult == 1
        and (month == 6 or month == 11)
        and day < 10
        and PlsDay > 1
        and PlsDay < 11
    ):
        S77_M2 = S77_Mult_i_PlsDay
    else:
        S77_M2 = S77_Mult_i
    return S77_M2


def Outlet1USRS(
    Release_Level: int,
    S77_RegRelRates_Zone_D1: float,
    S77avgL1: float,
    Pulses_S_77_L1: float,
    S77_Mult_2: float,
    C43RO: float,
    CE_SLE_turns_CEturn: float,
    S77_RegRelRates_Zone_D2: float,
    S77avgL2: float,
    Pulses_S_77_L2: float,
    Zone_Code: int,
    S77_RegRelRates_Zone_D3: float,
    S77avgL3: float,
    Pulses_S_77_L3: float,
    S77_RegRelRates_Zone_C: float,
    S77_RegRelRates_Zone_B: float,
    S77_RegRelRates_Zone_A: float,
    Pre_defined_Variables_Opt_Outlet1DSRG: float,
) -> float:
    """
    Calculate the Outlet1USRS value based on the given parameters.

    Parameters:
        Release_Level (int): The level of release.
        S77_RegRelRates_Zone_D1 (float): The S77 regulatory release rate for zone D1.
        S77avgL1 (float): The average level for zone D1.
        Pulses_S_77_L1 (float): The number of pulses for zone D1.
        S77_Mult_2 (float): The multiplier for S77.
        C43RO (float): The C43 runoff value.
        CE_SLE_turns_CEturn (float): The CE SLE turns value.
        S77_RegRelRates_Zone_D2 (float): The S77 regulatory release rate for zone D2.
        S77avgL2 (float): The average level for zone D2.
        Pulses_S_77_L2 (float): The number of pulses for zone D2.
        Zone_Code (int): The zone code.
        S77_RegRelRates_Zone_D3 (float): The S77 regulatory release rate for zone D3.
        S77avgL3 (float): The average level for zone D3.
        Pulses_S_77_L3 (float): The number of pulses for zone D3.
        S77_RegRelRates_Zone_C (float): The S77 regulatory release rate for zone C.
        S77_RegRelRates_Zone_B (float): The S77 regulatory release rate for zone B.
        S77_RegRelRates_Zone_A (float): The S77 regulatory release rate for zone A.
        Pre_defined_Variables_Opt_Outlet1DSRG (float): The predefined variable option for Outlet1DSRG.

    Returns:
        float: The calculated Outlet1USRS value based on the release level and other parameters.
    """
    S = 0
    if (Release_Level + 2) == 1:
        S = 0
    elif (Release_Level + 2) == 2:
        S = 0
    elif (Release_Level + 2) == 3:
        S = (
            max(
                0,
                S77_RegRelRates_Zone_D1
                / S77avgL1
                * Pulses_S_77_L1
                * S77_Mult_2
                - Pre_defined_Variables_Opt_Outlet1DSRG * C43RO,
            )
            * CE_SLE_turns_CEturn
        )
    elif (Release_Level + 2) == 4:
        S = (
            max(
                0,
                S77_RegRelRates_Zone_D2
                / S77avgL2
                * Pulses_S_77_L2
                * S77_Mult_2
                - Pre_defined_Variables_Opt_Outlet1DSRG * C43RO,
            )
            * CE_SLE_turns_CEturn
        )
    elif (Release_Level + 2) == 5 and Zone_Code < 7:
        S = (
            max(
                0,
                S77_RegRelRates_Zone_D3
                / S77avgL3
                * Pulses_S_77_L3
                * S77_Mult_2
                - Pre_defined_Variables_Opt_Outlet1DSRG * C43RO,
            )
            * CE_SLE_turns_CEturn
        )
    elif (Release_Level + 2) == 5 and Zone_Code >= 7:
        S = (
            max(
                0,
                S77_RegRelRates_Zone_D3
                / S77avgL3
                * Pulses_S_77_L3
                * S77_Mult_2,
            )
            * CE_SLE_turns_CEturn
        )
    elif (Release_Level + 2) == 6:
        S = S77_RegRelRates_Zone_C * S77_Mult_2 * CE_SLE_turns_CEturn
    elif (Release_Level + 2) == 7:
        S = S77_RegRelRates_Zone_B * S77_Mult_2 * CE_SLE_turns_CEturn
    elif (Release_Level + 2) == 8:
        S = S77_RegRelRates_Zone_A
    return S


def Sum_Outlet1USRS(day: int, Outlet1USRS: float) -> float:
    """
    Calculate the cumulative sum of Outlet1USRS.

    Args:
        day (int): The day number.
        Outlet1USRS (float): The Outlet1USRS value

    Returns:
        float: The cumulative sum of Outlet1USRS
    """
    Sigma = 0
    if day == 1:
        Sigma = Outlet1USRS
    else:
        Sigma = Outlet1USRS + Sigma
    return Sigma


def Outlet1USBK(
    Stage_LO: float,
    Outlet1USRS: int,
    Outlet1USBSAP: int,
    Outlet1USEWS: int,
    C43RO: float,
    S77BK_data: float,
    Pre_defined_Variables_Outlet1USBK_Switch: int,
    Pre_defined_Variables_Outlet1USBK_Threshold: float,
) -> float:
    """
    Calculate the Outlet1USBK value based on the given parameters.

    Args:
        Stage_LO (float): The stage level of Lake Okeechobee.
        Outlet1USRS (int): The Outlet1USRS value.
        Outlet1USBSAP (int): The Outlet1USBSAP value.
        Outlet1USEWS (int): The Outlet1USEWS value.
        C43RO (float): The C43 runoff value.
        S77BK_data (float): The S77BK data value.
        Pre_defined_Variables_Outlet1USBK_Switch (int): The predefined switch variable.
        Pre_defined_Variables_Outlet1USBK_Threshold (float): The predefined threshold variable.

    Returns:
        float: The calculated Outlet1USBK value.
    """
    Y: float = 0
    if Pre_defined_Variables_Outlet1USBK_Switch == 1:
        if (
            Stage_LO < Pre_defined_Variables_Outlet1USBK_Threshold
            and Outlet1USRS == 0
            and Outlet1USBSAP == 0
            and Outlet1USEWS == 0
        ):
            Y = max(0, 0.3 * C43RO - 5)
        else:
            Y = 0
    else:
        Y = S77BK_data
    return Y


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
def Outlet1DSBS(
    Release_Level: int,
    Sum_Outlet1USRS: float,
    VLOOKUP2_c: float,
    CalEst_baseflow: float,
    Pre_defined_Variables_Option_S77Baseflow: int,
) -> float:
    """
    Calculate the Outlet1DSBS value based on the given parameters.

    Args:
        Release_Level (int): The release level.
        Sum_Outlet1USRS (float): The cumulative sum of Outlet1USRS.
        VLOOKUP2_c (float): The VLOOKUP2 constant.
        CalEst_baseflow (float): The calculated estimated baseflow.
        Pre_defined_Variables_Option_S77Baseflow (int): The predefined option for S77 baseflow.

    Returns:
        float: The calculated Outlet1DSBS value.
    """
    BS = 0
    if Release_Level >= -1:
        if Pre_defined_Variables_Option_S77Baseflow == 0:
            BS = CalEst_baseflow
        elif Sum_Outlet1USRS / 31 > VLOOKUP2_c:
            BS = 0
        else:
            BS = VLOOKUP2_c
    else:
        BS = 0
    return BS


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
def Outlet1USBS(
    Outlet1DSBS: float, Outlet1USRS: float, ROwest: float, Pre_defined_Variables_Option_S77Baseflow: int
) -> float:
    """
    Calculate the Outlet1USBS value based on the given parameters.

    Args:
        Outlet1DSBS (float): The Outlet1DSBS value.
        Outlet1USRS (float): The Outlet1USRS value.
        ROwest (float): The ROwest value.
        Pre_defined_Variables_Option_S77Baseflow (int): The predefined option for S77 baseflow.

    Returns:
        float: The calculated Outlet1USBS value.
    """
    Z = 0
    if Pre_defined_Variables_Option_S77Baseflow == 0:
        Z = max(0, Outlet1DSBS - (Outlet1USRS + ROwest))
    else:
        Z = max(0, Outlet1DSBS - Outlet1USRS)
    return Z


def Outlet1USBSAP(
    Outlet1USBS: float,  # The Outlet1USBS value.
    Post_Ap_Baseflow: float,  # The Post-AP baseflow value.
    Pre_defined_Variables_Opt_AdapProt: int  # The option for AdapProt.
) -> float:
    """
    Calculate the Outlet1USBSAP value based on the given parameters.

    Args:
        Outlet1USBS (float): The Outlet1USBS value.
        Post_Ap_Baseflow (float): The Post-AP baseflow value.
        Pre_defined_Variables_Opt_AdapProt (int): The option for AdapProt.

    Returns:
        float: The calculated Outlet1USBSAP value.
    """
    S77 = 0
    if Pre_defined_Variables_Opt_AdapProt == 0:
        S77 = Outlet1USBS
    else:
        S77 = Post_Ap_Baseflow
    return S77


def Outlet1USEWS(
    Post_AP_EWS: float,  # The Post-AP EWS value.
    CAEST_data: float,  # The CAEST data value.
    Pre_defined_Variables_Outlet1USEWS_Switch: int,  # The switch variable.
    Pre_defined_Variables_Opt_AdapProt: int,  # The option for AdapProt.
) -> float:
    """
    Calculate the Outlet1USEWS value based on the given parameters.

    Args:
        Post_AP_EWS (float): The Post-AP EWS value.
        CAEST_data (float): The CAEST data value.
        Pre_defined_Variables_Outlet1USEWS_Switch (int): The switch variable.
        Pre_defined_Variables_Opt_AdapProt (int): The option for AdapProt.

    Returns:
        float: The calculated Outlet1USEWS value.
    """
    X = 0
    if Pre_defined_Variables_Outlet1USEWS_Switch == 1:
        if Pre_defined_Variables_Opt_AdapProt == 0:
            X = 0
        else:
            X = Post_AP_EWS
    else:
        X = CAEST_data
    return X


def Outlet1USREG(
    Outlet1USRS: float,  # The Outlet1USRS value.
    Outlet1USBSAP: float,  # The Outlet1USBSAP value.
    S77RG_data: float,  # The S77RG data value.
    Pre_defined_Variables_Outlet1USREG_Switch: int,  # The switch variable.
    Pre_defined_Variables_Option_RegS77S308: int,  # The option for regulating S77 and S308.
) -> float:
    """
    Calculate the Outlet1USREG value based on the given parameters.

    Args:
        Outlet1USRS (float): The Outlet1USRS value.
        Outlet1USBSAP (float): The Outlet1USBSAP value.
        S77RG_data (float): The S77RG data value.
        Pre_defined_Variables_Outlet1USREG_Switch (int): The switch variable.
        Pre_defined_Variables_Option_RegS77S308 (int): The option for regulating S77 and S308.

    Returns:
        float: The calculated Outlet1USREG value.
    """
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
    return W


def Outlet1DS(
    Outlet1USREG: float,  # The Outlet1USREG value.
    Outlet1USEWS: float,  # The Outlet1USEWS value.
    ROwest: float,  # The ROwest value.
    Outlet1DS_data: float,  # The Outlet1DS data value.
    Pre_defined_Variables_Outlet1DS_Switch: int  # The switch variable.
) -> float:
    """
    Calculate the Outlet1DS value based on the given parameters.

    Args:
        Outlet1USREG (float): The Outlet1USREG value.
        Outlet1USEWS (float): The Outlet1USEWS value.
        ROwest (float): The ROwest value.
        Outlet1DS_data (float): The Outlet1DS data value.
        Pre_defined_Variables_Outlet1DS_Switch (int): The switch variable.

    Returns:
        float: The calculated Outlet1DS value.
    """
    S = 0
    if Pre_defined_Variables_Outlet1DS_Switch == 1:
        S = Outlet1USREG + Outlet1USEWS + ROwest
    else:
        S = Outlet1DS_data
    return S


def Choose_WCA(
    RegWCA_data: float,  # The RegWCA data value.
    Pre_defined_Variables_Option_RegWCA: int,  # The option for regulating WCA.
    Pre_defined_Variables_Constant_RegWCA: float,  # The constant value for regulating WCA.
) -> float:
    """
    Calculate the chosen value based on the given parameters.

    Args:
        RegWCA_data (float): The RegWCA data value.
        Pre_defined_Variables_Option_RegWCA (int): The option for regulating WCA.
        Pre_defined_Variables_Constant_RegWCA (float): The constant value for regulating WCA.

    Returns:
        float: The chosen value.
    """
    Ch = 0
    if Pre_defined_Variables_Option_RegWCA == 1:
        Ch = 0
    elif Pre_defined_Variables_Option_RegWCA == 2:
        Ch = RegWCA_data
    elif Pre_defined_Variables_Option_RegWCA == 3:
        Ch = 0
    elif Pre_defined_Variables_Option_RegWCA == 4:
        Ch = Pre_defined_Variables_Constant_RegWCA
    return Ch


def Choose_L8C51(
    RegL8C51_data: float,  # The RegL8C51 data value.
    Pre_defined_Variables_Option_RegL8C51: int,  # The option for regulating L8C51.
    Pre_defined_Variables_Constant_RegL8C51: float,  # The constant value for regulating L8C51.
) -> float:
    """
    Calculate the chosen value based on the given parameters.

    Args:
        RegL8C51_data (float): The RegL8C51 data value.
        Pre_defined_Variables_Option_RegL8C51 (int): The option for regulating L8C51.
        Pre_defined_Variables_Constant_RegL8C51 (float): The constant value for regulating L8C51.

    Returns:
        float: The chosen value.
    """
    Ch = 0
    if Pre_defined_Variables_Option_RegL8C51 == 1:
        Ch = 0
    elif Pre_defined_Variables_Option_RegL8C51 == 2:
        Ch = RegL8C51_data
    elif Pre_defined_Variables_Option_RegL8C51 == 3:
        Ch = 0
    elif Pre_defined_Variables_Option_RegL8C51 == 4:
        Ch = Pre_defined_Variables_Constant_RegL8C51
    return Ch


def ET(
    et_dry_data: float,
    stage_to_area: float,
    et_litoral_data: float,
    stage_to_marsh: float,
    et_open_data: float,
    et_vol_data: float,
    predefined_variables_et_switch: int,
) -> float:
    """
    Calculate the evapotranspiration volume based on the provided parameters.

    Args:
        et_dry_data (float): Evapotranspiration data for dry areas.
        stage_to_area (float): Stage to area conversion value.
        et_litoral_data (float): Evapotranspiration data for litoral areas.
        stage_to_marsh (float): Stage to marsh conversion value.
        et_open_data (float): Evapotranspiration data for open areas.
        et_vol_data (float): Predefined ET volume data.
        predefined_variables_et_switch (int): Switch to determine calculation method.

    Returns:
        float: The calculated evapotranspiration volume.
    """
    evapotranspiration = 0.0
    if predefined_variables_et_switch == 1:
        evapotranspiration = (
            ((et_dry_data / 12) * (466000 - stage_to_area))
            + ((et_litoral_data / 12) * stage_to_marsh)
            + ((et_open_data / 12) * (stage_to_area - stage_to_marsh))
        )
    else:
        evapotranspiration = et_vol_data
    return evapotranspiration


def Choose_WSA_1(
    WSMs_WSM1: float,  # The WSMs WSM1 value.
    Pre_defined_Variables_Opt_WSA: int,  # The option for WSA.
    Pre_defined_Variables_WSAtrig2: float,  # The WSA trigger value.
    Pre_defined_Variables_WSAoff2: float,  # The WSA offset value.
) -> float:
    """
    Calculate the chosen WSA value based on the given parameters.

    Args:
        WSMs_WSM1 (float): The WSMs WSM1 value.
        Pre_defined_Variables_Opt_WSA (int): The option for WSA.
        Pre_defined_Variables_WSAtrig2 (float): The WSA trigger value.
        Pre_defined_Variables_WSAoff2 (float): The WSA offset value.

    Returns:
        float: The chosen WSA value.
    """
    Ch1 = 0
    if Pre_defined_Variables_Opt_WSA == 1:
        Ch1 = Pre_defined_Variables_WSAtrig2
    elif Pre_defined_Variables_Opt_WSA == 2:
        Ch1 = WSMs_WSM1 + Pre_defined_Variables_WSAoff2
    else:
        Ch1 = np.nan
    return Ch1


def Choose_WSA_2(
    WSMs_WSM1: float,  # The WSMs WSM1 value.
    Pre_defined_Variables_Opt_WSA: int,  # The option for WSA.
    Pre_defined_Variables_WSAtrig1: float,  # The WSA trigger value.
    Pre_defined_Variables_WSAoff1: float,  # The WSA offset value.
) -> float:
    """
    Calculate the chosen WSA value based on the given parameters.

    Args:
        WSMs_WSM1 (float): The WSMs WSM1 value.
        Pre_defined_Variables_Opt_WSA (int): The option for WSA.
        Pre_defined_Variables_WSAtrig1 (float): The WSA trigger value.
        Pre_defined_Variables_WSAoff1 (float): The WSA offset value.

    Returns:
        float: The chosen WSA value.
    """
    Ch2 = 0
    if Pre_defined_Variables_Opt_WSA == 1:
        Ch2 = Pre_defined_Variables_WSAtrig1
    elif Pre_defined_Variables_Opt_WSA == 2:
        Ch2 = WSMs_WSM1 + Pre_defined_Variables_WSAoff1
    else:
        Ch2 = np.nan
    return Ch2


def WSA_MIA(
    Are_WCA_stages_too_low: bool,  # A boolean indicating whether WCA stages are too low.
    LONINO_Seasonal_Classes: int,  # The LONINO seasonal classes.
    Stage_LO: float,  # The current lake stage.
    Choose_WSA_1: float,  # The chosen WSA value based on the parameters.
    MIA_cfs: float,  # The MIA value in cfs.
    S3PMP: float,  # The S3PMP value.
    Choose_WSA_2: float,  # The chosen WSA value based on the parameters.
    Pre_defined_Variables_Opt_WSA: int,  # The option for WSA.
    Pre_defined_Variables_WSA_THC: int,  # The WSA THC value.
    Pre_defined_Variables_MIAcap2: float,  # The MIA cap value 2.
    Pre_defined_Variables_MIAcap1: float,  # The MIA cap value 1.
) -> float:
    """
    Calculate the chosen WSA value based on the given parameters.

    Args:
        Are_WCA_stages_too_low (bool): A boolean indicating whether WCA stages are too low.
        LONINO_Seasonal_Classes (int): The LONINO seasonal classes.
        Stage_LO (float): The current lake stage.
        Choose_WSA_1 (float): The chosen WSA value based on the parameters.
        MIA_cfs (float): The MIA value in cfs.
        S3PMP (float): The S3PMP value.
        Choose_WSA_2 (float): The chosen WSA value based on the parameters.
        Pre_defined_Variables_Opt_WSA (int): The option for WSA.
        Pre_defined_Variables_WSA_THC (int): The WSA THC value.
        Pre_defined_Variables_MIAcap2 (float): The MIA cap value 2.
        Pre_defined_Variables_MIAcap1 (float): The MIA cap value 1.

    Returns:
        float: The chosen WSA value.
    """
    WSA = 0
    if (
        Pre_defined_Variables_Opt_WSA == 0
        or Are_WCA_stages_too_low
        or LONINO_Seasonal_Classes > Pre_defined_Variables_WSA_THC
    ):
        WSA = 0
    elif Stage_LO < Choose_WSA_1:
        WSA = max(0, min(Pre_defined_Variables_MIAcap2, MIA_cfs) - S3PMP)
    elif Stage_LO < Choose_WSA_2:
        WSA = max(0, min(Pre_defined_Variables_MIAcap1, MIA_cfs) - S3PMP)
    else:
        WSA = 0
    return WSA


def WSA_NNR(
    Are_WCA_stages_too_low: bool,  # Indicates if WCA stages are too low.
    LONINO_Seasonal_Classes: int,  # The LONINO seasonal class.
    Stage_LO: float,  # The current lake stage.
    Choose_WSA_1: float,  # The chosen WSA trigger value 1.
    NNR_cfs: float,  # The NNR value in cubic feet per second.
    S2PMP: float,  # The S2PMP value.
    Choose_WSA_2: float,  # The chosen WSA trigger value 2.
    Pre_defined_Variables_Opt_WSA: int,  # Option for WSA.
    Pre_defined_Variables_WSA_THC: int,  # WSA THC threshold.
    Pre_defined_Variables_NNRcap2: float,  # NNR cap value 2.
    Pre_defined_Variables_NNRcap1: float,  # NNR cap value 1.
) -> float:
    """
    Calculate the chosen WSA value based on the given parameters.

    Args:
        Are_WCA_stages_too_low (bool): Indicates if WCA stages are too low.
        LONINO_Seasonal_Classes (int): The LONINO seasonal class.
        Stage_LO (float): The current lake stage.
        Choose_WSA_1 (float): The chosen WSA trigger value 1.
        NNR_cfs (float): The NNR value in cubic feet per second.
        S2PMP (float): The S2PMP value.
        Choose_WSA_2 (float): The chosen WSA trigger value 2.
        Pre_defined_Variables_Opt_WSA (int): Option for WSA.
        Pre_defined_Variables_WSA_THC (int): WSA THC threshold.
        Pre_defined_Variables_NNRcap2 (float): NNR cap value 2.
        Pre_defined_Variables_NNRcap1 (float): NNR cap value 1.

    Returns:
        float: The calculated WSA value.
    """
    WSA = 0
    if (
        Pre_defined_Variables_Opt_WSA == 0
        or Are_WCA_stages_too_low
        or LONINO_Seasonal_Classes > Pre_defined_Variables_WSA_THC
    ):
        WSA = 0
    elif Stage_LO < Choose_WSA_1:
        WSA = max(0, min(Pre_defined_Variables_NNRcap2, NNR_cfs) - S2PMP)
    elif Stage_LO < Choose_WSA_2:
        WSA = max(0, min(Pre_defined_Variables_NNRcap1, NNR_cfs) - S2PMP)
    else:
        WSA = 0
    return WSA


def Storage(
    day_flags: str, storage_i: float, start_storage: float, next_day_storage: float, storage_change: float
) -> float:
    """
    Calculate the storage based on the given parameters.

    Args:
        day_flags (str): The day flags.
        storage_i (float): The current storage.
        start_storage (float): The start storage.
        next_day_storage (float): The storage for the next day.
        storage_change (float): The storage change.

    Returns:
        float: The calculated storage.
    """
    if day_flags == "CS start date":
        storage = storage_i
    elif day_flags == "SimDay1":
        storage = start_storage
    else:
        storage = next_day_storage + storage_change

    return storage


def Lake_Stage(
    stg2sto_Storage_iplus2: float, EOD_Stg_data: float, Pre_defined_Variables_Option_Stage: int
) -> float:
    """
    Calculate the lake stage based on the given parameters.

    Args:
        stg2sto_Storage_iplus2 (float): The storage for two days from now.
        EOD_Stg_data (float): The end of day stage data.
        Pre_defined_Variables_Option_Stage (int): The option for the stage calculation.

    Returns:
        float: The calculated lake stage.
    """
    Sta = 0
    if Pre_defined_Variables_Option_Stage + 1 == 1:
        Sta = stg2sto_Storage_iplus2
    elif Pre_defined_Variables_Option_Stage + 1 == 2:
        Sta = np.nan
    elif Pre_defined_Variables_Option_Stage + 1 == 3:
        Sta = EOD_Stg_data
    elif Pre_defined_Variables_Option_Stage + 1 == 4:
        Sta = np.nan
    elif Pre_defined_Variables_Option_Stage + 1 == 5:
        Sta = np.nan
    return Sta
