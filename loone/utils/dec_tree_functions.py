def Zone_D_Trib(tributary_condition: int, opt_new_tree: int) -> int:
    """
    Determine the zone D tributary condition based on the tributary condition
    and the option for the new tree.

    Args:
        tributary_condition(int):  The tributary condition code.
        opt_new_tree(int): The option for the new tree (0 or 1).

    Returns:
        int:  The zone D tributary condition code (1, 2, 3, or 4).
    """
    if opt_new_tree == 1:
        if tributary_condition in (1, 2):
            zone_d_trib = 1
        elif tributary_condition == 3:
            zone_d_trib = 2
        elif tributary_condition == 4:
            zone_d_trib = 3
        else:
            zone_d_trib = 4
    else:
        if tributary_condition in (1, 2):
            zone_d_trib = 1
        elif tributary_condition == 3:
            zone_d_trib = 2
        elif tributary_condition in (4, 5):
            zone_d_trib = 3
        else:
            zone_d_trib = 4
    return zone_d_trib


def Zone_D_stage(stage_lo: float, wsm_cb: float) -> int:
    """
    Determine the zone D stage based on the current lake stage and the C break point.

    Args:
        stage_lo(float): The current lake stage.
        wsm_cb(float): The stage at which the C break point is located.

    Returns:
        int: The zone D stage code (1 or 2).
    """
    return 2 if stage_lo > wsm_cb else 1


def Zone_D_Seas(lonino_seasonal_classes: int, zone_d_trib: int, opt_new_tree: int) -> int:
    """
    Determine the zone D seasonal class.

    Args:
        lonino_seasonal_classes(int): The LONINO seasonal class.
        zone_d_trib(int): The zone D tributary class.
        opt_new_tree(int): The option for the new tree (0 or 1).

    Returns:
        int: The zone D seasonal class (1 or 2).
    """
    if opt_new_tree == 0:
        zone_d_seas = 2 if lonino_seasonal_classes == 4 else 1
    else:
        if zone_d_trib == 2:
            zone_d_seas = 1 if lonino_seasonal_classes == 1 else 2
        else:
            zone_d_seas = 2 if lonino_seasonal_classes == 4 else 1
    return zone_d_seas


def Zone_D_MSeas(lonino_multiseasonal_classes: int) -> int:
    """
    Determine the zone D multi-seasonal class.

    Args:
        lonino_multiseasonal_classes (int): The LONINO multi-seasonal class.

    Returns:
        int: The zone D multi-seasonal class (1 or 2).
    """
    zone_d_mseas = 1 if lonino_multiseasonal_classes in (1, 2) else 2
    return zone_d_mseas


def Zone_D_Rel_Code(zone_d_branch_code: int, opt_dec_tree: int) -> int:
    """
    Determine the zone D release code based on the zone D branch code and
    the option for the decision tree.

    Args:
        zone_d_branch_code (int): The zone D branch code.
        opt_dec_tree (int): The option for the decision tree (0 or 1).

    Returns:
        int: The zone D release code (-1, 1, or 4).
    """
    if opt_dec_tree == 1:
        if zone_d_branch_code in (
            1111, 1112, 1121, 1122,
            1211, 1212, 1221, 1222,
            2111, 2112, 2121,
            2211, 2212, 2221,
            3111, 3211, 3121, 3221,
            4111, 4121
        ):
            zone_d_rel_code = -1
        elif zone_d_branch_code in (
            2122, 2222, 3112, 3122, 3212, 3222, 4112, 4122, 4211, 4212
        ):
            zone_d_rel_code = 1
        elif zone_d_branch_code in (4221, 4222):
            zone_d_rel_code = 4
        else:
            zone_d_rel_code = float("nan")
    else:
        zone_d_rel_code = 1

    return zone_d_rel_code


def Zone_C_Trib(tributary_condition: int, opt_new_tree: int) -> int:
    """Determine the zone C tributary condition based on the tributary condition and the option for the new tree.

    Args:
        tributary_condition (int): The tributary condition code.
        opt_new_tree (int): The option for the new tree (0 or 1).

    Returns:
        int: The zone C tributary condition code (1, 2, or 3).
    """
    if opt_new_tree == 1:
        if tributary_condition in (1, 2):
            zone_c_trib = 1
        elif tributary_condition in (3, 4):
            zone_c_trib = 2
        else:
            zone_c_trib = 3
    else:
        if tributary_condition in (1, 2, 3):
            zone_c_trib = 1
        elif tributary_condition in (4, 5):
            zone_c_trib = 2
        else:
            zone_c_trib = 3
    return zone_c_trib


def Zone_C_Seas(lonino_seasonal_classes: int, opt_new_tree: int) -> int:
    """Determine the zone C seasonal class based on the LONINO seasonal class and the option for the new tree.

    Args:
        lonino_seasonal_classes (int): The LONINO seasonal class.
        opt_new_tree (int): The option for the new tree (0 or 1).

    Returns:
        int: The zone C seasonal class (1 or 2).
    """
    if opt_new_tree == 1:
        zone_c_seasonal_class = 1 if lonino_seasonal_classes == 1 else 2
    else:
        zone_c_seasonal_class = 1 if lonino_seasonal_classes in (1, 2) else 2
    return zone_c_seasonal_class


def Zone_C_MSeas(lonino_multiseasonal_classes: int) -> int:
    """Determine the zone C multi-seasonal class.

    Args:
        lonino_multiseasonal_classes (int): The LONINO multi-seasonal class.

    Returns:
        int: The zone C multi-seasonal class (1 or 2).
    """
    return 1 if lonino_multiseasonal_classes == 1 else 2


def Zone_C_MetFcast(
    zone_c_seas: int, lonino_seas: int, opt_met_fcast: int
) -> int:
    """Determine the zone C meteorological forecast based on the zone C seasonal class
    and the option for the meteorological forecast.

    Args:
        zone_c_seas (int): The zone C seasonal class.
        lonino_seas (int): The LONINO seasonal class.
        opt_met_fcast (int): The option for the meteorological forecast (0 or 1).

    Returns:
        int: The zone C meteorological forecast (1 or 2).
    """
    if opt_met_fcast == 0:
        met_fcast = zone_c_seas
    else:
        met_fcast = 1 if lonino_seas in (1, 2) else 2
    return met_fcast


def Zone_C_Rel_Code(zone_c_branch_code: int, opt_dec_tree: int) -> int:
    """Determine the zone C release code based on the zone C branch code and the option for the decision tree.

    Args:
        zone_c_branch_code (int): The zone C branch code.
        opt_dec_tree (int): The option for the decision tree (0 or 1).

    Returns:
        int: The zone C release code (-1, 3, 4, or 5).
    """
    if opt_dec_tree == 1:
        if zone_c_branch_code in (1111, 1211):
            zone_c_rel_code = -1
        elif zone_c_branch_code in (1112, 1212):
            zone_c_rel_code = 3
        elif zone_c_branch_code in (
            1121,
            1122,
            1221,
            1222,
            2111,
            2112,
            2121,
            2122,
            2211,
            2212,
            2221,
            2222,
            3111,
            3112,
            3121,
            3122,
        ):
            zone_c_rel_code = 4
        elif zone_c_branch_code in (3211, 3212, 3221, 3222):
            zone_c_rel_code = 5
        else:
            zone_c_rel_code = float("nan")
    else:
        zone_c_rel_code = 4
    return zone_c_rel_code


def Zone_B_Trib(tributary_condition, opt_new_tree):
    """Determine the zone B tributary condition based on the tributary condition
    and the option for the new tree.

    Args:
        tributary_condition (int): The tributary condition code.
        opt_new_tree (int): The option for the new tree (0 or 1).

    Returns:
        int: The zone B tributary condition code (1, 2, or 3).
    """
    if opt_new_tree == 1:
        if tributary_condition in (1, 2):
            zone_b_trib = 1
        elif tributary_condition in (3, 4):
            zone_b_trib = 2
        else:
            zone_b_trib = 3
    else:
        if tributary_condition in (1, 2):
            zone_b_trib = 1
        elif tributary_condition == 6:
            zone_b_trib = 3
        else:
            zone_b_trib = 2
    return zone_b_trib


def Zone_B_Stage(stage_lo: float, season: int) -> int:
    """
    Determine the zone B stage based on the current lake stage and the season.

    Args:
        stage_lo (float): The current lake stage.
        season (int): The current season (1, 2, 3, or 4).

    Returns:
        int: The zone B stage code (1 or 2).
    """
    if stage_lo < 17.5:
        if season in (1, 2, 3, 4):
            zone_b_stage = 1
        else:
            zone_b_stage = 2
    else:
        zone_b_stage = 2
    return zone_b_stage


def Zone_B_Seas(lonino_seasonal_classes: int) -> int:
    """
    Determine the zone B seasonal class based on the LONINO seasonal class.

    Args:
        lonino_seasonal_classes (int): The LONINO seasonal class.

    Returns:
        int: The zone B seasonal class (1 or 2).
    """
    if lonino_seasonal_classes in (1, 2):
        zone_b_seasonal_class = 1
    else:
        zone_b_seasonal_class = 2
    return zone_b_seasonal_class


def Zone_B_Rel_Code(branch_code: int, opt_dec_tree: int) -> int:
    """
    Determine the zone B release code based on the branch code and option for the decision tree.

    Args:
        branch_code (int): The zone B branch code.
        opt_dec_tree (int): The option for the decision tree (0 or 1).

    Returns:
        int: The zone B release code (3, 5, 6, or NaN).
    """
    if opt_dec_tree == 1:
        if branch_code in (1111, 1211):
            return 3
        elif branch_code in (
            1112, 1121, 1122, 1212, 1221, 1222,
            2111, 2112, 2121, 2122, 2211, 2212,
            2221, 2222, 3111, 3112, 3121, 3122,
            3211, 3212, 3221, 3222, 1131, 1132,
            1141, 1142, 1231, 1232, 1241, 1242,
            2131, 2132, 2141, 2142, 2231, 2232,
            2241, 2242
        ):
            return 5
        elif branch_code in (
            3131, 3132, 3141, 3142, 3231, 3232,
            3241, 3242
        ):
            return 6
        else:
            return float('nan')
    else:
        return 5
