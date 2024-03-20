def Zone_D_Trib(Tributary_Condition, Pre_defined_Variables_Opt_NewTree):
    if Pre_defined_Variables_Opt_NewTree == 1:
        if Tributary_Condition == 1 or Tributary_Condition == 2:
            Z = 1
        elif Tributary_Condition == 3:
            Z = 2
        elif Tributary_Condition == 4:
            Z = 3
        else:
            Z = 4
    else:
        if Tributary_Condition == 1 or Tributary_Condition == 2:
            Z = 1
        elif Tributary_Condition == 3:
            Z = 2
        elif Tributary_Condition == 4 or Tributary_Condition == 5:
            Z = 3
        else:
            Z = 4
    return Z


def Zone_D_stage(Stage_LO, WSMs_Cb):
    if Stage_LO > WSMs_Cb:
        ZDS = 2
    else:
        ZDS = 1
    return ZDS


def Zone_D_Seas(
    LONINO_Seasonal_Classes, Zone_D_Trib, Pre_defined_Variables_Opt_NewTree
):
    if Pre_defined_Variables_Opt_NewTree == 0:
        if LONINO_Seasonal_Classes == 4:
            ZDSeas = 2
        else:
            ZDSeas = 1
    else:
        if Zone_D_Trib == 2:
            if LONINO_Seasonal_Classes == 1:
                ZDSeas = 1
            else:
                ZDSeas = 2
        else:
            if LONINO_Seasonal_Classes == 4:
                ZDSeas = 2
            else:
                ZDSeas = 1
    return ZDSeas


def Zone_D_MSeas(LONINO_MultiSeasonal_Classes):
    if LONINO_MultiSeasonal_Classes == 1 or LONINO_MultiSeasonal_Classes == 2:
        ZDMSeas = 1
    else:
        ZDMSeas = 2
    return ZDMSeas


def Zone_D_Rel_Code(Zone_D_Branch_Code, Pre_defined_Variables_Opt_DecTree):
    if Pre_defined_Variables_Opt_DecTree == 1:
        if Zone_D_Branch_Code in (
            1111,
            1112,
            1121,
            1122,
            1211,
            1212,
            1221,
            1222,
            2111,
            2112,
            2121,
            2211,
            2212,
            2221,
            3111,
            3211,
            3121,
            3221,
            4111,
            4121,
        ):
            ZDRC = -1
        elif Zone_D_Branch_Code in (
            2122,
            2222,
            3112,
            3122,
            3212,
            3222,
            4112,
            4122,
            4211,
            4212,
        ):
            ZDRC = 1
        elif Zone_D_Branch_Code in (4221, 4222):
            ZDRC = 4
        else:
            ZDRC = float("nan")
    else:
        ZDRC = 1
    return ZDRC


def Zone_C_Trib(Tributary_Condition, Pre_defined_Variables_Opt_NewTree):
    if Pre_defined_Variables_Opt_NewTree == 1:
        if Tributary_Condition == 1 or Tributary_Condition == 2:
            Z = 1
        elif Tributary_Condition == 3 or Tributary_Condition == 4:
            Z = 2
        else:
            Z = 3
    else:
        if (
            Tributary_Condition == 1
            or Tributary_Condition == 2
            or Tributary_Condition == 3
        ):
            Z = 1
        elif Tributary_Condition == 4 or Tributary_Condition == 5:
            Z = 2
        else:
            Z = 3
    return Z


def Zone_C_Seas(LONINO_Seasonal_Classes, Pre_defined_Variables_Opt_NewTree):
    if Pre_defined_Variables_Opt_NewTree == 1:
        if LONINO_Seasonal_Classes == 1:
            ZCSeas = 1
        else:
            ZCSeas = 2
    else:
        if LONINO_Seasonal_Classes == 1 or LONINO_Seasonal_Classes == 2:
            ZCSeas = 1
        else:
            ZCSeas = 2
    return ZCSeas


def Zone_C_MSeas(LONINO_MultiSeasonal_Classes):
    if LONINO_MultiSeasonal_Classes == 1:
        ZCMSeas = 1
    else:
        ZCMSeas = 2
    return ZCMSeas


def Zone_C_MetFcast(
    Zone_C_Seas,
    LONINO_Seasonal_Classes,
    Pre_defined_Variables_Zone_C_MetFcast_Indicator,
):
    if Pre_defined_Variables_Zone_C_MetFcast_Indicator == 0:
        MF = Zone_C_Seas
    else:
        if LONINO_Seasonal_Classes == 1 or LONINO_Seasonal_Classes == 2:
            MF = 1
        else:
            MF = 2
    return MF


def Zone_C_Rel_Code(Zone_C_Branch_Code, Pre_defined_Variables_Opt_DecTree):
    if Pre_defined_Variables_Opt_DecTree == 1:
        if Zone_C_Branch_Code in (1111, 1211):
            ZCRC = -1
        elif Zone_C_Branch_Code in (1112, 1212):
            ZCRC = 3
        elif Zone_C_Branch_Code in (
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
            ZCRC = 4
        elif Zone_C_Branch_Code in (3211, 3212, 3221, 3222):
            ZCRC = 5
        else:
            ZCRC = float("nan")
    else:
        ZCRC = 4
    return ZCRC


def Zone_B_Trib(Tributary_Condition, Pre_defined_Variables_Opt_NewTree):
    if Pre_defined_Variables_Opt_NewTree == 1:
        if Tributary_Condition == 1 or Tributary_Condition == 2:
            Z = 1
        elif Tributary_Condition == 3 or Tributary_Condition == 4:
            Z = 2
        else:
            Z = 3
    else:
        if Tributary_Condition == 1 or Tributary_Condition == 2:
            Z = 1
        elif Tributary_Condition == 6:
            Z = 3
        else:
            Z = 2
    return Z


def Zone_B_Stage(Stage_LO, Seasons_Season):
    if Stage_LO < 17.5:
        if Seasons_Season in (1, 2, 3, 4):
            St = 1
        else:
            St = 2
    else:
        St = 2
    return St


def Zone_B_Seas(LONINO_Seasonal_Classes):
    if LONINO_Seasonal_Classes == 1 or LONINO_Seasonal_Classes == 2:
        ZBSeas = 1
    else:
        ZBSeas = 2
    return ZBSeas


def Zone_B_Rel_Code(Zone_B_Branch_Code, Pre_defined_Variables_Opt_DecTree):
    if Pre_defined_Variables_Opt_DecTree == 1:
        if Zone_B_Branch_Code in (1111, 1211):
            ZBRC = 3
        elif Zone_B_Branch_Code in (
            1112,
            1121,
            1122,
            1212,
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
            3211,
            3212,
            3221,
            3222,
            1131,
            1132,
            1141,
            1142,
            1231,
            1232,
            1241,
            1242,
            2131,
            2132,
            2141,
            2142,
            2231,
            2232,
            2241,
            2242,
        ):
            ZBRC = 5
        elif Zone_B_Branch_Code in (3131, 3132, 3141, 3142, 3231, 3232, 3241, 3242):
            ZBRC = 6
        else:
            ZBRC = float("nan")
    else:
        ZBRC = 5
    return ZBRC
