def RF_Cls(Wkly_Trib_Cond_NetRF: float) -> int:
    """
    Returns a classification for the weekly tributary condition based on the
    NetRF value.

    Params:
        Wkly_Trib_Cond_NetRF(float): The NetRF value from the weekly tributary condition data.

    Returns:
        int: The classification of the weekly tributary condition.
    """
    if Wkly_Trib_Cond_NetRF >= 8:
        classification = 6
    elif 4 <= Wkly_Trib_Cond_NetRF < 8:
        classification = 5
    elif 2 <= Wkly_Trib_Cond_NetRF < 4:
        classification = 4
    elif -1 <= Wkly_Trib_Cond_NetRF < 2:
        classification = 3
    elif -3 <= Wkly_Trib_Cond_NetRF < -1:
        classification = 2
    elif -10 <= Wkly_Trib_Cond_NetRF < -3:
        classification = 1
    elif Wkly_Trib_Cond_NetRF == -999999:
        classification = 3
    else:
        classification = 0
    return classification


def MainTrib_Cls(Wkly_Trib_Cond_S65E: float) -> int:
    """
    Classifies the main tributary condition based on the S65E value.

    Args:
        Wkly_Trib_Cond_S65E (float): The S65E value from the weekly tributary condition data.

    Returns:
        int: The classification of the main tributary condition.
    """
    if Wkly_Trib_Cond_S65E >= 9000:
        classification = 6
    elif 6000 <= Wkly_Trib_Cond_S65E < 9000:
        classification = 5
    elif 3500 <= Wkly_Trib_Cond_S65E < 6000:
        classification = 4
    elif 500 <= Wkly_Trib_Cond_S65E < 3500:
        classification = 3
    elif 200 <= Wkly_Trib_Cond_S65E < 500:
        classification = 2
    elif 0 <= Wkly_Trib_Cond_S65E < 200:
        classification = 1
    elif Wkly_Trib_Cond_S65E == -999999:
        classification = 1
    else:
        classification = 0
    return classification


def Palmer_Cls(Wkly_Trib_Cond_Palmer: float) -> int:
    """
    Classify the Palmer value

    Args:
        Wkly_Trib_Cond_Palmer (float): The weekly tributary condition Palmer value.

    Returns:
        int: The classification of the Palmer value
    """
    if Wkly_Trib_Cond_Palmer >= 4:
        classification = 6
    elif 2.9999 <= Wkly_Trib_Cond_Palmer < 4:
        classification = 5
    elif 1.5 <= Wkly_Trib_Cond_Palmer < 2.9999:
        classification = 4
    elif -1.4999 <= Wkly_Trib_Cond_Palmer < 1.5:
        classification = 3
    elif -3 <= Wkly_Trib_Cond_Palmer < -1.4999:
        classification = 2
    elif -5 <= Wkly_Trib_Cond_Palmer < -3:
        classification = 1
    elif Wkly_Trib_Cond_Palmer == -999999:
        classification = 1
    else:
        classification = 0
    return classification


def NetInflow_Cls(Wkly_Trib_Cond_NetInf: float) -> int:
    """
    Classifies the weekly tributary condition based on the NetInflow value.

    Args:
        Wkly_Trib_Cond_NetInf (float): The NetInflow value from the weekly tributary condition data.

    Returns:
        int: The classification of the weekly tributary condition.
    """
    if Wkly_Trib_Cond_NetInf >= 15000:
        classification = 6
    elif 6000 <= Wkly_Trib_Cond_NetInf < 15000:
        classification = 5
    elif 2500 <= Wkly_Trib_Cond_NetInf < 6000:
        classification = 4
    elif 500 <= Wkly_Trib_Cond_NetInf < 2500:
        classification = 3
    elif -5000 <= Wkly_Trib_Cond_NetInf < 500:
        classification = 2
    elif -10000 <= Wkly_Trib_Cond_NetInf < -5000:
        classification = 1
    elif Wkly_Trib_Cond_NetInf == -999999:
        classification = 1
    else:
        classification = 0
    return classification


def LONINO_Seas_cls(LONINO_df_LONINO_Seas: float) -> int:
    """
    Classify the LONINO seasonal class.

    Args:
        LONINO_df_LONINO_Seas (float): The LONINO seasonal value.

    Returns:
        int: The classification of the LONINO seasonal class.
    """
    if LONINO_df_LONINO_Seas >= 2.0001:
        classification = 4
    elif 1.8 <= LONINO_df_LONINO_Seas < 2.0001:
        classification = 3
    elif 1.39 <= LONINO_df_LONINO_Seas < 1.8:
        classification = 2
    elif -10 <= LONINO_df_LONINO_Seas < 1.39:
        classification = 1
    elif LONINO_df_LONINO_Seas == -999999:
        classification = 1
    else:
        classification = 0
    return classification


def LONINO_M_Seas_cls(LONINO_df_LONINO_Mult_Seas: float) -> int:
    """
    Classify the LONINO multi-seasonal class.

    Args:
        LONINO_df_LONINO_Mult_Seas (float): The LONINO multi-seasonal value.

    Returns:
        int: The classification of the LONINO multi-seasonal class.
    """
    if LONINO_df_LONINO_Mult_Seas >= 2.0001:
        classification = 4
    elif 1.18 <= LONINO_df_LONINO_Mult_Seas < 2.0001:
        classification = 3
    elif 0.5001 <= LONINO_df_LONINO_Mult_Seas < 1.18:
        classification = 2
    elif -10 <= LONINO_df_LONINO_Mult_Seas < 0.5001:
        classification = 1
    elif LONINO_df_LONINO_Mult_Seas == -999999:
        classification = 1
    else:
        classification = 0
    return classification
