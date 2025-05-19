import os
from yaml import safe_load


def load_config(workspace: str) -> dict:
    """Load the configuration file.

    Args:
        workspace (str): Path to directory containing the config file and its associated input data.

    Returns:
        dict: dictionary containing the configuration parameters.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
    """
    for config_file_name in ["config.yaml", "config.yml"]:
        config_file = os.path.join(workspace, config_file_name)
        if os.path.exists(config_file):
            with open(config_file, "r") as f:
                config = safe_load(f)
            return config
    else:
        raise FileNotFoundError("Config file not found in the workspace.")


def leap_year(year: int) -> bool:
    """Determines if a given year is a leap year.

    A leap year is exactly divisible by 4 except for century years (years ending with 00).
    The century year is a leap year only if it is perfectly divisible by 400.

    Args:
        year (int): The year to check.

    Returns:
        bool: True if the year is a leap year, False otherwise.
    """
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def replicate(
    year: int, day_num: int, x: int, targ_stg: dict
) -> float:  # Where x is the percentage value (i.e., 10,20,30,40,50,60)
    """Retrieves the value from the target stage DataFrame corresponding to the given day and percentage.
    Takes leap years into account.

    Args:
        year (int): The year that day_num is in.
        day_num (int): The day number to get the value for (ranging from 0 to 365).
        x (int): The percentage value (i.e., 10, 20, 30, 40, 50, 60).
        targ_stg (pandas.DataFrame): The target stage DataFrame.

    Returns:
        float: The value from the target stage DataFrame corresponding to the day number and the given percentage.
    """
    leap_day_val = targ_stg[f"{x}%"].iloc[59]
    if leap_year(year):
        day_num_adj = day_num
    else:
        day_num_adj = day_num + (1 if day_num >= 60 else 0)
    day_value = (
        leap_day_val
        if day_num_adj == 60 and leap_year(year)
        else targ_stg[f"{x}%"].iloc[day_num_adj - 1]
    )
    return day_value
