def leap_year(y):
    """
    Determines if a given year is a leap year.

    A leap year is exactly divisible by 4 except for century years (years ending with 00). 
    The century year is a leap year only if it is perfectly divisible by 400. 

    Args:
        y (int): The year to check.

    Returns:
        bool: True if the year is a leap year, False otherwise.
    """
    if y % 400 == 0:
        return True
    if y % 100 == 0:
        return False
    if y % 4 == 0:
        return True
    else:
        return False


def replicate(
    year, day_num, x, targ_stg
):  # Where x is the percentage value (i.e., 10,20,30,40,50,60)
    """
    Retrieves the value from the target stage DataFrame corresponding to the given day and percentage.
    Takes leap years into account.

    Args:
        year (int): The year that day_num is in.
        day_num (int): The day number to get the value for (ranging from 0 to 365).
        x (int): The percentage value (i.e., 10, 20, 30, 40, 50, 60).
        targ_stg (pandas.DataFrame): The target stage DataFrame.

    Returns:
        float: The value from the target stage DataFrame corresponding to the day number and the given percentage.
    """
    leap_day_val = targ_stg["%d%s" % (x, "%")].iloc[59]
    if leap_year(year):
        day_num_adj = day_num
    else:
        day_num_adj = day_num + (1 if day_num >= 60 else 0)
    day_value = (
        leap_day_val
        if day_num_adj == 60 and leap_year(year) == True
        else targ_stg["%d%s" % (x, "%")].iloc[day_num_adj - 1]
    )
    return day_value
