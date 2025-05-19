# To return storage or surface area given stage (i=0) else (i = any number) will return stage given storage or surface area.
import pandas as pd
from scipy import interpolate


def stg2sto(v: float, i: int) -> float:
    """
    Calculate the Stage_Storage relationship for interpolation.

    Args:
        v (float): The stage value.
        i (int): The type of interpolation.

    Returns:
        float: The storage given the stage value, or the stage given the storage.
    """
    stgsto_data = pd.read_csv("StgSto_data.csv")
    # NOTE: We Can use cubic interpolation instead of linear
    x = stgsto_data["Stage"]
    y = stgsto_data["Storage"]
    if i == 0:
        # return storage given stage
        return interpolate.interp1d(
            x, y, fill_value="extrapolate", kind="linear"
        )(v)
    else:
        # return stage given storage
        return interpolate.interp1d(
            y, x, fill_value="extrapolate", kind="linear"
        )(v)


# Calculate the Stage_Area relationship for interpolation!
def stg2ar(v: float, i: int) -> float:
    """
    Calculate the Stage_Surface Area relationship for interpolation.

    Args:
        v(float): The stage value.
        i(int): The type of interpolation.

    Returns:
        float: The surface area given the stage value, or the stage given the surface area.
    """
    stgar_data = pd.read_csv("Stgar_data.csv")
    # NOTE: We Can use cubic interpolation instead of linear
    x = stgar_data["Stage"]
    y = stgar_data["Surf_Area"]
    if i == 0:
        # return surface area given stage
        return interpolate.interp1d(
            x, y, fill_value="extrapolate", kind="linear"
        )(v)
    else:
        # return stage given surface area
        return interpolate.interp1d(
            y, x, fill_value="extrapolate", kind="linear"
        )(v)


# Calculate the Stage_MarshArea relationship for interpolation!
def stg2mar(v: float, i: int) -> float:
    """
    Calculate the Stage_MarshArea relationship for interpolation.

    Args:
        v(float): The stage value.
        i(int): The type of interpolation.

    Returns:
        float: The interpolated value.
    """
    stgmar_data = pd.read_csv("Stgmar_data.csv")
    # NOTE: We Can use cubic interpolation instead of linear
    x = stgmar_data["Stage"]
    y = stgmar_data["Marsh_Ar"]
    if i == 0:
        # return marsh area given stage
        return interpolate.interp1d(
            x, y, fill_value="extrapolate", kind="linear"
        )(v)
    else:
        # return stage given marsh area
        return interpolate.interp1d(
            y, x, fill_value="extrapolate", kind="linear"
        )(v)
