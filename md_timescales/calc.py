import pandas as pd


def calc_mean_square_displacement(mass_center: pd.DataFrame) -> pd.DataFrame:
    """
    :param mass_center: contains [time_ns,x,y,z] columns
    :return: data frame with [time_ns, msd] columns
    """
