from typing import List


def fit_auto_correlation(time: List[float], acorr: List[float], order: int, with_constant: bool) -> List[float]:
    """
    Fit input data with :math:`\sum_n A_n \exp(-t/\\tau_n) + const`

    :param time: time data series
    :param acorr: auto-correlation data series
    :param order: number of exponentials in fit
    :param with_constant: if false const = 0
    :return: Fit curve parameters
    """
    ...


def fit_mean_square_displacement(time: List[float], msd: List[float]) -> List[float]:
    """
    Fit input data with :math:` a  t + b`

    :param time: time data series
    :param msd: mean square displacement data series
    :return: Fit curve parameters
    """
    ...
