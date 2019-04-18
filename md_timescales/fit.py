from typing import List, Tuple, Union
import numpy as np
from bionmr_utils.md import *
import os
import pandas as pd
import glob
from scipy.optimize import curve_fit


def fit_limit(data):
    def moving_average(data_set, periods=3):
        weights = np.ones(periods) / periods
        return np.convolve(data_set, weights, mode='valid')

    diff = np.diff(data)
    pos_diff = (diff > 0).astype(int)
    window_size = 50
    pos_diff_avg = moving_average(pos_diff, window_size)
    index = np.argmax(pos_diff_avg > 0.5) + window_size // 2
    return index


def __multi_exp_f(x, A, TAU, C):
    return np.sum(
        (a * np.exp(-x / tau)) for a, tau in zip(A, TAU)
    ) + C


def multi_exp(x, *args):
    TAU = args[1::2]

    if len(args) % 2 == 0:
        C = 0
        A = args[::2]
    else:
        C = args[-1]
        A = args[:-1:2]

    return __multi_exp_f(x, A, TAU, C)


def fit_auto_correlation(time: List[float], acorr: List[float], bounds) -> Tuple[int, List[float]]:
    """
    Fit input data with :math:`\sum_n A_n \exp(-t/\\tau_n) + const`

    :param time: time data series
    :param acorr: auto-correlation data series
    :param bounds: curve parameters bounds
    :return: Fit curve parameters
    """

    limit = fit_limit(acorr)
    p0 = np.mean(bounds, axis=0)
    popt, pcov = curve_fit(multi_exp,
                           time[:limit],
                           acorr[:limit],
                           p0=p0,
                           bounds=bounds)
    return limit, popt


def fit_mean_square_displacement(time: List[float], msd: List[float]) -> List[float]:
    """
    Fit input data with :math:` a  t + b`

    :param time: time data series
    :param msd: mean square displacement data series
    :return: Fit curve parameters
    """
    ...


def save_fit_auto_correlation(path_to_ref: str,
                              path_to_csv_acorr: str,
                              output_directory: str,
                              curve_bounds: List[List[List[Union[float, int]]]],
                              tumbling=False
                              ):
    traj, ref = traj_from_dir(path_to_ref, first=1, last=1)

    chain = ref.asChains[0]
    csv_files = sorted(glob.glob(os.path.join(path_to_csv_acorr, "*.csv")))

    for bounds in curve_bounds:
        with_constant = len(bounds[0]) % 2 == 1
        order = (len(bounds[0]) + 1) // 2
        tau_table = pd.DataFrame()
        for file in csv_files:

            df = pd.read_csv(file)
            limit, popt = fit_auto_correlation(df.time_ns,
                                               df.acorr,
                                               bounds=bounds)

            name = os.path.splitext(os.path.basename(file))[0]

            if tumbling:
                axis = name.split("_")[-1]
                D = {
                    'axis': axis, 'limit': limit
                }
            else:
                rid, aname = name.split("_")
                rid = ResidueId(int(rid))

                D = {
                    'rName': chain[rid].name.str, 'aName': aname, 'rId': rid.serial, 'limit': limit
                }

            if with_constant:
                D.update({"constant": popt[-1]})
                popt = popt[:-1]

            amplitudes = popt[::2]
            taus = popt[1::2]

            D.update(
                {("exp-%d-a%d" % (order, i + 1)): a for i, a in enumerate(amplitudes)}
            )
            D.update(
                {("exp-%d-tau%d" % (order, i + 1)): tau for i, tau in enumerate(taus)}
            )

            temp = pd.DataFrame(D, index=[0])

            tau_table = pd.concat([tau_table, temp])
        if tumbling:
            tau_table.sort_values(by=['axis'], inplace=True)
        else:
            tau_table.sort_values(by=['rId'], inplace=True)
        tau_table.to_csv(os.path.join(output_directory, 'tau_%d_exp.csv' % order), index=False)
