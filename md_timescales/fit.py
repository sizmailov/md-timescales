from typing import *
import numpy as np
from bionmr_utils.md import *
import os
import pandas as pd
import glob
from scipy.optimize import curve_fit


def fit_limit(data: List[Union[float, int]], window_size=50, pos_diff_ratio=0.5) -> int:
    """
    Returns minimum of
    1) number of points at the beginning of `data` for which in every
    `window_size` range there are no more than `pos_diff_ratio` positive
    derivatives
    2) first data point which cross zero

    :param data: array of function values
    :param window_size: number of points in window
    :param pos_diff_ratio: allowed ratio of positive diff in window
    :return: index separates region with small
             amount of negative derivatives
    """

    def moving_average(data_set, periods=3):
        weights = np.ones(periods) / periods
        return np.convolve(data_set, weights, mode='valid')

    diff = np.diff(data)
    pos_diff = (diff > 0).astype(int)
    pos_diff_avg = moving_average(pos_diff, window_size)
    index = np.argmax(pos_diff_avg >= pos_diff_ratio) + window_size // 2
    first_negative = np.array(data < 0).argmax()
    if data[first_negative] >= 0:
        first_negative = len(data)
    index = min(index, first_negative)
    return index


def __multi_exp_f(x: Union[float, int],
                  A: List[Union[float, int]],
                  TAU: List[Union[float, int]],
                  C: Union[float, int]) -> float:
    """
    :param x: argument of some exponential functions composition
    :param A: array of amplitudes
    :param TAU: array of time constants
    :param C: free element
    :return: sum exponential functions composition
    """
    return np.sum(
        (a * np.exp(-x / tau)) for a, tau in zip(A, TAU)
    ) + C


def multi_exp(x: Union[float, int],
              *args):
    """
    :param x: argument of some exponential functions composition
    :param args: array of amplitudes and time constants
    :return: callable __multi_exp
    """
    TAU = list(args[1::2])

    if len(args) % 2 == 0:
        C = 0
        A = list(args[::2])
    else:
        C = args[-1]
        A = args[:-1:2]

    return __multi_exp_f(x, A, TAU, C)


def mult_exp_sum_1(x, *args):
    TAU = args[0::2]

    if len(args) % 2 == 1:
        C = 0
        A = args[1:-1:2]
    else:
        C = args[-1]
        A = args[1:-1:2]
    A0 = 1 - sum(A) - C
    return __multi_exp_f(x, [A0] + list(A), TAU, C)


def fit_auto_correlation(time: List[float],
                         acorr: List[float],
                         bounds: List[List[List[Union[float, int]]]]) \
        -> Tuple[int, Union[np.ndarray, Iterable, int, float]]:
    """
    Fit input data with :math:`\sum_n A_n \exp(-t/\\tau_n) + const`

    :param time: time data series
    :param acorr: auto-correlation data series
    :param bounds: curve parameters bounds
    :return: Fit curve parameters
    """

    p0 = np.mean(bounds, axis=0)[1:]

    args, pcov = curve_fit(mult_exp_sum_1,
                           time,
                           acorr,
                           p0=p0,
                           bounds=np.array(bounds)[:, 1:])

    if len(args) % 2 == 1:
        C = 0
        A = args[1:-1:2]
    else:
        C = args[-1]
        A = args[1:-1:2]
    A0 = 1 - sum(A) - C

    return [A0] + list(args)


def decorated_fit_auto_correlation(time: List[float],
                                   acorr: List[float],
                                   bounds: List[List[Union[float, int]]],
                                   window_size=50,
                                   pos_diff_ratio=0.5,
                                   ) \
        -> Tuple[int, Union[np.ndarray, Iterable, int, float]]:
    def scale_times(args, scale):
        args[1::2] = np.array(args[1::2]) * scale

    scales = [1, 2, 3]

    limit = fit_limit(acorr,
                      window_size=window_size,
                      pos_diff_ratio=pos_diff_ratio
                      )

    time = time[:limit]
    acorr = acorr[:limit]

    R_square = []
    popt_all = []

    for i, scale in enumerate(scales):
        scale_times(bounds[0], scale)
        scale_times(bounds[1], scale)
        try:
            popt = fit_auto_correlation(time, acorr, bounds)

            R_square.append(np.sum((np.array(acorr) -
                                    np.array(multi_exp(time, *popt))) ** 2))
            popt_all.append(popt)
        except RuntimeError:
            print("Fit error n={}, scale={}".format(len(bounds[0]) // 2, scale))

        scale_times(bounds[0], 1.0 / scale)
        scale_times(bounds[1], 1.0 / scale)

    min_ind_r_square = np.argmin(R_square)

    return limit, popt_all[min_ind_r_square]


def fit_mean_square_displacement(time: List[float], msd: List[float]) -> List[float]:
    """
    Fit input data with :math:` a  t + b`

    :param time: time data series
    :param msd: mean square displacement data series
    :return: Fit curve parameters
    """
    ...


def get_fit_auto_correlation(ref_chain: Chain,
                             csv_files: List[str],
                             output_directory: str,
                             curve_bounds: List[List[List[Union[float, int]]]],
                             tumbling: bool = False,
                             window_size=50,
                             pos_diff_ratio=0.5
                             ) -> None:
    """

    :param ref_chain: chain of reference
    :param csv_files:  two-column .csv files [time_ns, acorr]
    :param output_directory: output directory for 
           a particular fit function (e.g. tau-2-exp.csv)
    :param curve_bounds: restriction of function parameters
    :param ca_alignment: flag of aligment frames by Ca atoms
    :param tumbling: flag of tumbling calculation
    :param window_size: window size to detect fit limit
    :param pos_diff_ratio: positive diff ratio to detect fit limit
    """
    for bounds in curve_bounds:
        with_constant = len(bounds[0]) % 2 == 1
        order = (len(bounds[0]) + 1) // 2
        tau_table = pd.DataFrame()
        for file in csv_files:

            df = pd.read_csv(file)
            limit, popt = decorated_fit_auto_correlation(df.time_ns,
                                                         df.acorr,
                                                         bounds=bounds,
                                                         window_size=window_size,
                                                         pos_diff_ratio=pos_diff_ratio
                                                         )

            assert np.isclose(np.sum(popt[::2]), 1.0), "Sum of exponential amplitudes != 1 (sum Ai = {:.3f})".format(
                np.sum(popt[::2]))

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
                    'rName': ref_chain[rid].name.str, 'aName': aname, 'rId': rid.serial, 'limit': limit
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


def save_fit_auto_correlation(path_to_ref: str,
                              path_to_csv_acorr: str,
                              output_directory: str,
                              curve_bounds: List[List[List[Union[float, int]]]],
                              ca_alignment: bool = False,
                              tumbling: bool = False,
                              window_size=50,
                              pos_diff_ratio=0.5
                              ) -> None:
    """

    :param path_to_ref: path to reference
    :param path_to_csv_acorr: path to two-column .csv files [time_ns, acorr]
    :param output_directory: output directory for 
           a particular fit function (e.g. tau-2-exp.csv)
    :param curve_bounds: restriction of function parameters
    :param ca_alignment: flag of alignment frames by Ca atoms
    :param tumbling: flag of tumbling calculation
    :param window_size: window size to detect fit limit
    :param pos_diff_ratio: positive diff ratio to detect fit limit
    """
    traj, ref = traj_from_dir(path_to_ref, first=1, last=1)
    ref_chain = ref.asChains[0]
    csv_files = sorted(glob.glob(os.path.join(path_to_csv_acorr, "*.csv")))
    get_fit_auto_correlation(ref_chain,
                             csv_files,
                             output_directory,
                             curve_bounds,
                             tumbling,
                             window_size=window_size,
                             pos_diff_ratio=pos_diff_ratio)

    if ca_alignment:
        path_to_csv_acorr = os.path.join(path_to_csv_acorr, "ca_alignment")
        output_directory = os.path.join(output_directory, "ca_alignment")
        os.makedirs(output_directory, exist_ok=True)
        csv_files = sorted(glob.glob(os.path.join(path_to_csv_acorr, "*.csv")))
        get_fit_auto_correlation(ref_chain,
                                 csv_files,
                                 output_directory,
                                 curve_bounds,
                                 tumbling,
                                 window_size=window_size,
                                 pos_diff_ratio=pos_diff_ratio)
