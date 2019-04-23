from typing import *
import numpy as np
from bionmr_utils.md import *
import os
import pandas as pd
import glob
from scipy.optimize import curve_fit


def fit_limit(data: List[Union[float, int]]) -> int:
    """

    :param data: array of function values
    :return: index separates region with small
             amount of negative derivatives
    """

    def moving_average(data_set, periods=3):
        weights = np.ones(periods) / periods
        return np.convolve(data_set, weights, mode='valid')

    diff = np.diff(data)
    pos_diff = (diff > 0).astype(int)
    window_size = 50
    pos_diff_avg = moving_average(pos_diff, window_size)
    index = np.argmax(pos_diff_avg > 0.5) + window_size // 2
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


def generate_bounds(path_to_csv_acorr: str, with_constant: bool = False):
    """
    Function finds boundaries for fit parametrs of autocorrelation function 
    for particular trajectory. Selection is based on fix for point of loosing 
    autocorelation values: 99%, 50%, 20% and 1%.
    Generated boundaries can use for two, three and for exponential fit

    :param path_to_csv_acorr: path to two-column .csv files [time_ns, acorr]
    :param with_constant: condition of fit with constant
    """

    df = pd.read_csv(path_to_csv_acorr)
    tau0 = 0
    tau1 = df.time_ns[df.acorr < 0.9999].iloc[0]
    tau2 = df.time_ns[df.acorr < 0.5].iloc[0]
    tau3 = df.time_ns[df.acorr < 0.2].iloc[0]
    tau4 = df.time_ns[df.acorr < 0.01].iloc[0]
    tau_left = [tau0, tau1, tau2, tau3]
    tau_right = [tau1, tau2, tau3, tau4]
    amp0 = [0, 0, 0, 0]
    amp1 = [1, 1, 1, 1]
    if with_constant:
        tau_left = [tau0, tau1, tau1]
        tau_right = [tau1, tau4, tau4]

    left = [None] * (len(amp0) + len(tau_left))
    left[::2] = amp0
    left[1::2] = tau_left

    right = [None] * (len(amp1) + len(tau_right))
    right[::2] = amp1
    right[1::2] = tau_right

    bounds_4 = [left, right]
    bounds_3 = [left[0:6 - int(with_constant)], right[0:6 - int(with_constant)]]
    bounds_2 = [left[0:4 - int(with_constant)], right[0:4 - int(with_constant)]]

    if with_constant:
        bounds_2[0][1] = tau1
        bounds_2[1][1] = tau4

    bounds = (bounds_2, bounds_3, bounds_4)
    return bounds


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


def get_fit_auto_correlation(ref_chain: Chain,
                             csv_files: List[str],
                             output_directory: str,
                             with_constant: bool = False,
                             tumbling: bool = False,
                             ) -> None:
    """

    :param ref_chain: chain of reference
    :param csv_files:  two-column .csv files [time_ns, acorr]
    :param output_directory: output directory for 
           a particular fit function (e.g. tau-2-exp.csv)
    :param curve_bounds: restriction of function parameters
    :param ca_alignment: flag of aligment frames by Ca atoms
    :param tumbling: flag of tumbling calculation
    """
    curve_bounds = generate_bounds(csv_files[0], with_constant)
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
                              with_constant: bool = False,
                              ca_alignment: bool = False,
                              tumbling: bool = False
                              ) -> None:
    """

    :param path_to_ref: path to reference
    :param path_to_csv_acorr: path to two-column .csv files [time_ns, acorr]
    :param output_directory: output directory for 
           a particular fit function (e.g. tau-2-exp.csv)
    :param curve_bounds: restriction of function parameters
    :param ca_alignment: flag of aligment frames by Ca atoms
    :param tumbling: flag of tumbling calculation
    """
    traj, ref = traj_from_dir(path_to_ref, first=1, last=1)
    ref_chain = ref.asChains[0]
    csv_files = sorted(glob.glob(os.path.join(path_to_csv_acorr, "*.csv")))
    get_fit_auto_correlation(ref_chain, csv_files, output_directory,
                             with_constant, tumbling)

    if ca_alignment:
        path_to_csv_acorr = os.path.join(path_to_csv_acorr, "ca_alignment")
        output_directory = os.path.join(output_directory, "ca_alignment")
        os.makedirs(output_directory, exist_ok=True)
        csv_files = sorted(glob.glob(os.path.join(path_to_csv_acorr, "*.csv")))
        get_fit_auto_correlation(ref_chain, csv_files, output_directory,
                                 with_constant, tumbling)
