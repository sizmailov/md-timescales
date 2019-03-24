from typing import List
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


def two_exp(x, a1, tau1, a2, tau2):
    return a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2)


def three_exp(x, a1, tau1, a2, tau2, a3, tau3):
    return a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2) + a3 * np.exp(-x / tau3)


def four_exp(x, a1, tau1, a2, tau2, a3, tau3, a4, tau4):
    return a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2) + a3 * np.exp(-x / tau3) + a4 * np.exp(-x / tau4)


def fit_auto_correlation(time: List[float], acorr: List[float], order: int, with_constant: bool,
                         curve_bounds) -> List[float]:
    """
    Fit input data with :math:`\sum_n A_n \exp(-t/\\tau_n) + const`

    :param time: time data series
    :param acorr: auto-correlation data series
    :param order: number of exponentials in fit
    :param with_constant: if false const = 0
    :param curve_bounds: curve bounds parametrs
    :return: Fit curve parameters
    """
    fit_func = {2: two_exp, 3: three_exp, 4: four_exp}
    limit = fit_limit(acorr)
    popt, pcov = curve_fit(fit_func[order], time[:limit], acorr[:limit],
                           bounds=curve_bounds)
    return popt


def save_fit_auto_correlation(ref: Frame, path_to_csv_accor: str, output_directory: str):
    NH_tau_table = pd.DataFrame()
    csv_files = sorted(glob.glob(os.path.join(path_to_csv_accor, "*.csv")))
    # columns = {4: {'rName': chain[int(rId)].name.str,
    #               'aName': chain[int(rId)].asAtoms[0].aName.str,
    #               'rId': int(rId),
    #               'exp3-a1': popt[0], 'exp3-tau1': popt[1],
    #               'exp3-a2': popt[2], 'exp3-tau2' : popt[3],
    #               'exp3-a3' : popt[4], 'exp3-tau3' : popt[5],
    #               'exp4-a4' : popt[6], 'exp4-tau4' : popt[7]}
    #            3: {'rName': chain[int(rId)].name.str,
    #               'aName': chain[int(rId)].asAtoms[0].aName.str,
    #               'rId': int(rId),
    #               'exp3-a1': popt[0], 'exp3-tau1': popt[1],
    #               'exp3-a2': popt[2], 'exp3-tau2' : popt[3],
    #               'exp3-a3' : popt[4], 'exp3-tau3' : popt[5]},
    #            2: {'rName': chain[int(rId)].name.str,
    #               'aName': chain[int(rId)].asAtoms[0].aName.str,
    #               'rId': int(rId),
    #               'exp3-a1': popt[0], 'exp3-tau1': popt[1],
    #               'exp3-a2': popt[2], 'exp3-tau2' : popt[3]}}

    curve_bounds = {
        2: (-np.inf, np.inf),  # TODO: set actual curve bounds
        3: (-np.inf, np.inf),
        4: (-np.inf, np.inf),
    }
    chain = ref.asChains[0]
    for order in range(2, 5):
        name = os.path.splitext(os.path.basename(file))[0]
        rId = ResidueId(name)
        df = pd.read_csv(file)
        for file in csv_files:
            popt = fit_auto_correlation(df.time_ns,
                                        df.acorr, order,
                                        with_constant=False,
                                        curve_bounds=curve_bounds[order])
            amplitudes = popt[::2]
            taus = popt[1::2]

            D = {
                'rName': chain[rId].name.str, 'aName': chain[rId].asAtoms[0].aName.str, 'rId': rId.serial
            }
            D.update(
                {"exp-%d-a%d" % (order, a) for a in amplitudes}
            )
            D.update(
                {"tau-%d-tau%d" % (order, tau) for tau in taus}
            )

            temp = pd.DataFrame(D)

            NH_tau_table = pd.concat([NH_tau_table, temp] )

        NH_tau_table = NH_tau_table.sort_values(by=['rId'])
        NH_tau_table.to_csv(os.path.join(output_directory, 'tau_NH_%d_exp.csv' % order), index=False)
