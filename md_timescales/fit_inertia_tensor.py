from typing import List
import argparse
import numpy as np
from bionmr_utils.md import *
import os
import pandas as pd
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
    return limit, popt


def save_fit_auto_correlation(path_to_csv_acorr: str, output_directory: str):

    curve_bounds = {
        2: ([[0, 0.1, 0, 1], [1, 1, 1, 10]]),
        3: ([[0, 0.01, 0, 0.1, 0, 1], [1, 0.1, 1, 1, 1, 10]]),
        4: ([[0, 0.001, 0, 0.01, 0, 1, 0, 10], [1, 0.01, 1, 0.1, 1, 10, 1, 100]]),
    }
    for order in range(2, 5):
        inertia_tensor_tau_table = pd.DataFrame()
        for inertia_axis in range(1,4):
	        df = pd.read_csv(os.path.join(path_to_csv_acorr, "overall_tumbling_%d.csv"%(inertia_axis)))
	        limit, popt = fit_auto_correlation(df.time_ns,
	                                    df.acorr, order,
	                                    with_constant=False,
	                                    curve_bounds=curve_bounds[order])
	        amplitudes = popt[::2]
	        taus = popt[1::2]

	        D = {
	            'inertia_axis': inertia_axis, 'limit': limit
	        }
	        D.update(
	            {("exp-%d-a%d" % (order, i + 1)): a for i, a in enumerate(amplitudes)}
	        )
	        D.update(
	            {("exp-%d-tau%d" % (order, i + 1)): tau for i, tau in enumerate(taus)}
	        )

	        temp = pd.DataFrame(D, index=[0])

        	inertia_tensor_tau_table = pd.concat([inertia_tensor_tau_table, temp] )

        inertia_tensor_tau_table = inertia_tensor_tau_table.sort_values(by=['inertia_axis'])
        inertia_tensor_tau_table.to_csv(os.path.join(output_directory, 'tau_inertia_tensor_%d_exp.csv' % order), index=False)


if __name__ == '__main__':

  #-accor "/home/legosta/bioinf/scripts/data/acorr"
  parser = argparse.ArgumentParser(description='fit inertia tensor autocorr')
  parser.add_argument('-accor', '--path_to_csv_acorr', required=True)
  parser.add_argument('-o', '--output_directory', default=os.getcwd())
  args = parser.parse_args()
  save_fit_auto_correlation(args.path_to_csv_acorr, args.output_directory)