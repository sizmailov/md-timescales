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


def mono_exp(x, a1, tau1):
  return a1 * np.exp(-x / tau1) 

def fit_auto_correlation(time: List[float], acorr: List[float], with_constant: bool,
                         curve_bounds) -> List[float]:
    limit = fit_limit(acorr)
    popt, pcov = curve_fit(mono_exp, time[:limit], acorr[:limit], bounds=curve_bounds)
    return limit, popt

def save_fit_auto_correlation(path_to_csv_acorr: str, output_directory: str):
    curve_bounds = {1: ([[0, 1], [1, 10]])}
    inertia_tensor_tau_table = pd.DataFrame()
    for inertia_axis in range(1,4):
        df = pd.read_csv(os.path.join(path_to_csv_acorr, "overall_tumbling_%d.csv"%(inertia_axis)))
        limit, popt = fit_auto_correlation(df.time_ns,df.acorr, with_constant=False, curve_bounds=curve_bounds[1])

        amplitudes = popt[::2]
        taus = popt[1::2]

        D = {'inertia_axis': inertia_axis, 'limit': limit}

        D.update({("exp-%d-a%d" % (1,1)): i for i in amplitudes})
        D.update({("exp-%d-tau%d" % (1,1)): j for j in taus})


        temp = pd.DataFrame(D, index=[0])

        inertia_tensor_tau_table = pd.concat([inertia_tensor_tau_table, temp])

    inertia_tensor_tau_table.to_csv(os.path.join(output_directory, 'tau_inertia_tensor_%d_exp.csv' % 1), index=False)

if __name__ == '__main__':

  #-accor "/home/legosta/bioinf/scripts/data/acorr"
  parser = argparse.ArgumentParser(description='fit inertia tensor autocorr')
  parser.add_argument('-acorr', '--path_to_csv_acorr', required=True)
  parser.add_argument('-o', '--output_directory', default=os.getcwd())
  args = parser.parse_args()
  save_fit_auto_correlation(args.path_to_csv_acorr, args.output_directory)
