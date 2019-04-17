import argparse
import numpy as np
import os
import pandas as pd
from md_timescales.calc import calc_mean_square_displacement

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate mean square displacement of center of mass (msd)')

    parser.add_argument('--mass-center-csv')
    parser.add_argument('--output-directory')
    parser.add_argument('--lag-spacing', choices=["linear", "log"], default="log")
    parser.add_argument('--n-lag-points', default=1000, type=str)
    args = parser.parse_args()

    cm = pd.read_csv(args.mass_center_csv)
    time = cm["time_ns"]
    xyz = cm[["cm_x", "cm_z", "cm_z"]].values

    if args.lag_spacing == "linear":
        lag_index = np.unique(np.linspace(1, len(time), args.n_lag_points, endpoint=False).astype(int))
    else:
        lag_index = np.unique(np.logspace(0, np.log10(len(time) - 1), args.n_lag_points).astype(int))

    lag, msd = calc_mean_square_displacement(time, xyz, lag_index=lag_index)

    pd.DataFrame({
        "time_ns": lag,
        "msd": msd,
    }).to_csv(os.path.join(args.output_directory, "msd.csv"), index=False)
