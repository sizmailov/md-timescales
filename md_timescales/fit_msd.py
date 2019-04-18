import argparse
import pandas as pd
import numpy as np
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit CH3 autocorr')
    parser.add_argument('--path-to-msd', required=True)
    parser.add_argument('--output-directory', default="./")
    parser.add_argument('--fit-fraction', type=float, default=0.25)
    args = parser.parse_args()

    df = pd.read_csv(os.path.join(args.path_to_msd, "msd.csv"))
    N = int(len(df)*args.fit_fraction)
    p_coeff = np.polyfit(df["time_ns"][:N], df["msd"][:N], 1)

    pd.DataFrame({
        "a1": [p_coeff[0]],
        "a0": [p_coeff[1]]
    }).to_csv(os.path.join(args.output_directory, "fit.csv"),index=False)
