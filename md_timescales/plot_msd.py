import argparse
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from md_timescales.plot import add_relpath_to_top_corner
import pandas as pd
import numpy as np
import os


def get_translation_fit(path_to_msd, path_to_fit, output_directory):
    assert os.path.isfile(path_to_msd)
    assert os.path.isfile(path_to_fit)
    df_msd = pd.read_csv(path_to_msd)

    fig = plt.figure()  # type:plt.Figure
    add_relpath_to_top_corner(fig)

    ax = fig.add_subplot(111)  # type: plt.Axes
    ax.set_xlabel('time, ns', fontsize=13)
    ax.set_ylabel('msd, A^2', fontsize=13)
    ax.set_title('Mean square displacement center mass (msd)')
    ax.plot(df_msd.time_ns, df_msd.msd, label="")
    fit = pd.read_csv(path_to_fit)
    p_coeff = [fit["a1"][0], fit["a0"][0]]
    P = np.poly1d(p_coeff)
    ax.plot(df_msd.time_ns, P(df_msd.time_ns), label="D = {:.3e} $m^2/s$".format(fit["a1"][0] * 1e-10 ** 2 / 1e-9 / 6))
    ax.axvline(fit["fit_limit_ns"][0], lw=1, ls="--", color="green")
    ax.legend()
    ax.grid(True)
    plt.savefig(os.path.join(output_directory, "fit_msd_plot.png"))
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot msd autocorrelation')
    parser.add_argument('--path-to-msd', required=True)
    parser.add_argument('--path-to-fit', required=True)
    parser.add_argument('--output-directory', default=os.getcwd())
    args = parser.parse_args()
    get_translation_fit(args.path_to_msd, args.path_to_fit, args.output_directory)
