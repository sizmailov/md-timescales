import argparse
from md_timescales.plot import plot_acorr_fit

if __name__ == '__main__':

    #-f_csv "/home/olebedenko/bioinf/scripts/md-timescales/md_timescales" -a_csv "/home/olebedenko/bioinf/handling/h4/tip4p-ew/NPT_gamma_ln_2/autocorr/CH3/data"
    parser = argparse.ArgumentParser(description="Plot autocorrelation function")
    parser.add_argument('--path-to-fit', required=True)
    parser.add_argument('--path-to-acorrs', required=True)
    parser.add_argument('--output-directory', default="./")
    args = parser.parse_args()
    plot_acorr_fit(args.path_to_fit, args.path_to_acorrs, args.output_directory)
