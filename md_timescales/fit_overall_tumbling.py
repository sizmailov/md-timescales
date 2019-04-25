import argparse
from md_timescales.fit import save_fit_auto_correlation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit CH3 autocorr')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--path-to-acorrs', required=True)
    parser.add_argument('--output-directory', default="./")
    args = parser.parse_args()

    bounds =[
        ([[0.99, 0.002], [1, 25]])
    ]

    save_fit_auto_correlation(args.path_to_trajectory,
                              args.path_to_acorrs,
                              args.output_directory,
                              curve_bounds=bounds,
                              tumbling=True,
                              window_size=50,
                              pos_diff_ratio=0.5,
                              )
