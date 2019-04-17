import argparse
from md_timescales.fit import save_fit_auto_correlation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit CH3 autocorr')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--path-to-acorrs', required=True)
    parser.add_argument('--output-directory', default="./")
    parser.add_argument('--vectors-group', choices=["NH", "CH3"], type=str)
    args = parser.parse_args()

    bounds = {
        "NH": [
            ([[0, 0.0001, 0, 0.002], [1, 0.002, 1, 5]]),
            ([[0, 0.0001, 0, 0.002, 0, 5], [1, 0.002, 1, 5, 1,10]]),
            ([[0, 0.0001, 0, 0.002, 0, 5, 0, 10], [  1, 0.002, 1, 10, 1, 10, 1, 20]]),
        ],
        "CH3": [
            ([[0, 0.002, 0], [1, 0.5, 1]]),
            ([[0, 0, 0, 0.002, 0], [1, 0.002, 1, 1.5, 1]]),
            ( [[0, 0, 0, 0.002, 0, 0.002, 0], [1, 0.002, 1, 0.5, 1, 0.5, 1]])

        ]
    }
    save_fit_auto_correlation(args.path_to_trajectory,
                              args.path_to_acorrs,
                              args.output_directory,
                              curve_bounds=bounds[args.vectors_group]
                              )
