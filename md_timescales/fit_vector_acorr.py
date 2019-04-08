import argparse
from fit import save_fit_auto_correlation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit CH3 autocorr')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--path-to-acorrs', required=True)
    parser.add_argument('--output-directory', default="./")
    parser.add_argument('--vectors-group', choices=["NH", "CH3"], type=str)
    args = parser.parse_args()

    bounds = {
        "NH": [
            ([[0.05, 0.01, 0.05, 1], [1, 1, 1, 10]]),
            ([[0.05, 0.01, 0.05, 0.1, 0.05, 1], [1, 0.1, 1, 1, 1, 10]]),
            ([[0.05, 0.001, 0.05, 0.01, 0.05, 1, 0.05, 10], [1, 0.01, 1, 0.1, 1, 10, 1, 100]]),
        ],
        "CH3": [
            ([[0, 1, 0], [1, 10, 1]]),
            ([[0, 0.01, 0, 0.1, 0], [1, 0.1, 1, 1, 1]]),
            ([[0, 0.001, 0, 0.01, 0, 1, 0], [1, 0.01, 1, 0.1, 1, 10, 1]])
        ]
    }
    save_fit_auto_correlation(args.path_to_trajectory,
                              args.path_to_acorrs,
                              args.output_directory,
                              curve_bounds=bounds[args.vectors_group]
                              )

