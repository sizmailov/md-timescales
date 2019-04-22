import argparse
from md_timescales.fit import save_fit_auto_correlation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit CH3 autocorr')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--path-to-acorrs', required=True)
    parser.add_argument('--output-directory', default="./")
    parser.add_argument('--vectors-group', choices=["NH", "CH3"], type=str)
    args = parser.parse_args()

    with_constant = {"NH": False,
                    "CH3": True}
    ca_alignment = {"NH": False,
                    "CH3": True}
    save_fit_auto_correlation(args.path_to_trajectory,
                              args.path_to_acorrs,
                              args.output_directory, 
                              with_constant=with_constant[args.vectors_group],
                              ca_alignment=ca_alignment[args.vectors_group]
                              )
