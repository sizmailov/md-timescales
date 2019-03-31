import argparse
from md_timescales.extract import extract_autocorr, get_NH_vectors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calc NH autocorr')
    parser.add_argument('-i', '--path_to_trajectory', required=True)
    parser.add_argument('-o', '--output_directory', default="./")
    parser.add_argument('-l', '--length_trajectory', default=1, type=int)
    args = parser.parse_args()

    extract_autocorr(args.path_to_trajectory,
                     args.output_directory,
                     trajectory_length=args.length_trajectory,
                     get_vectors=get_NH_vectors)
