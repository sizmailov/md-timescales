import argparse
from md_timescales.extract import extract_autocorr, get_NH_vectors, get_methyl_vectors

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calc autocorr for interactomic vectors')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--output-directory', required=True)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--vectors-group', choices=["NH", "CH3"], type=str)
    args = parser.parse_args()

    get_vectors = {
        "NH": get_NH_vectors,
        "CH3": get_methyl_vectors,
    }

    extract_autocorr(args.path_to_trajectory,
                     args.output_directory,
                     trajectory_length=args.trajectory_length,
                     get_vectors=get_vectors[args.vectors_group])
