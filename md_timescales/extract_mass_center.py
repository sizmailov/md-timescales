from bionmr_utils.md import *
import argparse
import os
import pandas as pd
from md_timescales.extract import extract_mass_center, extract_lattice_vectors_rst7

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract mass center')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--output-directory', default=os.getcwd())
    parser.add_argument('--trajectory-length', required=True, type=int)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--volume-file', type=str, default=None)
    group.add_argument('--lattice-rst7-file', type=str, default=None)

    args = parser.parse_args()

    if args.volume_file:
        volume = pd.read_csv(args.volume_file, sep=r"\s+").values[:, 1]
        lattice_vectors = extract_lattice_vectors_rst7(os.path.join(args.path_to_trajectory, "1_build", "box.inpcrd"))
    else:
        volume = None
        lattice_vectors = extract_lattice_vectors_rst7(args.lattice_rst7_file)

    time, cm = extract_mass_center(
        traj_from_dir(args.path_to_trajectory, last=args.trajectory_length)[0],
        dt=0.002,
        lattice_vectors=lattice_vectors,
        volume=volume
    )

    cm = cm.to_numpy()

    pd.DataFrame({
        "time_ns": time,
        "cm_x": cm[:, 0],
        "cm_y": cm[:, 1],
        "cm_z": cm[:, 2],
    }).to_csv(os.path.join(args.output_directory, "cm.csv"), index=False)
