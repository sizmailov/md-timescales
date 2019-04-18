from bionmr_utils.md import *
import argparse
import os
import pandas as pd
from tqdm import tqdm
from md_timescales.extract import extract_mass_center, extract_lattice_vectors_rst7, extract_time_step_ns

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract mass center')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--output-directory', default=os.getcwd())
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--volume-file', type=str, default=None)
    parser.add_argument('--lattice-rst7-file', type=str, default=None)

    args = parser.parse_args()

    if args.volume_file and args.volume_file != "None":
        volume = pd.read_csv(args.volume_file, sep=r"\s+").values[:, 1]
        lattice_vectors = extract_lattice_vectors_rst7(os.path.join(args.path_to_trajectory, "1_build", "box.inpcrd"))
    else:
        assert args.lattice_rst7_file != "None"
        volume = None
        lattice_vectors = extract_lattice_vectors_rst7(args.lattice_rst7_file)

    time, cm = extract_mass_center(
        tqdm(traj_from_dir(args.path_to_trajectory, last=args.trajectory_length)[0]),
        dt=extract_time_step_ns(args.path_to_trajectory),
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
