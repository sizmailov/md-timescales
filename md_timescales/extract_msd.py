from bionmr_utils.md import *
from tqdm import tqdm
import argparse
import numpy as np
import os
import pandas as pd


def decrease_points_in_data(array):
    # end_point = int(np.log2(len(array)))
    # index = np.logspace(0, end_point, num=end_point + 1, endpoint=True, base=2)
    # new_array = np.array([array[int(j)] for j in index])
    new_array = np.array([array[i] for i in range(0,len(array), 10)])
    return new_array


def get_msd(array):
    new_array = decrease_points_in_data(array)
    shifts = np.arange(len(new_array))
    msd = np.zeros(shifts.size)
    for i, shift in enumerate(tqdm(shifts)):
        diffs = new_array[:-shift if shift else None] - new_array[shift:]
        sqdist = np.square(diffs).sum(axis=1)
        msd[i] = sqdist.mean()
    return msd


def extract_center_mass_from_trajectory(path_to_trajectory: str, trajectory_length: int = 1):
    """
    :param path_to_trajectory:
    :param trajectory_length: by default equal 1
    """
    atomics_weight = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.973762, 'S': 32.06}
    traj, ref = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    masses = [atomics_weight[atom.name.str[0]] for atom in ref.asAtoms]
    center_all_without_shift = VectorXYZ()
    for frame in traj:
        atoms = frame.asAtoms
        center = atoms.mass_center(masses)
        center_all_without_shift.append(center)
    return center_all_without_shift


def extract_mass_center_with_shift(path_to_trajectory: str, trajectory_length: int, volume=None):
    lattice_path = os.path.join(path_to_trajectory, "1_build", "box.inpcrd")
    with open(lattice_path, "r") as in_file:
        last_line = in_file.readlines()[-1].split()
        vectors = [float(coordinate) for coordinate in last_line[0:3]]
        angles = [Degrees(float(coordinate)) for coordinate in last_line[3:]]
    latice = LatticeVectors(vectors[0], vectors[1], vectors[2],
                            angles[0], angles[1], angles[2])
    if volume:
        volume_0 = (latice[0].cross(latice[1])).dot(latice[2])
        factor_scale = (float(volume[0]) / float(volume_0)) ** (1 / 3)
        latice.scale_by(factor_scale)
    else:
        basis = BestShiftFinder(latice)

    center_all_with_shift = pd.DataFrame()
    center_all_without_shift = extract_center_mass_from_trajectory(path_to_trajectory, trajectory_length)
    prev_center = None
    center_all_with_shift = VectorXYZ()

    for ind, center in tqdm(enumerate(center_all_without_shift)):
        if prev_center is None:
            prev_center = center_all_without_shift[ind]
        else:
            if volume:
                factor_scale = (float(volume[ind]) / float(volume[ind - 1])) ** (1 / 3)
                latice.scale_by(factor_scale)
                basis = BestShiftFinder(latice)
            center_all_with_shift.append(prev_center)
            current_center = center
            length, shift = basis.find_best_shift(prev_center, current_center)
            prev_center = current_center + shift
    center_all_with_shift = center_all_with_shift.to_numpy()
    return center_all_with_shift


def extract_msd(path_to_trajectory: str, output_directory: str, trajectory_length: int, volume=None):
    center_all_with_shift = extract_mass_center_with_shift(
        path_to_trajectory, trajectory_length, volume=None)

    arr = np.array(center_all_with_shift)
    msd = get_msd(arr)
    T = np.linspace(0, len(msd) * 0.002, len(msd), endpoint=False)
    df = pd.DataFrame(np.array([T, msd]).T, columns=["time_ns", "delta_rsquare"])
    df.to_csv(os.path.join(output_directory, "msd.csv"))
    return msd


if __name__ == '__main__':
    # -i "/home/olebedenko/bioinf/1ubq/1ubq/"
    parser = argparse.ArgumentParser(description='Mean square displacement centr mass (msd)')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--output-directory', default=os.getcwd())
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--volume', default=False)
    args = parser.parse_args()
    extract_msd(args.path_to_trajectory, args.output_directory, args.trajectory_length, args.volume)
