from bionmr_utils.md import  *
from tqdm import tqdm
import argparse
import numpy as np
import os
import pandas as pd


def get_msd(array):
    shifts = np.arange(len(array))
    msd = np.zeros(shifts.size)
    for i, shift in enumerate(shifts):
        diffs = array[:-shift if shift else None] - array[shift:]
        sqdist = np.square(diffs).sum(axis=1)
        msd[i] = sqdist.mean()
    return msd


def extract_mass_center_without_shift(path_to_trajectory: str, output_directory: str, trajectory_length: int=1):
    atomics_weight = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.973762, 'S': 32.06}
    """

    :param path_to_trajectory:
    :param output_directory: name of output four-column .csv file [time_ns,x,y,z]
    """
    traj, ref = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    masses = [atomics_weight[atom.name.str[0]] for atom in ref.asAtoms]
    center_all_without_shift = []
    for frame in  tqdm(traj):
        atoms = frame.asAtoms
        center = atoms.mass_center(masses)
        center_all_without_shift.append(center)
    return center_all_without_shift


def extract_mass_center_with_shift(path_to_trajectory: str, output_directory: str, volume=None):
    latice_path = os.path.join(path_to_trajectory,"1_build", "box.inpcrd")
    with open(latice_path, "r") as in_file:
        last_line = in_file.readlines()[-1].split()
        vectors = [float(coordinate) for coordinate in last_line[0:3]]
        angles = [Degrees(float(coordinate)) for coordinate in last_line[3:]]
        print(angles)
    latice = LatticeVectors(vectors[0], vectors[1], vectors[2],
                            angles[0], angles[1], angles[2])
    if volume:
        volume_0 = np.dot(np.cross(latice[0].to_np, latice[1].to_np, axisc=0), latice[2].to_np)
        factor_scale = (float(volume[0]) / float(volume_0)) ** (1 / 3)
        latice.scale_by(factor_scale)
    else:
        basis = BestShiftFinder(latice)
    center_all_with_shift = pd.DataFrame()
    center_all_without_shift = extract_mass_center_without_shift(path_to_trajectory, output_directory)
    prev_center = None
    for center in center_all_without_shift:
        if prev_center is None:
            prev_center = center_all_without_shift[0]
        else:
            if volume:
                factor_scale = (float(volume[ind]) / float(volume[ind - 1])) ** (1 / 3)
                latice.scale_by(factor_scale)
                basis = BestShiftFinder(latice)
            temp = pd.DataFrame({'x': prev_center.x, 'y': prev_center.y, 'z': prev_center.z}, index=[0],
                                columns=['x', 'y', 'z'])
            center_all_with_shift = pd.concat([center_all_with_shift, temp])
            current_center = center
            length, shift = basis.find_best_shift(prev_center, current_center)
            prev_center = current_center + shift

    pd.DataFrame(np.array([np.linspace(0, len(center_all_with_shift) * 0.002, len(center_all_with_shift), endpoint=False),
                 center_all_with_shift.x, center_all_with_shift.y, center_all_with_shift.z]).T,
                 columns=["time_ns", "x", "y", "z"]).to_csv(os.path.join(output_directory, "center_mass_coordinates.csv"), index=False)

    return center_all_with_shift


def extract_msd(path_to_trajectory: str, output_directory: str, volume=None):
    center_all_with_shift = extract_mass_center_with_shift(
        path_to_trajectory, output_directory, volume = None)

    arr = np.array(center_all_with_shift)
    msd = get_msd(arr)
    T = np.linspace(0, len(msd) * 0.002, len(msd), endpoint=False)
    df = pd.DataFrame(np.array([T, msd]).T, columns=["time_ns", "delta_rsquare"])
    df.to_csv(os.path.join(output_directory, "msd.csv"))
    return msd


if __name__ == '__main__':
    # -i "/home/olebedenko/bioinf/1ubq/1ubq/"
    parser = argparse.ArgumentParser(description='Calc NH autocorr')
    parser.add_argument('-i', '--path_to_trajectory', required=True, )
    parser.add_argument('-o', '--output_directory', default=os.getcwd())
    parser.add_argument('-v', '--volume', default=None, type=int)
    args = parser.parse_args()
    extract_msd(args.path_to_trajectory, args.output_directory, args.volume)
   