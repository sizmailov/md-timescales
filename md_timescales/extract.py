import os
import pyxmolpp2
import pandas as pd
import numpy as np
from tqdm import tqdm

from typing import *
from bionmr_utils.md import *


def atom_name_mass_map(atom: Atom):
    atomics_weight = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.973762, 'S': 32.06}
    return atomics_weight[atom.name.str[0]]


def extract_lattice_vectors_rst7(lattice_path):
    with open(lattice_path, "r") as in_file:
        last_line = in_file.readlines()[-1].split()
        vectors = [float(coordinate) for coordinate in last_line[0:3]]
        angles = [Degrees(float(coordinate)) for coordinate in last_line[3:]]
    return LatticeVectors(vectors[0], vectors[1], vectors[2],
                          angles[0], angles[1], angles[2])


def extract_mass_center(traj: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice],
                        dt: float,
                        lattice_vectors: LatticeVectors,
                        volume: np.array = None,
                        atom_selector: Callable[[Atom], bool] = lambda _: True,
                        atom_mass_map: Callable[[Atom], float] = atom_name_mass_map
                        ) -> Tuple[np.array, VectorXYZ]:
    """
    Extract mass center

    :param traj: trajectory
    :param dt: time between trajectory frames
    :param lattice_vectors: translational vectors for periodic boundary conditions
    :param volume: cell volume along the trajectory. Must match *whole* trajectory in length.
                   If None no lattice_vectors rescale is performed
    :param atom_selector: selector for atom
    :param atom_mass_map: function to determinate mass of atom

    """
    time = []
    mass_centers = VectorXYZ()
    atoms = None
    prev_cm = None
    for frame in traj:
        if atoms is None:
            atoms = frame.asAtoms.filter(atom_selector)
            masses = [atom_mass_map(atom) for atom in atoms]
            bsf = BestShiftFinder(lattice_vectors)
            reference_volume = lattice_vectors[0].dot(lattice_vectors[1].cross(lattice_vectors[2]))

        if volume:
            factor_scale = (volume[frame.index] / reference_volume) ** (1 / 3)
            bsf.scale_lattice_by(factor_scale)

        current_cm = atoms.mass_center(masses)

        if prev_cm:
            current_cm += bsf.find_best_shift(prev_cm, current_cm)[1]

        time.append(frame.index * dt)
        mass_centers.append(current_cm)

        prev_cm = current_cm

        if volume:
            bsf.scale_lattice_by(1 / factor_scale)
    return np.array(time), mass_centers


def get_autocorr(trajectory: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice], get_vectors) -> dict:
    """
    :param trajectory
    """
    CH3S = None
    vectors = None
    for frame in tqdm(trajectory):
        if CH3S is None:
            CH3S = {}
            vectors = {}
            for A, B in get_vectors(frame):
                CH3S[A.rId, A.aName] = (A, B)
                vectors[A.rId, A.aName] = VectorXYZ()

        for (rid, aname), (A, B) in CH3S.items():
            vectors[rid, aname].append(A.r - B.r)

    autocorr = {
        (rid, aname): calc_autocorr_order_2(vector)
        for (rid, aname), vector in vectors.items()
    }

    return autocorr


def extract_autocorr(path_to_trajectory: str,
                     output_directory: str,
                     get_vectors: Callable[[Frame], List[Tuple[Atom, Atom]]],
                     trajectory_length: int
                     ) -> None:
    """

    :param get_vectors: returns list of atom pairs of interest
    :param path_to_trajectory:
    :param trajectory_length: number of .dat files to process
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    traj, ref = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    autocorr_CH3 = get_autocorr(traj, get_vectors=get_vectors)

    for (rid, aname), acorr in autocorr_CH3.items():
        outname = "%02d_%s.csv" % (rid.serial, aname,)
        pd.DataFrame(np.array([np.linspace(0, len(acorr) * 0.002, len(acorr), endpoint=False), acorr]).T,
                     columns=["time_ns", "acorr"]).to_csv(os.path.join(output_directory, outname), index=False)


def get_methyl_vectors(frame: Frame):
    """

    :param frame:
    :return: List[Tuple(Atom,Atom)]
    """
    CH3_dict = {('ALA', 'CB'): ['CB', 'HB1'],
                ('VAL', 'CG1'): ['CG1', 'HG11'],
                ('VAL', 'CG2'): ['CG2', 'HG21'],
                ('THR', 'CG2'): ['CG2', 'HG21'],
                ('LEU', 'CD1'): ['CD1', 'HD11'],
                ('LEU', 'CD2'): ['CD2', 'HD21'],
                ('ILE', 'CD1'): ['CD1', 'HD11'],
                ('ILE', 'CG2'): ['CG2', 'HG21'],
                ('MET', 'CE1'): ['CE1', 'HE1']}

    atom_pairs = []

    for r in frame.asResidues:
        for atom in r.asAtoms:
            if (r.name.str, atom.name.str) in CH3_dict.keys():
                C = r[AtomName(CH3_dict[(r.name.str, atom.name.str)][0])]
                H = r[AtomName(CH3_dict[(r.name.str, atom.name.str)][1])]
                atom_pairs.append((C, H))

    return atom_pairs


def get_NH_vectors(frame: Frame):
    """

    :param frame:
    :return: List[Tuple(Atom,Atom)]
    """

    atom_pairs = []

    for r in frame.asResidues:
        try:
            atom_pairs.append((r[AtomName("N")], r[AtomName("H")]))
        except pyxmolpp2.polymer.OutOfRangeResidue:
            pass

    return atom_pairs


def extract_inertia_tensor_vectors_autocorr(path_to_trajectory: str, output_directory: str) -> None:
    """

    :param path_to_trajectory:
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    ...
