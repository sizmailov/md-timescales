import os
import pyxmolpp2
import pandas as pd
import numpy as np
from tqdm import tqdm
from typing import *
from bionmr_utils.md import *


def extract_time_step_ns(path_to_trajectory: str) -> float:
    """
    Extract time step from trajectory

    :param path_to_trajectory: path to directory
    :return time between trajectory frames

    """
    path_to_firt_run_in = os.path.join(path_to_trajectory, "5_run", "run00001.in")
    with open(path_to_firt_run_in) as first_run_in:
        for line in first_run_in:
            row = line.strip().split()
            if row[0] == 'ntwx':
                ntwx = int(row[2].strip(","))
            if row[0] == 'dt':
                dt = float(row[2].strip(","))
    time_step_ns = dt / 1000 * ntwx
    return time_step_ns


def atom_name_mass_map(atom: Atom) -> float:
    """
    Mass map for H, C, O, N, P, S atoms

    :param atom 
    :return atomic weight

    """
    atomics_weight = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.973762, 'S': 32.06}
    return atomics_weight[atom.name.str[0]]


def extract_lattice_vectors_rst7(lattice_path: str) -> LatticeVectors:
    """
    Extract lattice vectors from Amber input coordinate (rst7) file

    :param lattice_path: path to input coordinate (rst7) file

    """
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

    :param traj: trajectory returns Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice]
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

        if volume is not None:
            factor_scale = (volume[frame.index] / reference_volume) ** (1 / 3)
            bsf.scale_lattice_by(factor_scale)

        current_cm = atoms.mass_center(masses)

        if prev_cm:
            current_cm += bsf.find_best_shift(prev_cm, current_cm)[1]

        time.append(frame.index * dt)
        mass_centers.append(current_cm)

        prev_cm = current_cm

        if volume is not None:
            bsf.scale_lattice_by(1 / factor_scale)
    return np.array(time), mass_centers


def get_autocorr(trajectory: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice],
                 ca_alignment: bool,
                 get_vectors: Callable[[Frame], List[Tuple[Atom, Atom]]]) -> Dict[Tuple, Float]:
    """
    Get autocorrelation from trajectory

    :param trajectory: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice]
    :param ca_alignment: condition of aligment frames by Ca atoms
    :param get_vectors: function to determinate tracked vector
    :return dict of (rid, aname): autocorrelation
    """
    ref = trajectory[0]
    ref_ca = ref.asAtoms.filter(aName == "CA")
    frame_ca = None
    CH3S = None
    vectors = None
    for frame in tqdm(trajectory):
        if ca_alignment:
            if frame_ca is None:
                frame_ca = frame.asAtoms.filter(aName == "CA")
                frame_ats.transform(frame_ca.aligment_to(ref_ca))
        if CH3S is None:
            CH3S = {}
            vectors = {}
            for A, B in get_vectors(new_frame):
                CH3S[A.rId, A.aName] = (A, B)
                vectors[A.rId, A.aName] = VectorXYZ()

        for (rid, aname), (A, B) in CH3S.items():
            vectors[rid, aname].append(A.r - B.r)

    autocorr = {
        (rid, aname): calc_autocorr_order_2(vector)
        for (rid, aname), vector in vectors.items()
    }
    return autocorr


def save_autocorr(autocorr: dict,
                  time_step_ns: float,
                  output_directory: str) -> None:
    """

    :param autocorr: returns dict of (rid, aname): autocorrelation
    :param time_step_ns: returns time between trajectory frames
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]

    """
    for (rid, aname), acorr in autocorr.items():
        outname = "%02d_%s.csv" % (rid.serial, aname,)
        pd.DataFrame(np.array([np.linspace(0, len(acorr) * time_step_ns, len(acorr), endpoint=False), acorr]).T,
                     columns=["time_ns", "acorr"]).to_csv(os.path.join(output_directory, outname), index=False)


def get_autocorr(trajectory: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice],
                  time_step_ns: float,
                  output_directory: str,
                  get_vectors: Callable[[Frame], List[Tuple[Atom, Atom]]],
                  ca_alignment: bool) -> None:
    """

    :param get_vectors: returns list of atom pairs of interest
    :param path_to_trajectory: path to directory
    :param trajectory_length: number of .dat files to process
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]

    """
    autocorr = get_autocorr(trajectory, ca_alignment, get_vectors=get_vectors)
    save_autocorr(autocorr, time_step_ns, output_directory)


def extract_autocorr(path_to_trajectory: str,
                     output_directory: str,
                     get_vectors: Callable[[Frame], List[Tuple[Atom, Atom]]],
                     trajectory_length: int,
                     ca_alignment: bool) -> None:
    """

    :param get_vectors: returns list of atom pairs of interest
    :param path_to_trajectory:
    :param trajectory_length: number of .dat files to process
    :param ca_alignment: condition of aligment frames by Ca atoms
    :param get_vectors: function to determinate tracked vector
    """
    trajectory, ref = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    time_step_ns = extract_time_step_ns(path_to_trajectory)
    get_autocorr(trajectory, time_step_ns, output_directory, get_vectors, ca_alignment=False)
    if ca_alignment:
        output_directory = os.path.join(output_directory, "ca_alignment")
        os.makedirs(output_directory, exist_ok=True)
        get_autocorr(trajectory, time_step_ns, output_directory, get_vectors, ca_alignment=True)


def get_methyl_vectors(frame: Frame) -> List[Tuple(Atom, Atom)]:
    """

    :param frame: Frame
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


def get_NH_vectors(frame: Frame) -> List[Tuple(Atom, Atom)]:
    """

    :param frame: Frame
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
