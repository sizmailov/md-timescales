import os
import pyxmolpp2
import pandas as pd
import numpy as np
from tqdm import tqdm

from typing import *
from bionmr_utils.md import *

def extract_mass_center(path_to_trajectory: str, output_filename: str) -> None:
    """

    :param path_to_trajectory:
    :param output_filename: name of output four-column .csv file [time_ns,x,y,z]
    """
    ...


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


def extract_inertia_tensor_vectors_autocorr(path_to_trajectory: str, output_directory: str) -> None:
    """

    :param path_to_trajectory:
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    ...
