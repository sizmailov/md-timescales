from typing import *
from bionmr_utils.md import *
import argparse
import pyxmolpp2
from tqdm import tqdm
import numpy as np
import os
import pandas as pd


def get_methyls_vectors(frame: Frame):
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
            for A,B in get_vectors(frame):
                CH3S[A.rId, A.aName] = (A, B)
                vectors[A.rId, A.aName] = VectorXYZ()

        for (rid, aname), (A, B) in CH3S.items():
            vectors[rid, aname].append(A.r - B.r)

    autocorr = {
        (rid, aname): calc_autocorr_order_2(vector)
        for (rid, aname), vector in vectors.items()
    }

    return autocorr


def extract_CH3_autocorr(path_to_trajectory: str, output_directory: str, trajectory_length: int = 1) -> None:
    traj, ref = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    """
    :param path_to_trajectory:
    :param trajectory_length: by default equal 1
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    autocorr_CH3 = get_autocorr(traj, get_vectors=get_methyls_vectors)
    chain = ref.asChains[0]
    for (rid, Cname), acorr in autocorr_CH3.items():
        res_serial = rid.serial
        outname = "%02d_%s.csv" % (res_serial, Cname,)
        pd.DataFrame(np.array([np.linspace(0, len(acorr) * 0.002, len(acorr), endpoint=False), acorr]).T,
                     columns=["time_ns", "acorr"]).to_csv(os.path.join(output_directory, outname), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calc NH autocorr')
    parser.add_argument('-i', '--path_to_trajectory', required=True, )
    parser.add_argument('-o', '--output_directory', default=os.getcwd())
    parser.add_argument('-l', '--length_trajectory', default=1, type=int)
    args = parser.parse_args()
    extract_CH3_autocorr(args.path_to_trajectory, args.output_directory, args.length_trajectory)
