from typing import *
from bionmr_utils.md import  *
import argparse
import pyxmolpp2
from tqdm import tqdm
import numpy as np
import os
import pandas as pd


def get_autocorr_CH3(trajectory: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice]) -> dict:
    """
    :param trajectory
    """
    CH3_dict = {('ALA','CB'): ['CB', 'HB1'],
         ('VAL','CG1'): ['CG1', 'HG11'],
         ('VAL','CG2'): ['CG2', 'HG21'],
         ('THR','CG2'): ['CG2', 'HG21'],
         ('LEU','CD1'): ['CD1', 'HD11'],
         ('LEU','CD2'): ['CD2', 'HD21'],
         ('ILE','CD1'): ['CD1', 'HD11'],
         ('ILE','CG2'): ['CG2', 'HG21'],
         ('MET','CE1'): ['CE1', 'HE1']}
    CH3S = None
    vectors = None
    for frame in tqdm(trajectory):
        if CH3S is None:
            CH3S = {}
            vectors = {}
            for r in frame.asResidues[1:]:
                for atom in r.asAtoms:
                    if (r.name.str, atom.name.str) in CH3_dict.keys():
                        C = r[AtomName(CH3_dict[(r.name.str, atom.name.str)][0])]
                        H = r[AtomName(CH3_dict[(r.name.str, atom.name.str)][1])]
                        CH3S[r.id, C.id] = (C, H)
                        vectors[r.id, C.id] = VectorXYZ()

        for (rid,Cid),(C,H) in CH3S.items():
            vectors[rid, Cid].append(C.r - H.r)

    autocorr = {
        (rid,Cid): calc_autocorr_order_2(vector)
        for (rid,Cid), vector in vectors.items()
    }

    return autocorr


def extract_CH3_autocorr(path_to_trajectory: str, output_directory: str, trajectory_length: int=1) -> None:

    traj, ref  = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    """
    :param path_to_trajectory:
    :param trajectory_length: by default equal 1
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    autocorr_CH3 = get_autocorr_CH3(traj)
    chain = ref.asChains[0]
    for (rid,Cid), acorr in autocorr_CH3.items():
        res_serial = rid.serial
        aName = chain.asAtoms[Cid].aName.str
        outname = "%02d_%s.csv"%(res_serial, aName,)
        pd.DataFrame(np.array([ np.linspace(0,len(acorr)*0.002, len(acorr), endpoint=False), acorr]).T, 
                    columns=["time_ns", "acorr"]).to_csv(os.path.join(output_directory, outname), index=False)

    
if __name__ == '__main__':

    #-i "/home/olebedenko/bioinf/1ubq/1ubq/" -l 20
    parser = argparse.ArgumentParser(description='Calc NH autocorr')
    parser.add_argument('-i', '--path_to_trajectory', required=True,)
    parser.add_argument('-o', '--output_directory',default=os.getcwd())
    parser.add_argument('-l', '--length_trajectory', default=1, type=int)
    args = parser.parse_args()
    extract_CH3_autocorr(args.path_to_trajectory, args.output_directory, args.length_trajectory)


