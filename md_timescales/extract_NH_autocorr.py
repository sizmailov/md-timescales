from typing import *
from bionmr_utils.md import  *
from tqdm import tqdm
import numpy as np
import os
import pandas as pd


def extract_NH_autocorr(trajectory: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice]) -> dict:
    """
    :param trajectory
    """
    rPRO = ResidueName('PRO')
    NHS = None
    vectors = None
    for frame in tqdm(trajectory):
        if NHS is None:
            NHS = {}
            vectors = {}
            for r in frame.asResidues[1:]:
                if r.name != rPRO:
                    H = r[AtomName("H")]
                    N = r[AtomName("N")]
                    NHS[r.id] = (H,N)
                    vectors[r.id] = VectorXYZ()
        for rid,(H,N) in NHS.items():
            vectors[rid].append(H.r-N.r)
    autocorr = {rid: calc_autocorr_order_2(vector) 
                for rid,vector in vectors.items()}
    return autocorr


def extract_NH_autocorr(path_to_trajectory: str, trajectory_length=1: int, output_directory: str) -> None:

    traj, ref  = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    """
    :param path_to_trajectory:
    :param trajectory_length: by default equal 1
    :param output_directory: output directory for two-column .csv files [time_ns, acorr]
    """
    autocorr_NH = get_autocorr_NH(traj)
    os.makedirs(output_directory, exist_ok=True)
    chain = ref.asChains[0]
    for rid, acorr in acorrs.items():
        res_serial = rid.serial
        outname = "%02d.csv"%(res_serial,)
        pd.DataFrame(np.array([ np.linspace(0,len(acorr)*0.002, len(acorr), endpoint=False), acorr]).T, 
                    columns=["time_ns", "acorr"]).to_csv(os.path.join(outdir, outname), index=False)

    
if __name__ == '__main__':
    path_to_trajectory = "/home/olebedenko/bioinf/1ubq/1ubq/"
    output_directory = "/home/olebedenko/bioinf/handling/h4/tip3p/NPT_gamma_ln_2/autocorr/NH/data"
    trajectory_length = 2000
    extract_NH_autocorr(path_to_trajectory, trajectory_length, output_directory)

