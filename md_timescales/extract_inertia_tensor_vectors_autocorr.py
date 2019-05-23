from bionmr_utils.md import *
from tqdm import tqdm
import numpy as np
import os
import pandas as pd
import argparse
from md_timescales.extract import extract_time_step_ns


def extract_inertia_tensor_vectors_autocorr(path_to_trajectory: str, output_directory: str, trajectory_length=1):
    os.makedirs(output_directory, exist_ok=True)
    traj, ref = traj_from_dir(path_to_trajectory, first=1, last=trajectory_length)
    ref = traj[0]
    ref_ca = ref.asAtoms.filter(aName == "CA")
    frame_ca = None

    nodes = pd.read_csv("/home/sergei/GB1/nodes/225.txt.csv", names=["x","y","z","w"])
    xyz = nodes[["x","y","z"]].values
    weights = nodes["w"]

    vectors = [VectorXYZ() for v in xyz]

    for frame in tqdm(traj,desc="extract vectors"):
        if frame_ca is None:
            frame_ca = frame.asAtoms.filter(aName == "CA")

        al = ref_ca.alignment_to(frame_ca)
        m = al.matrix3d()
        for vector, node in zip(vectors, xyz):
            v1 = m.dot(node)
            vector.append(XYZ(v1[0], v1[1], v1[2]))

    sum_acorr = None
    for i, vector in enumerate(tqdm(vectors, desc="calc autocorr")):
        w = weights[i]
        autocorr = np.array(calc_autocorr_order_2(vector, 1000 * 100)) * w
        if sum_acorr is None:
            sum_acorr = autocorr
        else:
            sum_acorr += autocorr
    avg_acorr = sum_acorr / 4 / np.pi
    T = np.linspace(0, len(avg_acorr) * extract_time_step_ns(path_to_trajectory), len(avg_acorr))
    df = pd.DataFrame(np.array([T, avg_acorr]).T, columns=["time_ns", "acorr"])
    df.to_csv(os.path.join(output_directory, "tumbling_avg.csv"), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate_overall_tumbling')
    parser.add_argument('--path-to-trajectory', required=True, )
    parser.add_argument('--output-directory', default=os.getcwd())
    parser.add_argument('--trajectory-length', default=1, type=int)
    args = parser.parse_args()
    extract_inertia_tensor_vectors_autocorr(args.path_to_trajectory, args.output_directory, args.trajectory_length)
