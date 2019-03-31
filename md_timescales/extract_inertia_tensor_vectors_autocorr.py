from bionmr_utils.md import *
from tqdm import tqdm
import numpy as np
import os
import pandas as pd
import argparse

def diagonalize(A):
    eigen_values, eigen_vectors = np.linalg.eig(A)
    idx = eigen_values.argsort()[::-1]
    eigen_vectors = eigen_vectors[:, idx]
    return eigen_vectors


def extract_inertia_tensor_vectors_autocorr(path_to_trajectory: str, output_directory: str, trajectory_length=1):
    os.makedirs(output_directory, exist_ok=True)
    traj, ref = traj_from_dir(path_to_trajectory,first=1,last=trajectory_length)
    ref = traj[0]
    ref_ca = ref.asAtoms.filter(aName == "CA")
    inertia = ref_ca.geom_inertia_tensor()
    eigen_vectors = diagonalize(inertia)
    frame_ca = None

    vectors1 = VectorXYZ()
    vectors2 = VectorXYZ()
    vectors3 = VectorXYZ()

    for frame in tqdm(traj):
        if frame_ca is None:
            frame_ca = frame.asAtoms.filter(aName == "CA")

        al = ref_ca.alignment_to(frame_ca)
        m = al.matrix3d()
        v1 = m.dot(eigen_vectors[0])
        vectors1.append(XYZ(v1[0], v1[1], v1[2]))
        v2 = m.dot(eigen_vectors[1])
        vectors2.append(XYZ(v2[0], v2[1], v2[2]))
        v3 = m.dot(eigen_vectors[2])
        vectors3.append(XYZ(v3[0], v3[1], v3[2]))

    autocorr1 = calc_autocorr_order_2(vectors1)
    autocorr2 = calc_autocorr_order_2(vectors2)
    autocorr3 = calc_autocorr_order_2(vectors3)
    T = np.linspace(0,len(autocorr1)*0.002, len(autocorr1))

    pd.DataFrame(autocorr1, columns=["C"]).to_csv(os.path.join(output_directory, "overall_tumbling_1.csv"), index = False)
    pd.DataFrame(autocorr2, columns=["C"]).to_csv(os.path.join(output_directory, "overall_tumbling_2.csv"), index = False)
    pd.DataFrame(autocorr3, columns=["C"]).to_csv(os.path.join(output_directory, "overall_tumbling_3.csv"), index = False)
    pd.DataFrame(np.array([T, autocorr1]).T, columns=["time_ns", "acorr"]).to_csv(os.path.join(output_directory, "overall_tumbling_1.csv"), index = False)
    pd.DataFrame(np.array([T, autocorr2]).T, columns=["time_ns", "acorr"]).to_csv(os.path.join(output_directory, "overall_tumbling_2.csv"), index = False)
    pd.DataFrame(np.array([T, autocorr3]).T, columns=["time_ns", "acorr"]).to_csv(os.path.join(output_directory, "overall_tumbling_3.csv"), index = False)

if __name__ == '__main__':
    # path_to_trajectory = "/home/sergei/UBI/case_wong_replica/1ubq_intolerant_shake_ewald_SPCE/"
    # output_directory = "data"
    parser = argparse.ArgumentParser(description='Calculate_overall_tumbling')
    parser.add_argument('-i', '--path_to_trajectory', required=True, )
    parser.add_argument('-o', '--output_directory', default=os.getcwd())
    parser.add_argument('-l', '--trajectory_length', default=1, type=int)
    args = parser.parse_args()
    extract_inertia_tensor_vectors_autocorr(args.path_to_trajectory, args.output_directory, args.trajectory_length)

