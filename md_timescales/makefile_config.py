import argparse
import os
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print values for makefile')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--n-residues', default=None, action='store_true')
    group.add_argument('--trajectory-length', default=None, action='store_true')
    group.add_argument('--trajectory-path', default=None, action='store_true')
    group.add_argument('--first-ch3-residue', default=None, action='store_true')
    group.add_argument('--last-ch3-residue', default=None, action='store_true')
    group.add_argument('--volume-file', default=None, action='store_true')
    group.add_argument('--lattice-rst7-file', default=None, action='store_true')

    args = parser.parse_args()

    for folder in os.path.abspath(os.curdir).split("/"):

        if folder == "ubq":
            n_residues = 76
            first_ch3_residue = "03_CD1"
            last_ch3_residue = "73_CD2"

        if folder == "h4":
            n_residues = 25
            first_ch3_residue = "10_CD1"
            last_ch3_residue = "22_CD2"

        if folder.startswith("NVE") or folder.startswith("NVT"):
            lattice_rst7_file = "data/run00001.rst7"
            volume_file = "None"

        if folder.startswith("NPT"):
            lattice_rst7_file = "None"
            volume_file = "data/stats/summary.VOLUME"

    path_fragments = os.path.abspath(os.curdir).split("/")
    for i in range(1, len(path_fragments)):
        p = os.path.join("/", *path_fragments[:-i]).replace('/handling/', '/trj/')
        if os.path.exists(p):
            trajectory_path = p
            break

    trajectory_length = len(glob.glob(os.path.join(trajectory_path,"5_run","*.dat")))

    if args.n_residues:
        print(n_residues)
    if args.trajectory_path:
        print(trajectory_path)
    if args.trajectory_length:
        print(trajectory_length)
    if args.volume_file:
        print(volume_file)
    if args.lattice_rst7_file:
        print(lattice_rst7_file)
    if args.last_ch3_residue:
        print(last_ch3_residue)
    if args.first_ch3_residue:
        print(first_ch3_residue)
