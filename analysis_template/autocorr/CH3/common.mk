SCRIPT_DIR=/home/olebedenko/bioinf/scripts/md-timescales/md_timescales/
N_RESIDUES=76
TRAJECTORY_LENGTH=10
TRAJECTORY_PATH=/home/olebedenko/bioinf/trj/ubq/spce/NVE/

.PHONY: all clean

all: data/fit/tau_NH_2_exp.csv data/fit/tau_NH_4_exp.csv

clean:
	rm -rf data 
	rm -rf figures