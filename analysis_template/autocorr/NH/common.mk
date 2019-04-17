SCRIPT_DIR=/home/olebedenko/bioinf/scripts/md-timescales/md_timescales
N_RESIDUES=76
TRAJECTORY_LENGTH=1000
TRAJECTORY_PATH=/home/olebedenko/bioinf/trj/ubq/spce/NVE/
PYTHONPATH=${SCRIPT_DIR}/..
.EXPORT_ALL_VARIABLES:

.PHONY: all clean

all: figures/tau_2_exp.pdf figures/tau_4_exp.pdf

clean:
	rm -rf data 
	rm -rf figures
