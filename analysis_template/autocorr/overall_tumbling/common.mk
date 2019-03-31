SCRIPT_DIR=/home/legosta/bioinf/scripts/md-timescales/md_timescales/
N_RESIDUES=76
TRAJECTORY_LENGTH=1000
TRAJECTORY_PATH=/home/legosta/bioinf/trj/ubq/spce/NVE/

.PHONY: all clean 

all: figures/tau_inertia_tensor_2_exp.pdf figures/tau_inertia_tensor_4_exp.pdf

clean:
	rm -rf data 
	rm -rf figures



