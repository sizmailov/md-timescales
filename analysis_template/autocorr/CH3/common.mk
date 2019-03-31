SCRIPT_DIR=/home/olebedenko/bioinf/scripts/md-timescales/md_timescales/
N_RESIDUES=76
FIRST_CH3_RESIDUES = 03_HD11
LAST_CH3_RESIDUES = 73_HD21
TRAJECTORY_LENGTH=10
TRAJECTORY_PATH=/home/olebedenko/bioinf/trj/ubq/spce/NVE/

.PHONY: all clean

all: figures/tau_CH3_2_exp.pdf figures/tau_CH3_4_exp.pdf

clean:
	rm -rf data 
	rm -rf figures