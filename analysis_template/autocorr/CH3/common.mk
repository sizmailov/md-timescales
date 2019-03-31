SCRIPT_DIR=/home/sergei/PycharmProjects/md-timescales/md_timescales
N_RESIDUES=76
FIRST_CH3_RESIDUES = 03_CD1
LAST_CH3_RESIDUES = 73_CD2
TRAJECTORY_LENGTH=100
TRAJECTORY_PATH=/home/sergei/bionmr/olebedenko/bioinf/trj/ubq/spce/NVE/

.PHONY: all clean

all: figures/tau_2_exp.pdf figures/tau_4_exp.pdf

clean:
	rm -rf data 
	rm -rf figures