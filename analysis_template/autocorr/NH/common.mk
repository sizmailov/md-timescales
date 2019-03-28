SCRIPT_DIR=/home/olebedenko/bioinf/scripts/md-timescales/md_timescales/
N_RESIDUES=76
TRAJECTORY_LENGTH=1000
TRAJECTORY_PATH=/home/olebedenko/bioinf/trj/ubq/spce/NVE/

.PHONY: all clean

all: graph/tau_NH_2_exp.pdf graph/tau_NH_4_exp.pdf

clean:
	rm -rf data 
	rm -rf graph
