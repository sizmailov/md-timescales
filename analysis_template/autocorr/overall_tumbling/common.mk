SCRIPT_DIR=/home/sergei/PycharmProjects/md-timescales/md_timescales/
N_RESIDUES=76
TRAJECTORY_LENGTH=100
TRAJECTORY_PATH=/home/sergei/bionmr/olebedenko/bioinf/trj/ubq/spce/NVE/
PYTHONPATH=${SCRIPT_DIR}/..

.PHONY: all clean
.EXPORT_ALL_VARIABLES:

all: figures/overall_tumbling.pdf

clean:
	rm -rf data 
	rm -rf figures



