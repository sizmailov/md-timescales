SCRIPT_DIR=/home/legosta/bioinf/scripts/md-timescales/md_timescales/
N_RESIDUES=76
TRAJECTORY_LENGTH=10
TRAJECTORY_PATH=/home/legosta/bioinf/trj/ubq/spce/NVE/

.PHONY: all clean 

all: data/acorr/overall_tumbling1.csv data/acorr/overall_tumbling3.csv

clean:
	rm -rf data 
	rm -rf figures



