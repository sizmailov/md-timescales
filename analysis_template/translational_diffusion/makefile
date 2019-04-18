SCRIPT_DIR:=${HOME}/bioinf/md-timescales/md_timescales
include ${SCRIPT_DIR}/../analysis_template/translational_diffusion/common.mk



figures/fit_msd_plot.png: ${SCRIPT_DIR}/plot_msd.py data/msd.csv
	mkdir -p figures
	python3 "$<" \
		--path-to-msd=data/msd.csv \
		--output-directory=figures/ \


data/stats/summary.VOLUME: ${TRAJECTORY_PATH}/5_run/*.out
	mkdir -p data/stats ;
	(cd data/stats; process_mdout.perl ${TRAJECTORY_PATH}/5_run/*.out)

data/run00001.rst7:
	mkdir -p data
	cpptraj \
	    -p "${TRAJECTORY_PATH}/1_build/box.prmtop" \
	    -y "${TRAJECTORY_PATH}/5_run/run00001.rst" \
	    -x "$@"

data/cm.csv: ${SCRIPT_DIR}/extract_mass_center.py ${VOLUME_OR_RST7}
	mkdir -p data/
	python3 "$<" \
		--output-directory=data/ \
		--trajectory-length="${TRAJECTORY_LENGTH}" \
		--path-to-trajectory="${TRAJECTORY_PATH}" \
		--lattice-rst7-file=data/run00001.rst7


data/msd.csv: ${SCRIPT_DIR}/calc_msd.py data/cm.csv
	mkdir -p data/
	python3 "$<" \
		--mass-center-csv=data/cm.csv \
		--output-directory=data
