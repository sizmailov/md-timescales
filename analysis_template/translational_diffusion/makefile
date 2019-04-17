include common.mk


figures/fit_msd_plot.png: ${SCRIPT_DIR}/plot_msd.py data/msd.csv
	mkdir -p figures
	python3 "$<" \
		--path-to-msd=data/msd.csv \
		--output-directory=figures/ \



data/cm.csv: ${SCRIPT_DIR}/extract_mass_center.py 
	mkdir -p data/
	python3 "$<" \
		--output-directory=data/ \
		--trajectory-length="${TRAJECTORY_LENGTH}" \
		--path-to-trajectory="${TRAJECTORY_PATH}" \
		--lattice-rst7-file=run00001.rst7


data/msd.csv: ${SCRIPT_DIR}/calc_msd.py data/cm.csv
	mkdir -p data/
	python3 "$<" \
		--mass-center-csv=data/cm.csv \
		--output-directory=data
