include common.mk


figures/fit_msd_plot.png: ${SCRIPT_DIR}/plot_msd.py data/msd.csv
	mkdir -p figures
	python3 "$<" \
		--path-to-msd=data/msd.csv \
		--output-directory=figures/ \



data/msd.csv: ${SCRIPT_DIR}/extract_msd.py 
	mkdir -p data/fit/
	python3 "$<" \
		--output-directory=data/ \
		--trajectory-length="${TRAJECTORY_LENGTH}" \
		--path-to-trajectory="${TRAJECTORY_PATH}" \
