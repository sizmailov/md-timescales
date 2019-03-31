include common.mk


figures/tau_CH3_2_exp.pdf figures/tau_CH3_4_exp.pdf: ${SCRIPT_DIR}/plot_CH3_autocorr.py data/fit/tau_CH3_2_exp.csv data/fit/tau_CH3_4_exp.csv
	mkdir -p figures
	python3 "$<" \
		--path_to_fit_csv=data/fit/ \
		--path_to_csv_acorr=data/acorr/ \
		--output_directory=figures \




data/fit/tau_CH3_2_exp.csv data/fit/tau_CH3_4_exp.csv: ${SCRIPT_DIR}/fit_CH3_autocorr.py data/acorr/${FIRST_CH3_RESIDUES}.csv data/acorr/${LAST_CH3_RESIDUES}.csv
	mkdir -p data/fit/
	python3 "$<" \
		--output_directory=data/fit/ \
		--path_to_csv_accor=data/acorr/ \
		--path_to_ref="${TRAJECTORY_PATH}"



data/acorr/${FIRST_CH3_RESIDUES}.csv data/acorr/${LAST_CH3_RESIDUES}.csv:  ${SCRIPT_DIR}/extract_CH3_autocorr.py
	mkdir -p data/acorr/
	python3 "$<" \
		--output_directory=data/acorr/ \
		--length_trajectory="${TRAJECTORY_LENGTH}" \
		--path_to_trajectory="${TRAJECTORY_PATH}"

