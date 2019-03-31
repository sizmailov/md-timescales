include common.mk


figures/tau_NH_2_exp.pdf figures/tau_NH_4_exp.pdf: ${SCRIPT_DIR}/plot_NH_autocorr.py data/fit/tau_NH_2_exp.csv data/fit/tau_NH_4_exp.csv
	mkdir -p figures
	python3 "$<" \
		--path_to_fit_csv=data/fit/ \
		--path_to_csv_acorr=data/acorr/ \
		--output_directory=figures \



data/fit/tau_NH_2_exp.csv data/fit/tau_NH_4_exp.csv: ${SCRIPT_DIR}/fit_NH_autocotrr.py data/acorr/02.csv data/acorr/${N_RESIDUES}.csv
	mkdir -p data/fit/
	python3 "$<" \
		--output_directory=data/fit/ \
		--path_to_csv_accor=data/acorr/ \
		--path_to_ref="${TRAJECTORY_PATH}"


data/acorr/02.csv data/acorr/${N_RESIDUES}.csv:  ${SCRIPT_DIR}/extract_NH_autocorr.py
	mkdir -p data/acorr/
	python3 "$<" \
		--output_directory=data/acorr/ \
		--length_trajectory="${TRAJECTORY_LENGTH}" \
		--path_to_trajectory="${TRAJECTORY_PATH}"


