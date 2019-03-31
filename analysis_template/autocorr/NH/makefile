include common.mk


figures/tau_2_exp.pdf figures/tau_4_exp.pdf: ${SCRIPT_DIR}/plot_NH_autocorr.py data/fit/tau_2_exp.csv data/fit/tau_4_exp.csv
	mkdir -p figures
	python3 "$<" \
		--path_to_fit_csv=data/fit/ \
		--path_to_csv_acorr=data/acorr/ \
		--output_directory=figures \



data/fit/tau_2_exp.csv data/fit/tau_4_exp.csv: ${SCRIPT_DIR}/fit_vector_acorr.py data/acorr/02.csv data/acorr/${N_RESIDUES}.csv
	mkdir -p data/fit/
	python3 "$<" \
		--output-directory=data/fit/ \
		--path-to-acorrs=data/acorr/ \
		--path-to-trajectory="${TRAJECTORY_PATH}" \
		--vectors-group=NH


data/acorr/02.csv data/acorr/${N_RESIDUES}.csv:  ${SCRIPT_DIR}/extract_autocorr.py
	mkdir -p data/acorr/
	python3 "$<" \
		--output-directory=data/acorr/ \
		--trajectory-length="${TRAJECTORY_LENGTH}" \
		--path-to-trajectory="${TRAJECTORY_PATH}" \
		--vectors-group=NH
