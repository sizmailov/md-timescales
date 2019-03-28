include common.mk


graph/tau_NH_2_exp.pdf graph/tau_NH_4_exp.pdf: ${SCRIPT_DIR}plot_NH_autocotrr.py graph/tau_NH_2_exp.csv ${SCRIPT_DIR}plot_NH_autocotrr.py graph/tau_NH_4_exp.csv
	mkdir -p graph
	python3 $< \
		--path_to_fit_csv=graph \
		--path_to_csv_acorr=data \
		-output_directory=graph 



graph/tau_NH_2_exp.csv graph/tau_NH_4_exp.csv: ${SCRIPT_DIR}fit_NH_autocotrr.py data/02.csv data/${N_RESIDUES}.csv
	mkdir -p graph
	python3 $< \
		--output_directory=graph \
		--path_to_csv_accor=data \
		--path_to_ref=${TRAJECTORY_PATH}


data/02.csv data/${N_RESIDUES}.csv:  ${SCRIPT_DIR}extract_NH_autocorr.py
	mkdir -p data
	python3 $< \
		--output_directory=data \
		--length_trajectory=${TRAJECTORY_LENGTH} \
		--path_to_trajectory=${TRAJECTORY_PATH}


