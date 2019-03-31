include common.mk

data/fit/tau_NH_2_exp.csv data/fit/tau_NH_4_exp.csv: "${SCRIPT_DIR}"fit_NH_autocotrr.py data/acorr/02.csv data/acorr/"${N_RESIDUES}".csv
	mkdir -p data/fit/
	python3 $< \
		--output_directory=data/fit/ \
		--path_to_csv_accor=data/acorr/ \
		--path_to_ref="${TRAJECTORY_PATH}"
