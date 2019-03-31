include common.mk


figures/overall_tumbling.pdf: ${SCRIPT_DIR}/plot_overall_tumbling.py data/fit/tau_inertia_tensor_1_exp.csv
	mkdir -p figures
	python3 $< \
		--path_to_fit_csv=data/fit/ \
		--path_to_csv_acorr=data/acorr/ \
		--output_directory=figures \



data/fit/tau_inertia_tensor_1_exp.csv: ${SCRIPT_DIR}/mono_fit.py data/acorr/overall_tumbling_1.csv data/acorr/overall_tumbling_3.csv
	mkdir -p data/fit/
	python3 $< \
		--output_directory=data/fit/ \
		--path_to_csv_acorr=data/acorr/ \


data/acorr/overall_tumbling_1.csv data/acorr/overall_tumbling_3.csv: ${SCRIPT_DIR}/extract_inertia_tensor_vectors_autocorr.py
	mkdir -p data/acorr/
	python3 $< \
		--output_directory=data/acorr/ \
		--trajectory_length=${TRAJECTORY_LENGTH} \
		--path_to_trajectory="${TRAJECTORY_PATH}"
