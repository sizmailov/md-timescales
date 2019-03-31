include common.mk


# figures/tau_inertia_tensor_2_exp.pdf figures/tau_inertia_tensor_4_exp.pdf: ${SCRIPT_DIR}plot_inertia_tensor_autocorr.py data/fit/tau_inertia_tensor_2_exp.csv data/fit/tau_inertia_tensor_4_exp.csv
# 	mkdir -p figures
# 	python3 $< \
# 		--path_to_fit_csv=data/fit/ \
# 		--path_to_csv_acorr=data/acorr/ \
# 		--output_directory=figures \



data/fit/tau_inertia_tensor_1_exp.csv: ${SCRIPT_DIR}mono_fit.py data/acorr/overall_tumbling_1.csv data/acorr/overall_tumbling_3.csv
	mkdir -p data/fit/
	python3 $< \
		--output_directory=data/fit/ \
		--path_to_csv_acorr=data/acorr/ \


data/acorr/overall_tumbling_1.csv data/acorr/overall_tumbling_3.csv: ${SCRIPT_DIR}extract_inertia_tensor_vectors_autocorr.py
	mkdir -p data/acorr/
	python3 $< \
		--output_directory=data/acorr/ \
		--length_trajectory=${TRAJECTORY_LENGTH} \
		--path_to_trajectory=${TRAJECTORY_PATH}
