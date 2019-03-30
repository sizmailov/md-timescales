include common.mk

data/acorr/overall_tumbling1.csv data/acorr/overall_tumbling3.csv: ${SCRIPT_DIR}extract_inertia_tensor_vectors_autocorr.py
	mkdir -p data/acorr/
	python3 $< \
		--output_directory=data/acorr/ \
		--length_trajectory=${TRAJECTORY_LENGTH} \
		--path_to_trajectory=${TRAJECTORY_PATH}
