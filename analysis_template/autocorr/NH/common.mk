include ${SCRIPT_DIR}/../analysis_template/common.mk

all: figures/tau_2_exp.pdf figures/tau_4_exp.pdf

clean:
	rm -rf data 
	rm -rf figures
