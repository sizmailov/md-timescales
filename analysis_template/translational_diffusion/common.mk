include ${SCRIPT_DIR}/../analysis_template/common.mk

all: figures/fit_msd_plot.png

clean:
	rm -rf data 
	rm -rf figures