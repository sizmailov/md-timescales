import os


def extract_time_step_ns(path_to_trajectory):
	path_to_firt_run_in = os.path.join(path_to_trajectory, "5_run", "run00001.in")
	with open(path_to_firt_run_in) as first_run_in:
		for line in first_run_in:
			row = line.strip().split()
			if row[0] == 'nstlim':
				nstlim = int(row[2])
			if row[0] == 'ntpr':
				ntpr = int(row[2])
			if row[0] == 'dt':
				dt = float(row[2])
	time_step_ns = (dt * 1000) / (nstlim / ntpr)
	return(time_step_ns)
