import argparse

from md_timescales.fit import save_fit_auto_correlation

if __name__ == '__main__':

  #-ref  "/home/olebedenko/bioinf/1ubq/1ubq/" -accor "/home/olebedenko/bioinf/scripts/md-timescales/md_timescales"
  parser = argparse.ArgumentParser(description='fit CH3 autocorr')
  parser.add_argument('-ref', '--path_to_ref', required=True,)
  parser.add_argument('-accor', '--path_to_csv_accor', required=True)
  parser.add_argument('-o', '--output_directory', default="./")
  args = parser.parse_args()

  NH_bounds = [
        ([[0, 0.1, 0, 1], [1, 1, 1, 10]]),
        ([[0, 0.01, 0, 0.1, 0, 1], [1, 0.1, 1, 1, 1, 10]]),
        ([[0, 0.001, 0, 0.01, 0, 1, 0, 10], [1, 0.01, 1, 0.1, 1, 10, 1, 100]]),
  ]
  save_fit_auto_correlation(args.path_to_ref,
                            args.path_to_csv_accor,
                            args.output_directory,
                            curve_bounds=NH_bounds
                            )
