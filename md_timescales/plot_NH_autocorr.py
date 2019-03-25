import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os
from matplotlib.backends.backend_pdf import PdfPages


def two_exp(x, a1, tau1, a2, tau2):
    return a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2)


def three_exp(x, a1, tau1, a2, tau2, a3, tau3):
    return a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2) + a3 * np.exp(-x / tau3)


def four_exp(x, a1, tau1, a2, tau2, a3, tau3, a4, tau4):
    return a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2) + a3 * np.exp(-x / tau3) + a4 * np.exp(-x / tau4)


def plot_acorr_fit(path_to_fit_csv, path_to_accor_csv, path_to_output_pdf):
    """
    Plot one pdf with a particular fit function (e.g. NH-acorr-2-exp.pdf)
     - [ ] Show acorr data
     - [ ] Show fit curve
     - [ ] Show fit parameters
     - [ ] Denote fit region by vertical line
    """
    fit_func = {2: two_exp, 3: three_exp, 4: four_exp}
    exp_order = {2: "tau_NH_2_exp", 3: "tau_NH_3_exp", 4: "tau_NH_4_exp"}
    csv_files = sorted(glob.glob(os.path.join(path_to_accor_csv, "*.csv")))
    for order in range(2,5):
        with PdfPages(os.path.join(path_to_output_pdf, exp_order[order] + ".pdf")) as pdf:
            csv_fit =  os.path.join(path_to_fit_csv, exp_order[order] + ".csv")
            fit = pd.read_csv(csv_fit)
            for ind, file in enumerate(csv_files):
                df = pd.read_csv(file)
                fit_line = fit.iloc[ind]
                amplitude_label = ["exp-%d-a%d" %  (order, i + 1) for i in range(order)]
                tau_label = ["exp-%d-tau%d" %  (order, i + 1) for i in range(order)]
                popt = {"a%d" % (i+1): fit_line[amplitude_label[i]] for i in range(order)}
                tau = {"tau%d" % (i+1): fit_line[tau_label[i]] for i in range(order)}
                popt.update(tau)

                coeff_a = ["a%d" % (i+1) for i in range(order)]
                coeff_tau = ["tau%d" % (i+1) for i in range(order)]
                union_a_tau = list(zip(coeff_a,
                               [" = " + str(round(popt[label],2)) + " " for label in coeff_a],
                               coeff_tau,
                               [" = " + str(round(popt[label],2)) + " ns\n" for label in coeff_tau]))
                graph_label = "".join(["".join(elem) for elem in union_a_tau])
                left, width = .40, .54
                bottom, height = .40, .54
                right = left + width
                top = bottom + height

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.text(right, top, graph_label,
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    multialignment='left',
                    bbox={'facecolor': 'moccasin', 'alpha': 0.5, 'pad': 6})
                ax.set_ylim(-0.1,1.1)
                ax.set_xlim(-1, 20)
                ax.set_xlabel('time, ns', fontsize = 13)
                ax.set_ylabel('autocorrelation', fontsize = 13)
                ax.set_title('NH autocorrelation plot %s exp %s%s'%(str(order), fit_line["rId"], fit_line["rName"]))
                ax.plot( df.time_ns, df.acorr)
                ax.plot(df.time_ns, fit_func[order](df.time_ns, **popt))
                #ax.set_axvline(x=df.time_ns[limit], color='g', linestyle='--', label="fit limit %s"%(limit))
                ax.grid(True)
                # ax.set_legend(("NH autocorrelation", "fit , loc='upper right')
                pdf.savefig()
                plt.close()

if __name__ == "__main__":
    #-f_csv "/home/olebedenko/bioinf/handling/h4/tip4p-ew/NPT_gamma_ln_2/autocorr/NH/graph" -a_csv "/home/olebedenko/bioinf/handling/h4/tip4p-ew/NPT_gamma_ln_2/autocorr/NH/data"
    parser = argparse.ArgumentParser(description="plot NH autocorrelation")
    parser.add_argument('-f_csv', '--path_to_fit_csv', required=True)
    parser.add_argument('-a_csv', '--path_to_accor_csv', required=True)
    parser.add_argument('-o', '--path_to_output_pdf', default=os.getcwd())
    args = parser.parse_args()
    plot_acorr_fit(args.path_to_fit_csv, args.path_to_accor_csv, args.path_to_output_pdf)
