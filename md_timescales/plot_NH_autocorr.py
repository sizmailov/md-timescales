import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnchoredText




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
        pdf = PdfPages(os.path.join(path_to_output_pdf, exp_order[order] + ".pdf"))
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
            plt.figure()
            fig, ax = plt.subplots()

            ax.text(1, 3, 'boxed italics text in data coords', style='italic',
                 bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})

            fig.tight_layout()
            plt.ylim(-0.3, 1.5)
            plt.xlim(0, 20)
            plt.xlabel('time, ns', fontsize = 13)
            plt.ylabel('autocorrelation', fontsize = 13)
            plt.title('NH autocorrelation plot %s exp %s%s'%(str(order), fit_line["rId"], fit_line["rName"]))
            plt.plot( df.time_ns, df.acorr, linestyle=':')
            plt.plot(df.time_ns, fit_func[order](df.time_ns, **popt))
            #plt.axvline(x=df.time_ns[limit], color='g', linestyle='--', label="fit limit %s"%(limit))
            plt.grid(True)
            #popt = [round(elem,2) for elem in popt]
            # plt.legend(("NH autocorrelation", "fit {0}*exp-(-t/{1}) + {2}*exp-(-t/{3}) + \n +{4}*exp-(-t/{5}) + {6}*exp-(-t/{7})".format(*popt)
            #             , "fit limit = %s ns"%(df.time_ns[limit])), loc='upper right')
            pdf.savefig()
            plt.close()
        pdf.close()

"""
test = "/home/olebedenko/bioinf/scripts/md-timescales/md_timescales/tau_NH_2_exp.csv"
fit = pd.read_csv(test)
line = fit.iloc[23]
a = ["exp-%d-a%d" %  (2, i + 1) for i in range(2)]
t = ["exp-%d-tau%d" %  (2, i + 1) for i in range(2)]
amp = {"a%d" % (i+1): line[a[i]] for i in range(2)}
tau = {"tau%d" % (i+1): line[t[i]] for i in range(2)}
amp.update(tau)

def r(a1, tau1, a2, tau2):
    print(a1,tau1,a2,tau2)
r(**amp)
print(line)

path_to_csv_accor = "/home/olebedenko/bioinf/handling/h4/tip4p-ew/NPT_gamma_ln_2/autocorr/NH/data"
csv_files = sorted(glob.glob(os.path.join(path_to_csv_accor, "*.csv")))
"""
if __name__ == "__main__":
    #-f_csv "/home/olebedenko/bioinf/handling/h4/tip4p-ew/NPT_gamma_ln_2/autocorr/NH/graph" -a_csv "/home/olebedenko/bioinf/handling/h4/tip4p-ew/NPT_gamma_ln_2/autocorr/NH/data"
    parser = argparse.ArgumentParser(description="plot NH autocorrelation")
    parser.add_argument('-f_csv', '--path_to_fit_csv', required=True)
    parser.add_argument('-a_csv', '--path_to_accor_csv', required=True)
    parser.add_argument('-o', '--path_to_output_pdf', default=os.getcwd())
    args = parser.parse_args()
    plot_acorr_fit(args.path_to_fit_csv, args.path_to_accor_csv, args.path_to_output_pdf)
