import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os
from matplotlib.backends.backend_pdf import PdfPages
from fit import __multi_exp_f

def __get_autocorr_graph_label(fit_line):
    amplitude = fit_line.filter(like='-a')
    tau = fit_line.filter(like='-tau')
    union_a_tau = ["{a_label:2s} = {a_value:5.3f} ; {tau_label:3s} = {tau_value:8.3e}".format(
        a_label=a_label,
        a_value=fit_line[a_label],
        tau_label=tau_label,
        tau_value=fit_line[tau_label])
        for a_label, tau_label in zip(amplitude.index.tolist(), tau.index.tolist())
    ]
    if fit_line['aName'][0] == "C":
        union_a_tau.append(" ; {constant_label:3s} = {constant:5.3f}".format(
            constant=fit_line['constant'],
            constant_label="C"))
    graph_label = "\n".join(union_a_tau)
    return graph_label


def settings_plot(graph_label):
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
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlim(-1, 20)
    ax.set_xlabel('time, ns', fontsize=13)
    ax.set_ylabel('C(t)', fontsize=13)
    ax.grid(True)
    return fig, ax


def plot_figure_autocorr(time, fit_line, acorr):

    amplitude = fit_line.filter(like='-a')
    tau = fit_line.filter(like='-tau')

    order = len(amplitude)
    rid = fit_line["rId"]
    rname = fit_line["rName"] 
    aname = fit_line["aName"]
    limit = fit_line["limit"]

    graph_label = __get_autocorr_graph_label(fit_line)
    fig, ax = settings_plot(graph_label)
    if aname == "N":
        ax.set_title('NH autocorrelation plot %d exp %s %s' % (order, rid, rname))
    elif aname[0] == "C":
        ax.set_title('CH3 autocorrelation plot %d exp %s %s %s' % (order, rid, rname, aname))
    ax.plot(time, acorr)
    ax.plot(time, __multi_exp_f(time, amplitude.values.flatten(),
                                tau.values.flatten(), C=0))
    ax.axvline(x=time[limit], color='g', linestyle='--')
    return fig, ax


def plot_acorr_fit(path_to_fit_csv, path_to_csv_acorr, output_directory):
    """
    Plot one pdf with a particular fit function (e.g. NH-acorr-2-exp.pdf)
     - [ ] Show acorr data
     - [ ] Show fit curve
     - [ ] Show fit parameters
     - [ ] Denote fit region by vertical line
    """
    exp_order = {2: "tau_2_exp", 3: "tau_3_exp", 4: "tau_4_exp"}
    csv_files = sorted(glob.glob(os.path.join(path_to_csv_acorr, "*.csv")))
    for order in range(2, 5):
        with PdfPages(os.path.join(output_directory, exp_order[order] + ".pdf")) as pdf:
            csv_fit = os.path.join(path_to_fit_csv, exp_order[order] + ".csv")
            fit = pd.read_csv(csv_fit)
            for ind, file in enumerate(csv_files):
                df = pd.read_csv(file)
                fit_line = fit.iloc[ind]
                fig, ax = plot_figure_autocorr(df.time_ns, fit_line, df.acorr)
                pdf.savefig(fig)
                plt.close(fig)

