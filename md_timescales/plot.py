import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os
from matplotlib.backends.backend_pdf import PdfPages
from typing import *
from .fit import __multi_exp_f


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


def add_relpath_to_top_corner(figure: plt.Figure):
    big_axis = figure.add_axes([0, 0, 1, 1], facecolor=(1, 1, 1, 0))

    big_axis.text(0.99, 0.99,
                  os.path.relpath(os.path.abspath(os.curdir), os.path.expanduser("~/bioinf/handling")),
                  color="#CCCCCC",
                  horizontalalignment='right',
                  verticalalignment='top')

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
    ax.set_xlabel('time, ns', fontsize=13)
    ax.set_ylabel('C(t)', fontsize=13)
    ax.grid(True)
    return fig, ax


def plot_figure_autocorr(time, fit_line, acorr):
    amplitude = fit_line.filter(like='-a')
    tau = fit_line.filter(like='-tau')
    if fit_line['aName'][0] == "C":
        constant = fit_line['constant']
    else:
        constant = 0
    order = len(amplitude)
    rid = fit_line["rId"]
    rname = fit_line["rName"]
    aname = fit_line["aName"]
    limit = fit_line["limit"]

    graph_label = __get_autocorr_graph_label(fit_line)
    fig, ax = settings_plot(graph_label)
    if aname == "N":
        ax.set_title('NH autocorrelation plot %d exp %s %s' % (order, rid, rname))
        ax.set_xlim(-1, 20)
    elif aname[0] == "C":
        ax.set_title('CH3 autocorrelation plot %d exp %s %s %s' % (order, rid, rname, aname))
        ax.set_xlim(-0.05, 2)
    ax.set_ylim(-0.1, 1.1)
    ax.plot(time, acorr)
    ax.plot(time, __multi_exp_f(time, amplitude,
                                tau, C=constant))
    ax.axvline(x=time[limit], color='g', linestyle='--')
    return fig, ax


def get_plot_acorr_fit(path_to_fit_csv: str,
                       path_to_csv_acorr: str,
                       output_directory: str) -> None:
    """
    Plot one pdf with a particular fit function (e.g. tau-2-exp.pdf)

    :param path_to_fit_csv: path to particular fit
           function values (e.g. tau-2-exp.csv) files
    :param path_to_csv_acorr: path to two-column .csv files [time_ns, acorr]
    :param output_directory: path to fit pdf with
           a particular fit function (e.g. tau-2-exp.pdf)
    """
    exp_order = {2: "tau_2_exp", 3: "tau_3_exp", 4: "tau_4_exp"}
    for order in range(2, 5):
        with PdfPages(os.path.join(output_directory, exp_order[order] + ".pdf")) as pdf:
            csv_fit = os.path.join(path_to_fit_csv, exp_order[order] + ".csv")
            fit = pd.read_csv(csv_fit)
            for _, fit_line in fit.iterrows():
                file = "{}/{:02d}_{}.csv".format(path_to_csv_acorr, fit_line["rId"], fit_line["aName"])
                df = pd.read_csv(file)
                fig, ax = plot_figure_autocorr(df.time_ns, fit_line, df.acorr)
                add_relpath_to_top_corner(fig)
                pdf.savefig(fig)
                plt.close(fig)


def plot_acorr_fit(path_to_fit_csv: str,
                   path_to_csv_acorr: str,
                   output_directory: str,
                   ca_alignment: bool = False) -> None:
    """

    :param path_to_fit_csv: path to particular fit
           function values (e.g. tau-2-exp.csv) files
    :param path_to_csv_acorr: path to two-column .csv files [time_ns, acorr]
    :param output_directory: path to fit pdf with
           a particular fit function (e.g. tau-2-exp.pdf)
    :param ca_alignment: flag of aligment frames by Ca atoms
    """
    get_plot_acorr_fit(path_to_fit_csv, path_to_csv_acorr, output_directory)
    if ca_alignment:
        path_to_fit_csv = os.path.join(path_to_fit_csv, "ca_alignment")
        path_to_csv_acorr = os.path.join(path_to_csv_acorr, "ca_alignment")
        output_directory = os.path.join(output_directory, "ca_alignment")
        os.makedirs(output_directory, exist_ok=True)
        get_plot_acorr_fit(path_to_fit_csv, path_to_csv_acorr, output_directory)
