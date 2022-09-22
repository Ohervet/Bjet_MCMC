"""
File name: testing.py

Informal testing code, NOT unit tests
"""
import pickle
from datetime import datetime

import emcee
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interpolate

import blazar_model
import blazar_plots
import blazar_report
import blazar_utils
from blazar_properties import BASE_PATH

discard = 200
folder = "results/J1010_2022-02-04-23:58:58"
time_string = folder.split('/')[-1].split('_')[-1]
data = blazar_utils.read_data("real_data/J1010_SED_reduced.dat")
configs = blazar_utils.read_configs()
minima, maxima = blazar_utils.min_max_parameters(alpha2_limits=configs["alpha2_limits"])
reader = emcee.backends.HDFBackend(folder + "/backend.h5")
flat_log_probs = reader.get_log_prob(flat=True, discard=discard)
log_probs = reader.get_log_prob(discard=discard)
flat_samples = reader.get_chain(flat=True, discard=discard)
samples = reader.get_chain(discard=discard)
indices_within_1sigma = blazar_report.get_indices_within_1sigma(flat_log_probs, 9)
best_chi_sq_index = np.argmax(flat_log_probs)
best_chi_sq = -2 * flat_log_probs[best_chi_sq_index]
best_params = flat_samples[best_chi_sq_index]
min_1sigma, max_1sigma = blazar_report.min_max_params_1sigma(flat_samples, flat_log_probs)

cerruti_params = [96.83, 2.40, 2.0, 4.0, 2, 6.70, 4.725, -1.82, 16.11]


def test_mcmc():
    p0_string = "random"
    p0_data = blazar_utils.random_defaults(configs, minima, maxima)

    now = datetime.now()


def test_extrapolation(max_num_plot=200):
    # TODO before using, put back the modified chi sq in chi_squ func
    # load params from saved
    with open(BASE_PATH + "other/extrapolation_test_params.txt", 'r') as f:
        params = []
        current = ""
        for line in f:
            if line == "\n":
                continue
            if line.count('[') != 0:
                if len(current) != 0:
                    new_param = [float(i) for i in
                                 current.replace(']', '').replace('[', '').replace('\n', ' ').strip().split()]
                    params.append(new_param)
                current = ""
            current = current + " " + line

    num_plots = min(max_num_plot, len(params))

    axes = AxesSequence(num_plots=num_plots)
    if num_plots != len(params):
        indices = np.random.choice(np.arange(0, len(params)), size=num_plots, replace=False)
    else:
        indices = np.arange(0, len(params))

    for i, ax in zip(indices, axes):
        if not np.isfinite(blazar_utils.log_prior(params[i], minima, maxima)):
            continue
        model = blazar_model.make_model(params[i], "test_run", parameter_file="parameter_files/test_param.txt")
        logv_all = model[0]
        logvFv_all = model[1]

        func = interpolate.interp1d(logv_all, logvFv_all, fill_value='extrapolate')

        greater = np.where(data[0] > np.max(model[2]))[0]
        x = np.log10(data[0])
        y = func(x)

        plt.plot(x, y, 'o')
        x = np.log10(data[0][:greater[0]])
        plt.plot(x, func(x), '-o')

        chi_sq = blazar_utils.chi_squared_from_model(model, data[0], data[1], data[2])

        plt.title("SED from Model (blue points are extrapolated)\n" + r"$\chi^2$: " + "{:.3e}".format(chi_sq))
        plt.xlabel(r"$\nu$")
        plt.ylabel(r"$\nu F_{\nu}$")
    return axes, num_plots


class AxesSequence(object):
    """
    From https://stackoverflow.com/questions/13443474/matplotlib-sequence-of-figures-in-the-same-window
    By Joe Kington, Nov 18, 2012

    Creates a series of axes in a figure where only one is displayed at any
    given time. Which plot is displayed is controlled by the arrow keys.
    """

    def __init__(self, num_plots=50, window_name="Models with Extrapolated Points"):
        self.fig = plt.figure(window_name)
        self.axes = []
        self._i = 0  # Currently displayed axes index
        self._n = 0  # Last created axes index
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)
        self.num_plots = num_plots

    def __iter__(self):
        while True:
            yield self.new()

    def new(self):
        # The label needs to be specified so that a new axes will be created
        # instead of "add_axes" just returning the original one.
        ax = self.fig.add_axes([0.15, 0.1, 0.8, 0.8],
                               visible=False, label=self._n)
        self._n += 1
        self.axes.append(ax)
        return ax

    def on_keypress(self, event):
        if event.key == 'right':
            self.next_plot()
        elif event.key == 'left':
            self.prev_plot()
        else:
            return
        self.fig.canvas.draw()

    def next_plot(self):
        if self._i + 1 < len(self.axes):  # i added +1, think was error
            self.axes[self._i].set_visible(False)
            self.axes[self._i + 1].set_visible(True)
            self._i += 1
            print(" Plot number (out of " + str(self.num_plots) + "): " + str(self._i + 1) + "   ", end='\r')

    def prev_plot(self):
        if self._i > 0:
            self.axes[self._i].set_visible(False)
            self.axes[self._i - 1].set_visible(True)
            self._i -= 1
            print(" Plot number (out of " + str(self.num_plots) + "): " + str(self._i + 1) + "   ", end='\r')

    def show(self):
        self.axes[0].set_visible(True)
        print(" Plot number (out of " + str(self.num_plots) + "): " + str(self._i + 1) + "   ", end='\r')
        plt.show()


def load_plot():
    f = open(BASE_PATH + "results/J1010_2022-02-04-23:58:58/plot_with__params.pickle", 'rb')
    fig = pickle.load(f)
    f.close()
    plt.title("J1010 MCMC Results")
    c_model = blazar_model.make_model(cerruti_params, "test_run", parameter_file="parameter_files/test_param.txt")
    print(blazar_utils.log_prob_from_model(c_model, data[0], data[1], data[2]) * -2)
    print(blazar_utils.log_prob_from_model(c_model, data[0], data[1], data[2]) * -2 / 24)
    plt.plot(c_model[0], c_model[1], 'r--', label="Cerruti model", alpha=.75, linewidth=1)
    plt.legend(loc='best')
    f.close()
    return fig


def run_plot(extreme=True):
    blazar_plots.plot_1sigma(data[0], data[1], data[2], indices_within_1sigma, flat_samples, best_chi_sq_index,
                             both=False, extreme=False, folder="testing", save=False, show=False, serialize=False,
                             lower_adjust_multiplier=1.3, max_num_lines_to_graph=200)


def run_corner():
    blazar_plots.corner_plot(flat_samples, minima, maxima, best_params, min_1sigma, max_1sigma,
                             title="J1010 Corner Plot (1000 num, 150 walkers, discard=200)",
                             file_name=folder + "/new_augmented_corner.svg", save=True, show=False)


def load_extreme_plot():
    with open(BASE_PATH + folder + "/plot_with_extreme_params.pickle", 'rb') as f:
        fig = pickle.load(f)
    plt.savefig(BASE_PATH + folder + "/plot_with_extreme_params." + blazar_plots.image_type)
    # plt.show()


def get_tau_var_value(params, tau_var=24):
    # this value should be <= 1
    tau_var = tau_var * 60 * 60
    c = 2.997924 * 1.0e+10
    z = 0.143
    delta = params[0]
    R = np.power(10, params[-1])
    return tau_var >= R * (1 + z) / (c * delta)


def run_save_plots():
    blazar_report.save_plots_and_info(configs, data, minima, maxima, folder=folder, reader=reader, p0_source="random",
                                      acceptance_frac=0.1589575)


def plot_log_plot_model():
    model = blazar_model.make_model(best_params, "test_run", parameter_file="parameter_files/test_param.txt")
    v = model[2]
    vFv = model[3]
    fig, ax = plt.subplots()
    ax.set_xlabel(r"Frequency $\nu$ (Hz)")
    ax.set_ylabel(r"Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")
    # ax.set_title("RXSJ101015.9-311900 Data")
    plt.yscale("log")
    plt.xscale("log")
    ax.plot(model[2], model[3], '-')
    # ax.plot(v, vFv, '.')
    # plt.figtext(0.99, 0.01, 'Data from Cerruti (2013)', horizontalalignment='right', fontsize=8)
    ax.set_xlim(v[0], v[-1])
    plt.errorbar(data[0], data[1], yerr=(data[2], data[2]), fmt='.', ecolor='k')
    new_min, new_max = blazar_plots.scale_to_values(data[1], lower_adjust_multiplier=15,
                                                    upper_adjust_multiplier=3)
    ax.set_ylim(new_min, new_max)

    h = 4.135667662E-15

    def v_to_e(val):
        return val * h

    def e_to_v(val):
        return val / h

    secax = ax.secondary_xaxis('top', functions=(v_to_e, e_to_v))
    secax.set_xlabel("Energy (eV)")

    plt.vlines(3 * 10 ** 11, new_min, new_max, linestyles='dashed', colors='k', linewidth=.5)
    l, r = plt.xlim()
    plt.tight_layout()
    plt.subplots_adjust(top=0.84)

    # plt.savefig(BASE_PATH + "other/presentation_plots/data_line.png", dpi=400)


if __name__ == '__main__':
    """
    axis_sequence, models_plotted = test_extrapolation()
    print("Plot number (out of " + str(models_plotted) + "):")
    axis_sequence.show()
    """
    # blazar_plots.plot_chi_squared(log_probs, discard, plot_type='best', save=True, show=False)
    # plt.show()
    # output = blazar_mcmc.text_report(best_params, min_1sigma, max_1sigma, best_chi_sq)
    # print(output)
    plot_log_plot_model()
    plt.show()
