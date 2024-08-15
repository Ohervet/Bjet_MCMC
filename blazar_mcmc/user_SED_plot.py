#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: user_SED_plot.py

Script to produce custom SED plot. The user can change anything in this script without impacting 
the standard behavior of Bjet_MCMC. 

All standard outputs of Bjet_MCMC can be re-created with the option Full_outputs = True

"""
import os
import emcee
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

import blazar_model
import blazar_report
import blazar_utils
import blazar_plots
from blazar_mcmc.blazar_properties import *

__all__ = ["e_to_v", "v_to_e"]

cmap = plt.get_cmap("tab10")
filled_markers = ("o", "^", "<", ">", "8", "s", "p", "*", "h", "H", "D", "d", "P", "X")
marker_size = 4
params = {  #'backend': 'ps',
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
}
plt.rcParams.update(params)
plt.close()


# -----FLAGS-----#

# recreate all text files and plots outputs
Full_outputs = True

# user defined SED plot boundaries
boundaries = "default"
# if you want to set boundaries (nu in Hz, nuFnu in erg.cm-2.s-1):
# boundaries = [nu_min,nu_max, nuFnu_min, nuFnu_max]

# type of residual plot to be displayed
residual = "sigma"


backend_file = input("Enter relative path to backend file: ")
# example: local_results/J1010/J1010_2023-07-04-23:03:45/backend.h5

# local_results/J1010_2024-08-12-16:47:25/backend.h5


# frequency/energy conversion
h = 4.135667662e-15


def v_to_e(val):
    return val * h


def e_to_v(val):
    return val / h


folder = backend_file[: backend_file.rfind("/")]

if os.path.exists(FOLDER_PATH + folder + "/info.txt"):
    info = blazar_report.parse_info_doc(folder + "/info.txt")
    configs = info["configs"]
    print(configs)
else:
    configs = blazar_utils.read_configs()
data = blazar_utils.read_data(configs["data_file"], instrument=True)


v_data, vFv_data, err_data, instrument_data, nubin_data = data
# setting limits
uplims = [False] * len(v_data)
lolims = [False] * len(v_data)
for i in range(len(err_data[1])):
    if err_data[0][i] == 0:
        uplims[i] = True
        err_data[0][i] = vFv_data[i] / 4
    if err_data[0][i] == -1:
        lolims[i] = True
        err_data[1][i] = vFv_data[i] / 4

redshift = configs["redshift"]
discard = configs["discard"]
eic = configs["eic"]
fixed_params = configs["fixed_params"]


param_min_vals, param_max_vals = blazar_utils.min_max_parameters(
    eic=eic, fixed_params=fixed_params
)
reader = emcee.backends.HDFBackend(FOLDER_PATH + backend_file)
chain = reader.get_chain(discard=discard)
# all parameters values:
flat_samples = reader.get_chain(flat=True, discard=discard)

log_probs = reader.get_log_prob(discard=discard)
flat_log_probs = reader.get_log_prob(flat=True, discard=discard)

best_log_prob, best_params = blazar_report.get_best_log_prob_and_params(
    log_probs=flat_log_probs, flat_chain=flat_samples
)
indices_within_1sigma = blazar_report.get_indices_within_1sigma(
    flat_log_probs, eic=eic, ndim=modelProperties(eic).NUM_DIM
)

name_stem = "plot"
command_params_1, command_params_2 = blazar_model.command_line_sub_strings(
    name_stem=name_stem, redshift=redshift, prev_files=False, eic=eic
)
command_params_2[3] = "99"  # number of points used to make SED

best_model = blazar_model.make_model(
    best_params,
    name_stem=name_stem,
    redshift=redshift,
    command_params_1=command_params_1,
    command_params_2=command_params_2,
    eic=eic,
    fixed_params=fixed_params,
)


min_per_param, max_per_param = blazar_plots.get_params_1sigma_ranges(
    flat_samples, indices_within_1sigma, eic=eic, fixed_params=fixed_params
)
to_plot = np.concatenate((np.array(min_per_param), np.array(max_per_param)))
logv = best_model[0].copy()
lowest_per_point, highest_per_point = blazar_plots.get_min_max_per_point(
    logv,
    to_plot,
    name_stem=name_stem,
    redshift=redshift,
    command_params_1=command_params_1,
    command_params_2=command_params_2,
    eic=eic,
    fixed_params=fixed_params,
)


if residual:
    fig, (ax, ax1) = plt.subplots(
        2, sharex="all", gridspec_kw={"height_ratios": [6, 1.5]}, figsize=(9, 7)
    )
else:
    fig, ax = plt.subplots()
plt.xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")

ax.set_yscale("log")
plt.xscale("log")

ax.fill_between(
    np.power(10, logv),
    np.power(10, lowest_per_point),
    np.power(10, highest_per_point),
    alpha=0.5,
    label=r"Within 1$\sigma$",
)

ax.plot(best_model[2], best_model[3], label="Best model", zorder=2.5, alpha=0.8)

# This line re-run the model to import the good associated .dat files to plot
blazar_model.make_model(
    best_params,
    name_stem=name_stem,
    redshift=redshift,
    command_params_1=command_params_1,
    command_params_2=command_params_2,
    eic=eic,
    fixed_params=fixed_params,
)

synchrotron_model = np.loadtxt(
    FOLDER_PATH + "sed_calculations/" + name_stem + "_ss.dat", delimiter=" "
)
logv_synchrotron = synchrotron_model[:, 0]
logvFv_synchrotron = synchrotron_model[:, 2]
v_synchrotron = np.power(10, logv_synchrotron)
vFv_synchrotron = np.power(10, logvFv_synchrotron)

compton_model = np.loadtxt(
    FOLDER_PATH + "sed_calculations/" + name_stem + "_cs.dat", delimiter=" "
)
logv_compton = compton_model[:, 0]
logvFv_compton = compton_model[:, 2]
v_compton = np.power(10, logv_compton)
vFv_compton = np.power(10, logvFv_compton)

cs2 = np.loadtxt(
    FOLDER_PATH + "sed_calculations/" + name_stem + "_cs2.dat", delimiter=" "
)
logv_cs2 = cs2[:, 0]
logvFv_cs2 = cs2[:, 2]
v_cs2 = np.power(10, logv_cs2)
vFv_cs2 = np.power(10, logvFv_cs2)

if boundaries == "default":
    ax.set_xlim(5e8, 1e28)
    new_min, new_max = blazar_plots.scale_to_values(
        data[1], lower_adjust_multiplier=20, upper_adjust_multiplier=15
    )
    ax.set_ylim(new_min, new_max)
else:
    ax.set_xlim(boundaries[0], boundaries[1])
    ax.set_ylim(boundaries[2], boundaries[3])

ax.plot(v_synchrotron, vFv_synchrotron, "k--", label="Synchrotron", alpha=0.5)
ax.plot(v_compton, vFv_compton, "k-", label="Self Compton", alpha=0.5)
if max(vFv_cs2) > new_min:
    ax.plot(v_cs2, vFv_cs2, "-", label="2nd Order SC")


if eic:
    ecs = np.loadtxt(
        FOLDER_PATH + "sed_calculations/" + name_stem + "_ecs.dat", delimiter=" "
    )
    logv_ecs = ecs[:, 0]
    logvFv_ecs = ecs[:, 2]
    v_ecs = np.power(10, logv_ecs)
    vFv_ecs = np.power(10, logvFv_ecs)

    nuc = np.loadtxt(
        FOLDER_PATH + "sed_calculations/" + name_stem + "_nuc.dat", delimiter=" "
    )
    logv_nuc = nuc[:, 0]
    logvFv_nuc = nuc[:, 2]
    v_nuc = np.power(10, logv_nuc)
    vFv_nuc = np.power(10, logvFv_nuc)

    if max(vFv_ecs) > new_min:
        ax.plot(v_ecs, vFv_ecs, "-", label="EIC")
    ax.plot(v_nuc, vFv_nuc, "-", label="nucleus")


list_intruments = [instrument_data[0]]
v_data_inst = [v_data[0]]
vFv_data_inst = [vFv_data[0]]
err_data_inst_down = [err_data[0][0]]
err_data_inst_up = [err_data[1][0]]
nubin_data_inst_low = [nubin_data[0][0]]
nubin_data_inst_high = [nubin_data[1][0]]
uplims_inst = [uplims[0]]
lolims_inst = [lolims[0]]
tmp = 0
tmp2 = 0

for i in range(1, len(instrument_data)):
    # cycle colors & markers
    if len(list_intruments) - tmp >= cmap.N:
        tmp += cmap.N
    color_index = len(list_intruments) - tmp
    if len(list_intruments) - tmp2 >= len(filled_markers):
        tmp2 += len(filled_markers)
    marker_index = len(list_intruments) - tmp2

    if instrument_data[i] != list_intruments[-1]:
        ax.errorbar(
            v_data_inst,
            vFv_data_inst,
            xerr=(nubin_data_inst_low, nubin_data_inst_high),
            yerr=(err_data_inst_down, err_data_inst_up),
            uplims=uplims_inst,
            lolims=lolims_inst,
            label=str(list_intruments[-1]),
            markersize=marker_size,
            elinewidth=1,
            color=cmap(color_index),
            fmt=filled_markers[marker_index - 1],
        )
        list_intruments.append(instrument_data[i])
        v_data_inst = [v_data[i]]
        vFv_data_inst = [vFv_data[i]]
        err_data_inst_down = [err_data[0][i]]
        err_data_inst_up = [err_data[1][i]]
        nubin_data_inst_low = [nubin_data[0][i]]
        nubin_data_inst_high = [nubin_data[1][i]]
        uplims_inst = [uplims[i]]
        lolims_inst = [lolims[i]]
    else:
        v_data_inst.append(v_data[i])
        vFv_data_inst.append(vFv_data[i])
        err_data_inst_down.append(err_data[0][i])
        err_data_inst_up.append(err_data[1][i])
        nubin_data_inst_low.append(nubin_data[0][i])
        nubin_data_inst_high.append(nubin_data[1][i])
        uplims_inst.append(uplims[i])
        lolims_inst.append(lolims[i])
    if i == len(instrument_data) - 1:
        ax.errorbar(
            v_data_inst,
            vFv_data_inst,
            xerr=(nubin_data_inst_low, nubin_data_inst_high),
            yerr=(err_data_inst_down, err_data_inst_up),
            uplims=uplims_inst,
            lolims=lolims_inst,
            label=str(list_intruments[-1]),
            markersize=marker_size,
            elinewidth=1,
            color=cmap(color_index),
            fmt=filled_markers[marker_index - 1],
        )


ax.legend(loc="upper center", ncol=4)


secax = ax.secondary_xaxis("top", functions=(v_to_e, e_to_v))
secax.set_xlabel("Energy (eV)")

if residual:

    f_best = interpolate.interp1d(
        best_model[2], best_model[3], fill_value="extrapolate"
    )
    f_low = interpolate.interp1d(
        best_model[2], np.power(10, lowest_per_point), fill_value="extrapolate"
    )
    f_high = interpolate.interp1d(
        best_model[2], np.power(10, highest_per_point), fill_value="extrapolate"
    )

    ax1.hlines(0, best_model[2][0], best_model[2][-1])

    if residual == "delta":
        ax1.fill_between(
            best_model[2],
            f_high(best_model[2]) / f_best(best_model[2]) - 1,
            -f_best(best_model[2]) / f_low(best_model[2]) + 1,
            alpha=0.5,
            zorder=1,
        )
    if residual == "sigma":
        for i in range(len(err_data[1])):
            if err_data[1][i] == 0:
                err_data[0][i] = 0

        err_low = err_data[0]
        err_high = err_data[1]

        diff = vFv_data - f_best(v_data)
        # check if data is above or below model flux points, True if above
        sign = diff > 0
        # transform err_data to consider ULs (P=95% when below, P=5% when above)
        if_UL = err_data[0] == 0
        err_data[0] += if_UL * blazar_utils.chi_squared_Limit_to_err(
            0.95, f_best(v_data), vFv_data
        )
        err_data[1] += if_UL * blazar_utils.chi_squared_Limit_to_err(
            0.05, f_best(v_data), vFv_data
        )
        # transform err_data to consider LLs (P=95% when above, P=5% when below)
        if_LL = err_data[0] == -1
        err_data[0] += if_LL * blazar_utils.chi_squared_Limit_to_err(
            0.05, f_best(v_data), vFv_data
        )
        err_data[1] += if_LL * blazar_utils.chi_squared_Limit_to_err(
            0.95, f_best(v_data), vFv_data
        )

        vFv_sigma = (
            sign * (vFv_data - f_best(v_data)) / err_low
            - ~sign * (f_best(v_data) - vFv_data) / err_high
        )
        vFv_sigma_inst = [vFv_sigma[0]]

    list_intruments = [instrument_data[0]]
    v_data_inst = [v_data[0]]
    vFv_data_inst = [vFv_data[0]]
    err_data_inst_down = [err_low[0]]
    err_data_inst_up = [err_high[0]]
    uplims_inst = [uplims[0]]
    lolims_inst = [lolims[0]]
    tmp = 0
    tmp2 = 0

    for i in range(1, len(instrument_data)):
        # cycle colors & markers
        if len(list_intruments) - tmp >= cmap.N:
            tmp += cmap.N
        color_index = len(list_intruments) - tmp
        if len(list_intruments) - tmp2 >= len(filled_markers):
            tmp2 += len(filled_markers)
        marker_index = len(list_intruments) - tmp2

        if instrument_data[i] != list_intruments[-1]:
            if residual == "delta":
                ax1.errorbar(
                    v_data_inst,
                    vFv_data_inst / f_best(v_data_inst) - 1,
                    yerr=(
                        err_data_inst_down / f_best(v_data_inst),
                        err_data_inst_up / f_best(v_data_inst),
                    ),
                    uplims=uplims_inst,
                    lolims=lolims_inst,
                    markersize=4,
                    elinewidth=1,
                    color=cmap(color_index),
                    fmt=filled_markers[marker_index - 1],
                )
            if residual == "sigma":
                ax1.errorbar(
                    v_data_inst,
                    vFv_sigma_inst,
                    yerr=(1),
                    uplims=uplims_inst,
                    lolims=lolims_inst,
                    markersize=4,
                    elinewidth=1,
                    color=cmap(color_index),
                    fmt=filled_markers[marker_index - 1],
                )
                vFv_sigma_inst = [vFv_sigma[i]]

            list_intruments.append(instrument_data[i])
            v_data_inst = [v_data[i]]
            vFv_data_inst = [vFv_data[i]]
            err_data_inst_down = [err_low[i]]
            err_data_inst_up = [err_high[i]]
            uplims_inst = [uplims[i]]
            lolims_inst = [lolims[i]]
        else:
            if residual == "sigma":
                vFv_sigma_inst.append(vFv_sigma[i])
                v_data_inst.append(v_data[i])
                vFv_data_inst.append(vFv_data[i])
                err_data_inst_down.append(err_low[i])
                err_data_inst_up.append(err_high[i])
                uplims_inst.append(uplims[i])
                lolims_inst.append(lolims[i])

        if i == len(instrument_data) - 1:
            if residual == "delta":
                ax1.errorbar(
                    v_data_inst,
                    vFv_data_inst / f_best(v_data_inst) - 1,
                    yerr=(
                        err_data_inst_down / f_best(v_data_inst),
                        err_data_inst_up / f_best(v_data_inst),
                    ),
                    uplims=uplims_inst,
                    lolims=lolims_inst,
                    markersize=4,
                    elinewidth=1,
                    color=cmap(color_index),
                    fmt=filled_markers[marker_index - 1],
                )
            if residual == "sigma":
                ax1.errorbar(
                    v_data_inst,
                    vFv_sigma_inst,
                    yerr=(1),
                    uplims=uplims_inst,
                    lolims=lolims_inst,
                    markersize=4,
                    elinewidth=1,
                    color=cmap(color_index),
                    fmt=filled_markers[marker_index - 1],
                )
    if residual == "delta":
        ax1.set_ylim(
            min(vFv_data / f_best(v_data) - 1) * 2,
            max(vFv_data / f_best(v_data) - 1) * 2,
        )
        ax1.set_ylabel("Residuals")
    if residual == "sigma":
        ax1.set_ylabel(r"$\Delta \sigma$")
        ax1.set_ylim(min(vFv_sigma) * 1.5, max(vFv_sigma) * 1.5)
    ax1.yaxis.get_ticklocs(minor=True)
    ax1.xaxis.get_ticklocs(minor=True)
    ax1.minorticks_on()

plt.tight_layout()
# plt.show()
plt.savefig(BASE_PATH + folder + "/user_plot_SED.svg")


if Full_outputs:
    blazar_report.save_plots_and_info(
        configs,
        (v_data, vFv_data, err_data),
        param_min_vals,
        param_max_vals,
        folder=folder,
        samples=(chain, flat_samples, log_probs, flat_log_probs),
        use_samples=True,
        redshift=redshift,
        eic=eic,
        verbose=True,
    )
