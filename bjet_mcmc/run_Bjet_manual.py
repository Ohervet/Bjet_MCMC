"""

"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

import blazar_model
import blazar_utils
from blazar_properties import *
cmap = plt.get_cmap("tab10")
filled_markers = ("o", "^", "<", ">", "8", "s", "p", "*", "h", "H", "D", "d", "P", "X")
marker_size = 4
linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),
     ('long dash with offset', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]


h = 4.135667662E-15

def v_to_e(val):
    return val * h

def e_to_v(val):
    return val / h


Source = "B3_2247"

#-----input files and plotting parameters-----#
if Source == "J1010":
    Ylim = [1e-15,1e-9] #erg cm-2 s-1
    Xlim = [1e8,1e28] #Hz
    eic = False
    jet = False
    parameter_file = "parameter_files/J1010.par"
    instrument_data_file = "real_data/J1010_SED_reduced.dat"
    

if Source == "B3_2247":
    Ylim = [1e-15,1e-9] #erg cm-2 s-1
    Xlim = [1e8,1e28] #Hz
    eic = False
    jet = False
    parameter_file = "parameter_files/B3_2247+381_updated.par"
    instrument_data_file = "real_data/B32247_sed_data.dat"




#-----run Bjet-----#
model = blazar_model.file_make_SED(parameter_file=parameter_file, data_folder=None, executable=None, prev_files=False, 
                                    verbose=True)


#-----read instrumental SED-----#
v_data, vFv_data, err_data, instrument_data, nubin_data = blazar_utils.read_data(instrument_data_file, instrument=True)
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


#-----read model SED-----#
name_stem = 'test_bj'
data_folder = DATA_FOLDER

#blob synchrotron
blob_sync = np.loadtxt(BASE_PATH + data_folder + "/" + name_stem + "_ss.dat", delimiter=' ')
logv_blob_sync = blob_sync[:, 0]
logvFv_blob_sync= blob_sync[:, 2]
v_blob_sync = np.power(10, logv_blob_sync)
vFv_blob_sync = np.power(10, logvFv_blob_sync)

#blob SSC
blob_ssc = np.loadtxt(BASE_PATH + data_folder + "/" + name_stem + "_cs.dat", delimiter=' ')
logv_blob_ssc = blob_ssc[:, 0]
logvFv_blob_ssc= blob_ssc[:, 2]
v_blob_ssc = np.power(10, logv_blob_ssc)
vFv_blob_ssc = np.power(10, logvFv_blob_ssc)
logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_blob_sync, logvFv_blob_sync, v_blob_sync, vFv_blob_sync),
                                file_suffix='cs', name_stem=name_stem, data_folder=data_folder)

#blob SSC 2nd order
blob_ssc2 = np.loadtxt(FOLDER_PATH + data_folder + "/" + name_stem + "_cs2.dat", delimiter=' ')
logv_blob_ssc2 = blob_ssc2[:, 0]
logvFv_blob_ssc2 = blob_ssc2[:, 2]
v_blob_ssc2 = np.power(10, logv_blob_ssc2)
vFv_blob_ssc2 = np.power(10, logvFv_blob_ssc2)
logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                file_suffix='cs2', name_stem=name_stem, data_folder=data_folder)

if eic:
    # accretion disk
    disk = np.loadtxt(FOLDER_PATH + data_folder + "/" + name_stem + "_nuc.dat", delimiter=' ')
    logv_disk = disk[:, 0]
    logvFv_disk = disk[:, 2]
    v_disk = np.power(10, logv_disk)
    vFv_disk = np.power(10, logvFv_disk)
    logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                    file_suffix='nuc', name_stem=name_stem, data_folder=data_folder)
    
    # EIC blob-BLR
    eic_blr = np.loadtxt(FOLDER_PATH + data_folder + "/" + name_stem + "_ecs.dat", delimiter=' ')
    logv_eic_blr = eic_blr[:, 0]
    logvFv_eic_blr = eic_blr[:, 2]
    v_eic_blr = np.power(10, logv_eic_blr)
    vFv_eic_blr = np.power(10, logvFv_eic_blr)
    logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                    file_suffix='ecs', name_stem=name_stem, data_folder=data_folder)

if jet:
    # #jet synchrotron
    jet_sync = np.loadtxt(FOLDER_PATH + data_folder + "/F_jet_syn.dat", delimiter=' ')
    v_jet_sync = jet_sync[:, 0]
    vFv_jet_sync = jet_sync[:, 2]
    logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                    file_suffix='syn', name_stem="F_jet", data_folder=data_folder)
    
    #jet SSC
    jet_ssc = np.loadtxt(FOLDER_PATH + data_folder + "/F_jet_com.dat", delimiter=' ')
    v_jet_ssc = jet_ssc[:, 0]
    vFv_jet_ssc = jet_ssc[:, 2]
    logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                    file_suffix='com', name_stem="F_jet", data_folder=data_folder)
    
    # EIC blob-jet (emitted from the blob)
    eic_jet = np.loadtxt(FOLDER_PATH + data_folder + "/" + name_stem + "_ecs_jet.dat", delimiter=' ')
    logv_eic_jet = eic_jet[:, 0]
    logvFv_eic_jet = eic_jet[:, 2]
    v_eic_jet = np.power(10, logv_eic_jet)
    vFv_eic_jet = np.power(10, logvFv_eic_jet)
    logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                    file_suffix='ecs_jet', name_stem=name_stem, data_folder=data_folder)
    
    # EIC jet-blob (emitted from the jet)
    jet_eic = np.loadtxt(FOLDER_PATH + data_folder + "/F_jet_eic.dat", delimiter=' ')
    v_jet_eic = jet_eic[:, 0]
    vFv_jet_eic = jet_eic[:, 2]
    logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                    file_suffix='eic', name_stem="F_jet", data_folder=data_folder)



#-----plot instrumental SED-----#
fig, ax = plt.subplots(figsize=(7.25,5))

list_intruments=[instrument_data[0]]
v_data_inst = [v_data[0]]
vFv_data_inst = [vFv_data[0]]
err_data_inst_down = [err_data[0][0]]
err_data_inst_up = [err_data[1][0]]
nubin_data_inst_low= [nubin_data[0][0]]
nubin_data_inst_high= [nubin_data[1][0]]
uplims_inst = [uplims[0]]
lolims_inst = [lolims[0]]
tmp = 0
tmp2 = 0
p1 = []
marker1 = []
labels1 = []
n=0
for i in range(1,len(instrument_data)):
    #cycle colors & markers
    if len(list_intruments)-tmp >= cmap.N:
        tmp += cmap.N
    color_index = len(list_intruments) - tmp
    if len(list_intruments)-tmp2 >= len(filled_markers):
        tmp2 += len(filled_markers)
    marker_index = len(list_intruments) - tmp2

    if instrument_data[i] != list_intruments[-1]:
        ax.errorbar(v_data_inst, vFv_data_inst, xerr=(nubin_data_inst_low, nubin_data_inst_high), yerr=(err_data_inst_down, err_data_inst_up), uplims = uplims_inst,
                    lolims=lolims_inst, label=str(list_intruments[-1]), markersize=marker_size, elinewidth=1, color = cmap(color_index), 
                    fmt=filled_markers[marker_index-1])
        list_intruments.append(instrument_data[i])
        v_data_inst = [v_data[i]]
        vFv_data_inst = [vFv_data[i]]
        err_data_inst_down = [err_data[0][i]]
        err_data_inst_up = [err_data[1][i]]
        nubin_data_inst_low= [nubin_data[0][i]]
        nubin_data_inst_high= [nubin_data[1][i]]
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
    if i == len(instrument_data)-1:
        ax.errorbar(v_data_inst, vFv_data_inst, xerr=(nubin_data_inst_low, nubin_data_inst_high), yerr=(err_data_inst_down, err_data_inst_up), uplims = uplims_inst,
                    lolims=lolims_inst, label=str(list_intruments[-1]), markersize=marker_size, elinewidth=1, color = cmap(color_index), 
                    fmt=filled_markers[marker_index-1])

legend1 = ax.legend(loc="upper center", ncol=4)
ax.add_artist(legend1)

p2 = []
labels2 = []
ax.plot(v_blob_sync, vFv_blob_sync, color=cmap(0), label=None, alpha=.7)
labels2.append("Blob Sync. & SSC")
p2.append(ax.plot(v_blob_ssc, vFv_blob_ssc, color=cmap(0), label="Blob Sync. & SSC", alpha=.7))
p2[-1] = p2[-1][0]
if max(vFv_blob_ssc2)>Ylim[0]:
    labels2.append("Blob SSC 2nd order")
    p2.append(ax.plot(v_blob_ssc2, vFv_blob_ssc2, color=cmap(0), label="Blob SSC 2nd order", alpha=.7, ls=':'))
    p2[-1] = p2[-1][0]
labels2.append("Disk")
if eic:
    p2.append(ax.plot(v_disk, vFv_disk, color=cmap(2), label="disk", alpha=.7, ls='-.'))
    p2[-1] = p2[-1][0]
    if max(vFv_eic_blr)>Ylim[0]:
        labels2.append("EIC blob-BLR")
        p2.append(ax.plot(v_eic_blr, vFv_eic_blr, color=cmap(9), label="EIC blob-BLR", alpha=.7, ls=(0, (3, 1, 1, 1))))
        p2[-1] = p2[-1][0]
if jet:
    ax.plot(v_jet_sync, vFv_jet_sync, color='r', alpha=.7, ls=(0, (3, 1, 1, 1, 1, 1)))
    labels2.append("Jet Sync. & SSC")
    p2.append(ax.plot(v_jet_ssc, vFv_jet_ssc, color='r', label="Jet Sync. & SSC", alpha=.7, ls=(0, (3, 1, 1, 1, 1, 1))))
    p2[-1] = p2[-1][0]
    
    labels2.append("EIC Blob-Jet")
    p2.append(ax.plot(v_eic_jet, vFv_eic_jet, color=cmap(4), label="EIC blob-jet", alpha=.7, ls='--'))
    p2[-1] = p2[-1][0]
    
    labels2.append("EIC Jet-Blob")
    p2.append(ax.plot(v_jet_eic, vFv_jet_eic, color=cmap(5), label="EIC jet-blob", alpha=.7, ls='--'))
    p2[-1] = p2[-1][0]
    labels2.append("Sum")
    p2.append(ax.plot(v_sum, vFv_sum, color='brown', alpha=.7))
    p2[-1] = p2[-1][0]
    
if Source[:5] == "OJ287":
    labels2.append("8.6 pc core Sync. & SSC")
    p2.append(ax.plot(v_core_sync, vFv_core_sync, color='r', alpha=.7, ls=(0, (1, 1))))
    p2[-1] = p2[-1][0]
    ax.plot(v_core_SSC, vFv_core_SSC, color='r', alpha=.7, ls=(0, (1, 1)))

lines, labels = ax.get_legend_handles_labels()
legend2 = ax.legend(handles=p2, loc='lower center', labels=labels2, ncol=3)

ax.loglog()
ax.set_ylim(Ylim[0],Ylim[1])
ax.set_xlim(Xlim[0],Xlim[1])
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")
secax = ax.secondary_yaxis('right')
secax = ax.secondary_xaxis('top', functions=(v_to_e, e_to_v))
secax.set_xlabel("Energy (eV)")
plt.show()

