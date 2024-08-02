"""

"""

import glob
import os
import random
import sys

import emcee
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy import interpolate
import corner
from astropy.io import ascii

import blazar_model
import blazar_report
import blazar_utils
import blazar_run_mcmc
import blazar_plots
import blazar_initialize
import blazar_clean
from blazar_properties import *
cmap = plt.get_cmap("tab10")
filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
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

Fermi_spectrum = False
VERITAS_spectrum = False

Source = "NGC1275_Jan2"

#-----input files and plotting parameters-----#
if Source == "NGC1275_Jan1":
    eic = True
    jet = True
    Ylim = [1e-14,1e-8] #erg cm-2 s-1
    Xlim = [2e8,1e28] #Hz
    parameter_file = "parameter_files/NGC_1275/Jan1/NGC_1275_2017_Jan1.par"
    #parameter_file = "parameter_files/test_params.par"
    instrument_data_file = "real_data/NGC_1275/NGC1275_2017Jan1_SED.dat"
    
if Source == "NGC1275_Jan2":
    eic = True
    jet = True
    Ylim = [1e-14,1e-8] #erg cm-2 s-1
    Xlim = [2e8,1e28] #Hz
    parameter_file = "parameter_files/NGC_1275/Jan2/NGC_1275_2017_Jan2.par"
    instrument_data_file = "real_data/NGC_1275/NGC1275_2017Jan2_SED.dat"


if Source == "PKS1222":
    Ylim = [1e-15,1e-9] #erg cm-2 s-1
    Xlim = [1e8,1e28] #Hz
    #parameter_file = "parameter_files/PKS_1222/PKS_1222_EIC_v4.par"
    parameter_file = "parameter_files/PKS_1222/PKS_1222_BjetMCMC_test.par"
    instrument_data_file = "real_data/PKS_1222+216/PKS1222+216_SED.dat"

    
if Source == "J1010":
    Ylim = [1e-15,1e-9] #erg cm-2 s-1
    Xlim = [1e8,1e28] #Hz
    parameter_file = "parameter_files/OJ_287/Low/OJ287_LowState.par"
    instrument_data_file = "real_data/J1010_SED_reduced.dat"
    
if Source == "OJ287_low":
    eic = True
    jet = True
    Ylim = [1e-15,1e-9] #erg cm-2 s-1
    Xlim = [1e8,1e28] #Hz
    #Low state
    parameter_file = "parameter_files/OJ_287/Low/OJ287_LowState.par"
    instrument_data_file = "real_data/OJ_287/Low/OJ287_Low_SED.dat"
    Fermi_spectrum = "real_data/OJ_287/Low/Fermi_spectrum_Low.dat"
    core_sync = "parameter_files/OJ_287/Low/F_core_syn_LowState_v7.dat"
    core_SSC = "parameter_files/OJ_287/Low/F_core_com_LowState_v7.dat"
if Source == "OJ287_flare":
    eic = True
    jet = True
    Ylim = [1e-15,1e-9] #erg cm-2 s-1
    Xlim = [1e8,1e28] #Hz
    parameter_file = "parameter_files/OJ_287/Flare/OJ287_Flare.par"
    instrument_data_file = "real_data/OJ_287/Flare/OJ287_Flare_SED.dat"
    Fermi_spectrum = "real_data/OJ_287/Flare/Fermi_spectrum_Flare.dat"
    core_sync = "parameter_files/OJ_287/Low/F_core_syn_LowState_v7.dat"
    core_SSC = "parameter_files/OJ_287/Low/F_core_com_LowState_v7.dat"
if Source == "OJ287_Postflare":
    eic = True
    jet = True
    Ylim = [1e-15,1e-9] #erg cm-2 s-1
    Xlim = [1e8,1e28] #Hz
    parameter_file = "parameter_files/OJ_287/PostFlare/OJ287_PostFlare.par"
    instrument_data_file = "real_data/OJ_287/PostFlare/OJ287_PostFlare_SED.dat"
    Fermi_spectrum = "real_data/OJ_287/PostFlare/Fermi_spectrum_PostFlare.dat"
    core_sync = "parameter_files/OJ_287/PostFlare/F_core_syn_PostFlare_v7.dat"
    core_SSC = "parameter_files/OJ_287/PostFlare/F_core_com_PostFlare_v7.dat"
    #VERITAS_spectrum = "real_data/OJ_287/PostFlare/VERITAS_spectrum_PostFlare.dat"





#-----run Bjet-----#
model = blazar_model.file_make_SED(parameter_file=parameter_file, data_folder=None, executable=None, prev_files=False, 
                                    verbose=True)


#-----read instrumental SED-----#
v_data, vFv_data, err_data, instrument_data, nubin_data = blazar_utils.read_data(instrument_data_file, instrument=True)
#setting upper limits
uplims = [False]*len(v_data)
for i in range(len(err_data[1])):
    if err_data[0][i] == 0:
        uplims[i] = True
        err_data[0][i] = vFv_data[i]/4
        
if Source[:5] == "OJ287":
    for i in range(len(instrument_data)):
        if instrument_data[i] in ('KUEHR','Planck_PCCS','VLSS'):
            instrument_data[i] = 'Archives'



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
    
    # EIC blob-jet
    eic_jet = np.loadtxt(FOLDER_PATH + data_folder + "/" + name_stem + "_ecs_jet.dat", delimiter=' ')
    logv_eic_jet = eic_jet[:, 0]
    logvFv_eic_jet = eic_jet[:, 2]
    v_eic_jet = np.power(10, logv_eic_jet)
    vFv_eic_jet = np.power(10, logvFv_eic_jet)
    logv_sum, logvFv_sum, v_sum, vFv_sum = blazar_model.add_data((logv_sum, logvFv_sum, v_sum, vFv_sum),
                                    file_suffix='ecs_jet', name_stem=name_stem, data_folder=data_folder)

if Source[:5] == "OJ287":
    #radio core
    #sync
    core_sync = np.loadtxt(FOLDER_PATH + core_sync, delimiter=' ')
    v_core_sync = core_sync[:, 0]
    vFv_core_sync = core_sync[:, 2]
    #SSC
    core_SSC = np.loadtxt(FOLDER_PATH + core_SSC, delimiter=' ')
    v_core_SSC = core_SSC[:, 0]
    vFv_core_SSC = core_SSC[:, 2]



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

p1 = []
marker1 = []
labels1 = []
markersize = 4
n=0
for i in range(1,len(instrument_data)):
    if instrument_data[i] != list_intruments[-1]:
        marker1.append(filled_markers[len(list_intruments)-1])
        labels1.append(str(list_intruments[-1]))
        if Source[:5] == "OJ287" and len(list_intruments) == 1:
            n=1
            p1.append(ax.errorbar(v_data_inst, vFv_data_inst, yerr=(err_data_inst_down, err_data_inst_up), label=labels1[-1], 
                        markersize=markersize, elinewidth=1, color = '0.7', fmt=marker1[-1],zorder=0))
        else:
            p1.append(ax.errorbar(v_data_inst, vFv_data_inst, xerr=(nubin_data_inst_low, nubin_data_inst_high), yerr=(err_data_inst_down, err_data_inst_up), 
                                  uplims = uplims_inst, label=labels1[-1],  markersize=markersize, elinewidth=1, color = cmap(len(list_intruments)-n), fmt=marker1[-1],zorder=0))
        p1[-1] = p1[-1][0]
        list_intruments.append(instrument_data[i])
        v_data_inst = [v_data[i]]
        vFv_data_inst = [vFv_data[i]]
        err_data_inst_down = [err_data[0][i]]
        err_data_inst_up = [err_data[1][i]]
        nubin_data_inst_low= [nubin_data[0][i]]
        nubin_data_inst_high= [nubin_data[1][i]]
        uplims_inst = [uplims[i]]
    else:   
        v_data_inst.append(v_data[i])
        vFv_data_inst.append(vFv_data[i])
        err_data_inst_down.append(err_data[0][i])
        err_data_inst_up.append(err_data[1][i])
        nubin_data_inst_low.append(nubin_data[0][i])
        nubin_data_inst_high.append(nubin_data[1][i])
        uplims_inst.append(uplims[i])
    if i == len(instrument_data)-1:
        marker1.append(filled_markers[len(list_intruments)-1])
        labels1.append(str(list_intruments[-1]))
        p1.append(ax.errorbar(v_data_inst, vFv_data_inst, xerr=(nubin_data_inst_low, nubin_data_inst_high), yerr=(err_data_inst_down, err_data_inst_up), 
                              uplims = uplims_inst, label=str(list_intruments[-1]), markersize=markersize, elinewidth=1, color = cmap(len(list_intruments)), fmt=marker1[-1],zorder=0))
        p1[-1] = p1[-1][0]        
        
legend1 = ax.legend(p1, loc="upper center", labels=labels1, ncol=4,fontsize=9)
ax.add_artist(legend1)

if Fermi_spectrum:
    Spectrum_table =  ascii.read(FOLDER_PATH + Fermi_spectrum, format='csv',data_start = 2,delimiter='\s')
    E = Spectrum_table["E"] #[eV]
    v = e_to_v(E) #[Hz]
    F = Spectrum_table["E**2.dN/dE"] #[erg cm-2 s-1]
    F_err_low = Spectrum_table["E**2.dN/dE_ErrLow"] #[erg cm-2 s-1]
    F_err_high = Spectrum_table["E**2.dN/dE_ErrHigh"] #[erg cm-2 s-1]
    ax.fill_between(v,F_err_high,F_err_low,color = cmap(list_intruments.index('Fermi-LAT')), alpha=0.4,zorder=0)
    
if VERITAS_spectrum:
    Spectrum_table =  ascii.read(FOLDER_PATH + VERITAS_spectrum, format='csv',data_start = 2,delimiter='\s')
    v = Spectrum_table["nu"] #[Hz]
    #v = e_to_v(E) #[Hz]
    F = Spectrum_table["E**2.dN/dE"] #[erg cm-2 s-1]
    F_err_low = Spectrum_table["E**2.dN/dE_ErrLow"] #[erg cm-2 s-1]
    F_err_high = Spectrum_table["E**2.dN/dE_ErrHigh"] #[erg cm-2 s-1]
    ax.fill_between(v,F_err_high,F_err_low,color = cmap(list_intruments.index('VERITAS')), alpha=0.4,zorder=0)
        
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
    labels2.append("Sum")
    p2.append(ax.plot(v_sum, vFv_sum, color='brown', alpha=.7))
    p2[-1] = p2[-1][0]
    
if Source[:5] == "OJ287":
    labels2.append("8.6 pc core Sync. & SSC")
    p2.append(ax.plot(v_core_sync, vFv_core_sync, color='r', alpha=.7, ls=(0, (1, 1))))
    p2[-1] = p2[-1][0]
    ax.plot(v_core_SSC, vFv_core_SSC, color='r', alpha=.7, ls=(0, (1, 1)))

lines, labels = ax.get_legend_handles_labels()
legend2 = ax.legend(handles=p2, loc='lower center', labels=labels2, ncol=3,fontsize=9)

ax.loglog()
ax.set_ylim(Ylim[0],Ylim[1])
ax.set_xlim(Xlim[0],Xlim[1])
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")
secax = ax.secondary_yaxis('right')
secax = ax.secondary_xaxis('top', functions=(v_to_e, e_to_v))
secax.set_xlabel("Energy (eV)")


#Calculate gamma-ray Chi2 for OJ 287
model_result = [np.log10(v_sum), np.log10(vFv_sum)]
#remove ULs
v_data_noUL = []
vFv_data_noUL = []
err_data_low_noUL = []
err_data_up_noUL = []
instrument_data_noUL = []
for i in range(len(v_data)):
    if err_data[1][i] != 0:
        v_data_noUL.append(v_data[i])
        vFv_data_noUL.append(vFv_data[i])
        err_data_low_noUL.append(err_data[0][i])
        err_data_up_noUL.append(err_data[1][i])    
        instrument_data_noUL.append(instrument_data[i])  
        
        
n = list(instrument_data_noUL).index('Fermi-LAT')
err_cut = [err_data_low_noUL[n:],err_data_up_noUL[n:]]
Chi2 = blazar_utils.chi_squared_from_model(model_result, v_data_noUL[n:], vFv_data_noUL[n:], err_cut)
dof = len(v_data_noUL[n:])-1

print('reduce Chi2 gamma=',Chi2/dof)

#distaces blob to be checked (0.5 to 20 pc in cm):
# array([1.54300000e+18, 2.00813893e+18, 2.61349447e+18, 3.40133506e+18,
#        4.42667101e+18, 5.76109552e+18, 7.49778368e+18, 9.75799826e+18,
#        1.26995568e+19, 1.65278513e+19, 2.15101890e+19, 2.79944575e+19,
#        3.64334155e+19, 4.74163059e+19, 6.17100000e+19])


