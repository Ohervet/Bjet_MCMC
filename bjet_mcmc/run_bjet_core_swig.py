#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from bjet_core import bj_core
import blazar_model
import blazar_utils
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import os
import emcee


plot = 1
test_hdf5 = True
general_check = False
eic = 0
#now run it by entering all parameters in inputs
model_type = 0
n_components = 6
redshift = 0.143
theta = 0.57
min_freq=1e+8
max_freq=1e+28

if test_hdf5:
    #fetch list of params and Chi2 from hdf5 file
    results_directory = "/test_debug"
    #swig:
    #backend_path = "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/local_results/J1010_2024-09-12-22:10:28/backend.h5"
    backend_path = "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/local_results/test_debug/backend.h5"
    #Classic:
    #backend_path = "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/local_results/J1010_2023-07-04-23:03:45/backend.h5"
    #backend_path = "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/local_results/test_debug_classic/backend.h5"
    fixed_params= [83.8, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf]
    reader = emcee.backends.HDFBackend(backend_path)
    flat_samples = reader.get_chain(flat=True)
    flat_log_probs = reader.get_log_prob(flat=True)
    
    for j in range(len(flat_samples)-100, len(flat_samples)):
        print(j)
        sample =  np.array([fixed_params[0]]+ list(flat_samples[j]))
        #print(sample)
        
        #check best model
        # sample = np.array([83.79739284,3.6269726,  2.55506677,  3.75140817,  0.96470574 , 6.29819794, 5.32818923, -2.76763863,
        #                    17.14558664])
            
        list_params = [redshift, 69.6, 0.57] + [sample[0]] + [10**sample[1]] + list(sample[2:4]) + list(10**sample[4:]) + [0, 1, 9.0e+17,
                      99, min_freq,max_freq]
        
        chi2_hdf5 = -2*flat_log_probs[j]
        #print("Chi2 (hdf5)=",chi2_hdf5)
        
        #EIC 
        if eic:
            model_type = 1
            n_components = 10
            list_params = [0.069, 69.6, 0.57, 50.298,  2.24e+03,  2.024,  4.166,  5.82e+01, 1.03e+06, 8.46e+03, 2.13e-02, 4.34e+16, 0, 1,
                           1.54e+17, 1.51e+05,1.0e+03, 3.13e+41, 1.40e-03, 1.0e+20, 0.1,99, 1e+8,1e+29]
                       
        
        #./bj_core 3 /home/olivier/Bjet_MCMC_root/Bjet_MCMC/sed_calculations 1 0.143 69.6 0.57 8.571702988940395 11283.043346621951 2.4733878206948723 5.541684562450467 4.145882718669621 55842.50406406743 928.8571578806423 0.0015206064774071673 3.9790008915879885e+18 0 1 6.309289739332383e+18 328500.0765626769 2.0e+4 4.17521574983951e+44 0.007025610341282167 5.5e+20 9.0e-5 99 100000000.0 1e+28 plot
        
        
        #change into strings
        list_params = list((map(str, list_params))) 
        string_params =""
        for i in range(len(list_params)):
            string_params += list_params[i]+" "
            
        #print(string_params)   
        #string_params = "0.143 69.6 0.57 8.3800e+01 4.2624e+03 2.3971e+00 3.0938e+00 1.0951e+01 7.7916e+05 3.2333e+04 4.2967e-03 4.9208e+16 0 1 9.0e+17 99 100000000.0 1e+28 "
        
        NU_DIM = 99
        full_size = int(NU_DIM*n_components)
        #aa = bj_core.main_swig2("/home/olivier/Bjet_MCMC_root/Bjet_MCMC/parameter_files/J1010.par")
        #aa = bj_core.main_swig(*list_params, 3, "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/sed_calculations", "temp_stem")
        aa = bj_core.main_swig(string_params, model_type, "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/sed_calculations", "temp_stem")
        p = (ctypes.c_double * full_size).from_address(int(aa))
        p = list(p)
        #print(p[3*99:4*99])
        
        
        # #this work for 1D array
        start = 0
        stop = NU_DIM
        sliced_output = [0]*n_components
        for i in range(n_components):
            sliced_output[i] = p[start:stop+start]
            start += NU_DIM
            if i%2:
                #remove zeros in frequency array
                while 0 in sliced_output[i-1]:
                    tmp = sliced_output[i-1].index(0)
                    del sliced_output[i][tmp]
                    del sliced_output[i-1][tmp]   
                if i == 1:
                    logv = np.array(sliced_output[0])
                    logvFv = np.array(sliced_output[1])
                    v = 10**logv
                    vFv = 10**logvFv
                    # print(v)
                    # print(vFv,"\n")
                else:
                    logv_new = np.array(sliced_output[i-1])
                    logvFv_new = np.array(sliced_output[i])
                    v_new = 10**logv_new
                    vFv_new = 10**logvFv_new   
                    # if i == 3:
                    #     print("v_new:",v_new)
                    #     print(vFv_new,"\n")
                    logv, logvFv, v, vFv = blazar_model.add_SED_component((logv, logvFv, v, vFv), (logv_new, logvFv_new, v_new, vFv_new)) 
                    # print("Sum")
                    # print(v)
                    # print(vFv)
                    
                    
        #calculate Chi2
        model_results = (logv, logvFv, v, vFv)
        v_data, vFv_data, err_data, instrument_data, nubin_data = blazar_utils.read_data('real_data/J1010_SED_reduced.dat', instrument=True)
        Chi2 = blazar_utils.chi_squared_from_model(model_results, v_data, vFv_data, err_data)
        #print("Chi2_from_model (manual check)=",Chi2)
        
        # #print("nu", v)
        
        
        # #test chi_squared
        # name_stem="test"
        # data_folder="sed_calculations"
        # prev_files=None
        # eic=False
        # command_params_1, command_params_2 = blazar_model.command_line_sub_strings(
        #     name_stem=name_stem,
        #     theta=theta,
        #     redshift=redshift,
        #     min_freq=min_freq,
        #     max_freq=max_freq,
        #     data_folder=data_folder,
        #     prev_files=prev_files,
        #     eic=eic,
        # )
        
        # Chi_squared = blazar_utils.chi_squared(
        #     sample,
        #     v_data,
        #     vFv_data,
        #     err_data,
        #     name_stem=name_stem,
        #     theta=theta,
        #     redshift=redshift,
        #     min_freq=min_freq,
        #     max_freq=max_freq,
        #     data_folder=data_folder,
        #     command_params_1=command_params_1,
        #     command_params_2=command_params_2,
        #     prev_files=prev_files,
        #     eic=eic,
        # )
        # #print("chi_squared (manual check) = ", Chi_squared)
        
        
        # #print("Test log_probability...")
        # Chi2_full_method = -2*blazar_utils.log_probability(
        #     sample,
        #     v_data,
        #     vFv_data,
        #     err_data,
        #     name_stem=None,
        #     theta=theta,
        #     redshift=redshift,
        #     min_freq=min_freq,
        #     max_freq=max_freq,
        #     data_folder=data_folder,
        #     command_params_1=command_params_1,
        #     command_params_2=command_params_2,
        #     prev_files=prev_files,
        #     eic=eic,
        #     use_param_file=False
        # )       
        # #print("Chi2 from log_probability =",Chi2_full_method,"\n")
    
        if np.abs(1-Chi2/chi2_hdf5)>1e-2:
            print("Chi2 (hdf5)=",chi2_hdf5)
            print("Chi2 from model =",Chi2,"\n")
    
     
if plot:
    plt.plot(sliced_output[0],sliced_output[1], label= "Synchrotron")
    plt.plot(sliced_output[2],sliced_output[3], label= "SSC")
    plt.plot(sliced_output[4],sliced_output[5], label= "SSC2")
    if eic:
        plt.plot(sliced_output[6],sliced_output[7], label= "EIC")
        plt.plot(sliced_output[8],sliced_output[9], label= "Disk")
    plt.plot(logv, logvFv, label= "All")
    plt.ylim([-16,-9])

    plt.errorbar(np.log10(v_data),np.log10(vFv_data), fmt ="o")
    
    plt.legend()
    plt.show()




if general_check:
        
    params = np.array([8.3800e+01, 3.629654203154188, 2.3971e+00, 3.0938e+00, 1.0394537789617364, 5.891626648920381, 4.509646002260803,
              -2.366864968143627, 16.69203571401569])
    
    list_params = [redshift, 69.6, 0.57] + [params[0]] + [10**params[1]] + list(params[2:4]) + list(10**params[4:]) + [0, 1, 9.0e+17,
                  99, min_freq, max_freq]
    
    
    list_params = list((map(str, list_params))) 
    string_params =""
    for i in range(len(list_params)):
        string_params += list_params[i]+" "
        
    #print(string_params)   
    #string_params = "0.143 69.6 0.57 8.3800e+01 4.2624e+03 2.3971e+00 3.0938e+00 1.0951e+01 7.7916e+05 3.2333e+04 4.2967e-03 4.9208e+16 0 1 9.0e+17 99 100000000.0 1e+28 "
    
    NU_DIM = 99
    full_size = int(NU_DIM*n_components)
    #aa = bj_core.main_swig2("/home/olivier/Bjet_MCMC_root/Bjet_MCMC/parameter_files/J1010.par")
    #aa = bj_core.main_swig(*list_params, 3, "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/sed_calculations", "temp_stem")
    aa = bj_core.main_swig(string_params, model_type, "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/sed_calculations", "temp_stem")
    p = (ctypes.c_double * full_size).from_address(int(aa))
    p = list(p)
    
    v_data, vFv_data, err_data, instrument_data, nubin_data = blazar_utils.read_data('real_data/J1010_SED_reduced.dat', instrument=True)
    
    
    
    print("test process_model...")
    
    logv, logvFv, v, vFv = blazar_model.process_model(p)
    
    plt.plot(logv, logvFv, ls="--", color='0', label= "Process_model_test")
    
    
    print("test make_model...")
    
    name_stem="test"
    data_folder="sed_calculations"
    prev_files=None
    eic=False
    
    command_params_1, command_params_2 = blazar_model.command_line_sub_strings(
        name_stem=name_stem,
        theta=theta,
        redshift=redshift,
        min_freq=min_freq,
        max_freq=max_freq,
        data_folder=data_folder,
        prev_files=prev_files,
        eic=eic,
    )
    print("command_params_1",command_params_1)
    print("command_params_2",command_params_2)
    
    
    logv, logvFv, v, vFv =  blazar_model.make_model(
              params,
              name_stem=name_stem,
              theta=theta,
              redshift=redshift,
              min_freq=min_freq,
              max_freq=max_freq,
              torus_temp=None,
              torus_luminosity=None,
              torus_frac=None,
              data_folder=None,
              executable=None,
              command_params_full=None,
              command_params_1=command_params_1,
              command_params_2=command_params_2,
              parameter_file=None,
              prev_files=False,
              use_param_file=False,
              verbose=False,
              eic=False,
              fixed_params=None,
              folder=None,
          )
    
    print("test chi_squared...")
    Chi_squared = blazar_utils.chi_squared(
        params,
        v_data,
        vFv_data,
        err_data,
        name_stem=name_stem,
        theta=theta,
        redshift=redshift,
        min_freq=min_freq,
        max_freq=max_freq,
        data_folder=data_folder,
        command_params_1=command_params_1,
        command_params_2=command_params_2,
        prev_files=prev_files,
        eic=eic,
    )
    print("chi_squared= ", Chi_squared)
    
    plt.plot(logv, logvFv, ls=":", label= "make_model_test")
    plt.errorbar(np.log10(v_data),np.log10(vFv_data), fmt ="o")
    plt.ylim([-16,-9])
    plt.legend()
    plt.show()
    
