Outputs
=======
.. _outputs:

By default, all outputs are written to a folder in ``local_results/`` entitled ``folder_label_yyyy-mm-dd-hh:mm:ss``

``basic_info.txt``
------------------
Contains information from the configuations file, some paths and the running time.
e.g. ::

 folder name: local_results/J1010_2023-07-04-23:03:45
 report description: J1010

 config file: mcmc_config.txt
 prev_files: False, use_param_file: False

 configurations:
 {'description': 'J1010', 'folder_label': 'J1010', 'eic': False, 'data_file': 'real_data/J1010_SED_reduced.dat', 'n_steps': 5000, 'n_walkers': 100, 'discard': 200, 'parallel': True, 'cores': 15, 'use_variability': True, 'tau_variability': 24.0, 'redshift': 0.143, 'custom_alpha2_limits': False, 'bb_temp': 'null', 'l_nuc': 'null', 'tau': 'null', 'blob_dist': 'null', 'alpha2_limits': [1.5, 7.5], 'fixed_params': [-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf]}

 p0: random

 time: 5:45:22.151980

``info.txt``
------------
Main output containing the best fit parameters with associated errors, and statistical information.

