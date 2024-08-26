bjet_core
=========

``bjet_core`` contains an updated version of the multi-zone Bjet code developed in `Hervet et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015A%26A...578A..69H/abstract>`_. This is the physical core of Bjet_MCMC and can be run in standalone mode outside the MCMC frame.

Manual installation
-------------------

The code is in the ``bjet_core`` folder. The necessary files are ``bj_core.cpp``, ``bj_core.h``, ``processes_supp_core.cpp``, and ``processes_supp_core.h``.

Compile the executable with 

.. code-block:: console

   $ make bj_core

Usage
-----

Options for calling the executable:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Usage 1: ``bj_core --help``
- Usage 2: ``bj_core <parameter file>``
- Usage 3: ``bj_core 0 <parameter file>`` same execution as usage 2, ``_prev files`` made
- Usage 4: ``bj_core 1 <parameter file>``                           no _prev files made
- Usage 5: ``bj_core 2 <data folder> <parameter file>``            _prev files made
- Usage 6: ``bj_core 3 <data folder> <parameter file>``               no _prev files made
- Usage 7: ``bj_core 0/1 <model type> <list of parameters>``, 0 = yes ``_prev files``, 1 = no
- Usage 8: ``bj_core 2/3 <data folder> <model type> <list of parameters>``, 2 = yes ``_prev files``, 3 = no

Parameter files
^^^^^^^^^^^^^^^

Parameter files specified in the command line should be paths relative to the ``bj_core`` executable. See ``params.txt`` in ``bjet_core`` for an example parameter file.


Modes:
^^^^^^

A `_prev`` file is a file with the previous model SED. When this option is activated, the model will not erase the SED from the former execution of the code. This is especially useful when doing "fit-by-eye"> to see the effects of a change in the model parameters on the same SED plot.
- 0: ``_prev`` files are made, save in default data folder
- 1: no ``_prev`` files are made, save in default data folder
- 2: ``_prev`` files are made, save in specified data folder
- 3: no ``_prev`` files are made, save in specified data folder

(The default data directory is in the same directory as the executable, named "data")

Valid values for model type:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- 0: model with just blob
- 1: blob + EIC parameters
- Other models not yet implemented. The full multi-zone model can be used only with a parameter file.

Order of parameters in command line:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Model type 0:

.. code-block:: console
 
    $ bj_core <0, 1, 2, 3> <data folder if applicable> 0 (model type) z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src, L_src, IIR_level, D_b, NU_DIM, NU_STR, NU_END, prefix
 
argc should be 22 or 23 depending on if data folder is listed

- Model type 1 [*]_:

.. code-block:: console

    $ bj_core <0, 1, 2, 3> <data folder if applicable> 1 (model type) z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src, L_src, IIR_level, D_b, T_BB, TBB_tor, L_nuc, tau, L_tor, tau, NU_DIM, NU_STR, NU_END, prefix

.. [*] Note that tau is present twice, this is a slight error in the bjet code. The second tau value is not used for anything, but it must be inputted.

argc should be 28 or 29 depending on if data folder is listed

- Example:

.. code-block:: console

   $ ./bj_core 3 /Users/sed_calculations 1 0.34 69.6 0.57 50.0 612.1 2.28 3.74 2816.9 1803000 44806 0.00236 5.94e+17 0 1 3.8e+15 2013 2.0e+4 1.7e+21 1.5e-10 5.5e+20 9.0e-5 99 50000000.0 1e+29 run

The 3 indicates that the data folder is specified and no prev file is made. 1 is the EIC model type. Then 0.34 is z (redshift) and the rest of the parameters are enumerated.
