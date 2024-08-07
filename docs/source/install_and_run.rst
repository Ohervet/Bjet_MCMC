Installation and running Bjet_MCMC
==================================

.. _installation:

Installation
------------
You first need to clone the Bjet_MCMC Github repo in your local computer with

.. code-block:: console

   $ git clone https://github.com/Ohervet/Bjet_MCMC

Ensure all dependencies are installed, see :doc:`dependencies`. Recommended: create the conda env "bjet-mcmc" from ``environment.yml`` using

.. code-block:: console 

   $ conda env create -f environment.yml

Then load your environment with

.. code-block:: console

   $ conda activate bjet-mcmc



Running Bjet_MCMC for the first time
------------------------------------

1. Create a copy of ``mcmc_config_template.txt`` called `mcmc_config.txt`. For information and customized configuration files, see :doc:`configuration_file`. 
2. Ensure the data file is well formatted and is at the relative path specified in ``mcmc_config.txt``, see :doc:`data_format`.
3. Navigate inside folder ``blazar_mcmc``. Currently, code **must** be run from inside this folder or the imports won't work.
4. First time use: execute ``blazar_initialize.py`` with ``python blazar_initialize.py`` or ``python3 blazar_initialize.py`` depending on your Python setup. This creates all necessary folders and compiles the C++ code (see :doc:`bjet_core`).
5. Execute ``blazar_run_mcmc.py`` with ``python blazar_run_mcmc.py`` or ``python3 blazar_run_mcmc.py`` depending on your Python setup. 
6. Retrieve and check the results (see :doc:`outputs`).

