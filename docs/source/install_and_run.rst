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

If you use an updated version of Bjet_MCMC from a previously installed version, you might need to update your conda environment, with

.. code-block:: console 

   $ conda env update --file environment.yml  --prune

Then load your environment with

.. code-block:: console

   $ conda activate bjet-mcmc



Running Bjet_MCMC for the first time
------------------------------------

1. Create a copy of ``mcmc_config_template.txt`` called `mcmc_config.txt`. For information and customized configuration files, see :doc:`configuration_file`.
2. Ensure the data file is well formatted and is at the relative path specified in ``mcmc_config.txt``, see :doc:`data_format`.
3. Install Bjet-mcmc with ``pip install -e .`` (in top directory, where setup.py and pyproject.toml are located). This also compiles the C++ code (see :doc:`bjet_core`).
4. First time use: execute ``blazar_initialize`` with ``python blazar_initialize.py`` or ``python3 blazar_initialize.py`` depending on your Python setup. This creates all necessary folders. If using ``python blazar_initialize.py``, you must be in the directory `bjet_mcmc`.
5. Execute ``blazar_run_mcmc`` with ``blazar_run_mcmc``, ``python blazar_run_mcmc.py`` or ``python3 blazar_run_mcmc.py`` depending on your Python setup.  If using ``python blazar_run_mcmc.py``, you must be in the directory `bjet_mcmc`.
6. Retrieve and check the results (see :doc:`outputs`).

