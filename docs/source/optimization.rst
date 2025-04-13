Computing Optimization Advice
=============================
Users can configure the underlying MCMC algorithm. Most importantly, users
must define the number of walkers, and the number of steps each walker will 
take. We have tested many different combinations of steps and walkers, to 
come up with the following guidelines. 

Configuration
-------------
- Emcee moves

We are yet to find any computational improvements from modifying the walker 
moves passed to the emcee library. After significant testing, we found the 
default ``StretchMove()`` with parameter ``a=2`` to still be the best.

- Walkers and Steps

Intuitvely, having more walkers and steps will increase computation time. 
More iterations over the MCMC algorithm will result in better results. As to
the ratio between steps and walkers, we suggest having **120 times as many 
steps as walkers.**

Testing
-------
Tests were conducted such that num_steps*num_walkers = 300,000. We observed 
minimum convergence time when we tested 6000 steps and 50 walkers. And find 
our conclusion that 6000/50 = 120 steps/walker. We assume that this ratio 
remains close to optimal as users scale the number of walkers and steps.

 .. image:: ../figures/convergence_time.png
   :width: 800

On the y-axis, :math:`\tau` is the time constant of the exponential decay fitted to the chi squared as a function of steps. The total steps and total runtime are used to convert its units to hours. As a general rule, increasing the number of steps and decreasing the number of walkers increases the total computation time but ensures convergence. 
Decreasing the number of steps and increasing the number of walkers reduces 
the computation time but risks a failure to converge.

OpenHPC - Slurm
---------------
Bjet_MCMC benefits greatly from running on multiple parallel CPUs, and takes 
hours to run. Users may opt to utilize a computing cluster to decrease runtime
and free up their PC. Below is a script that can be used to submit Bjet_MCMC
jobs to the Slurm batch scheduling system used on many compute clusters such
as UCSC's Hummingbird.

.. code-block:: console

  #!/bin/bash 
  #SBATCH --job-name=bjet	#Job Name
  #SBATCH --mail-type=ALL	#What events you want emailed to you (ALL, BEGIN, END, REQUEUE, FAIL)
  #SBATCH --mail-user=<email>	#Email to send job status. Input your email into <email>
  #SBATCH -p 128x24       #Name of which partition is being used 128x24
  #SBATCH --nodes=1       #Number of nodes requested
  #SBATCH --ntasks=1	#Number of tasks requested
  #SBATCH --cpus-per-task=24 #Number of CPUs requested. Make sure that cores requested here matches the amount of cores requested in the config file.
  #SBATCH --output=bjet_%j.out	#Output log. When job runs, output will be written to this file.
  #SBATCH --error=bjet_%j.err	#Error log. When job runs, error will be written to this file.
  #SBATCH --mem=10G         #Amount of memory per node
  #SBATCH --time=5-00:00:00 #the job will be automatically terminated after 5 days if it is still	running

  #Code to run Bjet
  module load miniconda3/3.12 # load miniconda on the node

  # init conda and activate the bjet-mcmc environment
  conda init
  conda activate bjet-mcmc

  # start Bjet_MCMC
  python3 /absolute/path/to/Bjet_MCMC/bjet_mcmc/blazar_run_mcmc.py


