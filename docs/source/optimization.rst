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

