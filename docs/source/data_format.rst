Data format
===========

Data by default is assumed to be in data files with the following columns as labels in the first row::
  
!E(eV)		F(ergcm-2s-1)	delta_E(-)	delta_E(+)	delta_F(-)	delta_F(+)	instrument


For upper limits, set ``delta_F(-)``	and ``delta_F(+)`` at 0. The value set for ``F(ergcm-2s-1)`` should be the flux U.L. at 95% confidence level. Upper limits are considered in the fitting algorithm as likelihood step functions, with a flat probability of 95% below the U.L. and 5% above, extending on both sides to infinity.

Example of an SED data file
----------------------

This is the default data file, already available in Bjet_MCMC : ``real_data/J1010_SED_reduced.dat``.
This dataset of the blazar 1RXS J101015.9-311909 comes from  `Cerruti et al. (2013) <https://www.aanda.org/articles/aa/full_html/2013/10/aa20963-12/aa20963-12.html>`_ ::


  !E(eV)          F(ergcm-2s-1)   delta_E(-)      delta_E(+)      delta_F(-)      delta_F(+)	instrument
  1.5548		8.48E-012	0.17621		0.16843		4.11E-013	4.11E-013	ATOM
  1.9356		8.97E-012	0.16314		0.27997		2.03E-012	3.52E-012	ATOM
  2.8327		1.05E-011	0.35125		0.43236		3.03E-012	5.25E-012	ATOM
  2.2558		1.02E-011	0.18799		0.22558		2.87E-013	2.87E-013	Swift-UVOT
  2.8198		8.63E-012	0.33838		0.44523		1.96E-013	1.96E-013	Swift-UVOT
  3.5449		7.52E-012	0.44311		0.59081		1.83E-013	1.83E-013	Swift-UVOT
  4.0023		9.83E-012	0.90052		1.6373		2.36E-013	2.36E-013	Swift-UVOT
  5.1696		1.30E-011	0.73852		1.0339		6.32E-013	6.32E-013	Swift-UVOT
  5.6396		1.10E-011	0.86763		1.2532		2.14E-013	2.14E-013	Swift-UVOT
  395		9.65E-012	85		85		1.14E-012	1.14E-012	Swift-XRT
  520		1.13E-011	40		40		1.42E-012	1.42E-012	Swift-XRT
  625		7.70E-012	65		65		8.18E-013	8.18E-013	Swift-XRT
  750		7.95E-012	60		60		7.81E-013	7.81E-013	Swift-XRT
  855		7.92E-012	45		45		8.42E-013	8.42E-013	Swift-XRT
  955		6.19E-012	55		55		6.44E-013	6.44E-013	Swift-XRT
  1065		7.15E-012	55		55		7.01E-013	7.01E-013	Swift-XRT
  1185		6.79E-012	65		65		6.37E-013	6.37E-013	Swift-XRT
  1310		7.11E-012	60		60		7.02E-013	7.02E-013	Swift-XRT
  1440		5.75E-012	70		70		6.10E-013	6.10E-013	Swift-XRT
  1585		5.82E-012	75		75		6.36E-013	6.36E-013	Swift-XRT
  1760		5.62E-012	100		100		6.04E-013	6.04E-013	Swift-XRT
  2030		4.69E-012	170		170		5.15E-013	5.15E-013	Swift-XRT
  2420		5.45E-012	220		220		6.26E-013	6.26E-013	Swift-XRT
  3040		4.14E-012	400		400		4.87E-013	4.87E-013	Swift-XRT
  4280		3.41E-012	840		840		4.10E-013	4.10E-013	Swift-XRT
  5820		2.81E-012	700		700		6.79E-013	6.79E-013	Swift-XRT
  1698600000	1.29E-012	698650000	1186800000	0.00E+000	0.00E+000	Fermi-LAT
  4901300000	1.92E-012	2015900000	3424300000	5.59E-013	5.59E-013	Fermi-LAT
  14142000000	1.43E-012	5816600000	9880400000	7.64E-013	7.64E-013	Fermi-LAT
  40806000000	2.27E-012	16783000000	28509000000	1.41E-012	1.41E-012	Fermi-LAT
  117740000000	9.22E-012	48426000000	82259000000	0.00E+000	0.00E+000	Fermi-LAT
  296236003888.79	1.57E-012	88136003888.790	125463996111.20	4.25E-013	4.40E-013	H.E.S.S.
  600355719552.99	5.25E-013	178655719552.99	254344280447.00	2.14E-013	2.22E-013	H.E.S.S.
  1216727524961.9	2.22E-013	362027524961.93	515372475038.06	1.71E-013	1.81E-013	H.E.S.S.
  2465769133556.5	4.30E-013	733669133556.50	1044430866443.5	1.75E-013	1.90E-013	H.E.S.S.
