Data format
===========

Data by default is assumed to be in data files with the following columns as labels in the first row:
`!E(eV)		F(ergcm-2s-1)	delta_E(-)	delta_E(+)	delta_F(-)	delta_F(+)	instrument`.
The first row is therefore skipped when reading data.

`E`, `F`, and `delta F(-)` are the columns used by the program, which by default are in columns 0, 1, and 4.
As an optional parameter, these can be changed. The data must be space-delimited.

For upper limits, set `delta_F(-)`	and `delta_F(+)` at 0. The value set for `F(ergcm-2s-1)` should be the flux U.L. at 95% confidence level. Upper limits are considered in the fitting algorithm as likelihood step functions, with a flat probability of 95% below the U.L. and 5% above, extending on both sides to infinity.

All rows of data **must have the same number of rows**.
