#CHANGELOG
###last updated: 2022-05-30

**2022-05-30** (27ca2c1) HUGE ERROR found. When make_model called
process_model, it did not specify EIC, making the EIC component do
practically nothing. This was really bad.     

     Plotting functions in plot_data_with sources modified to make
     results for the paper.

**2022-05-28** (604e45e) Fixed big in plotting chi squared, added
options for parsing past info.txt files and processing backends.     

     Previously, the plotting the chi squared directly modified the
     array instead of making a copy, which meant that the wrong plots
     were created. info.txt and basic_info.txt files can be parsed to
     read the configurations. This combined with h5 files means that
     data can be re-loaded with the configurations used when the code
     was run. There is now a save function in blazar_report that allows
     just inputting a folder instead of having to directly pass all the
     function results like indices within one sigma. Code was written in
     plot_data_with_sources.py to plot the different sources of the
     3C66A block 1 data in different colors. (This will probably work
     for other blocks too, but I just eyeballed the number cutoffs for
     values that worked with the data, so this may not be true.)

**2022-05-25** (31a1fe6) Delete Bjet_02_2022_03_28.tar.xz     


**2022-05-25** (fcad12e) Add files via upload     


**2022-05-24** (d9e03b7) Changelog updated     


**2022-05-24** (b8ac193) mkdtemp used instead of TemporaryDirectory to
prevent being deleted prematurely with nohup. blazar_run.sh runs
caffeinate only on Darwin and not Linux.     


**2022-05-19** (81ba653) blazar_run.sh runs the command to run the
program and saves the process number. Output of compiling bjet is now
suppressed by default in blazar_initialize     


**2022-05-19** (6c3ceaa) 3C66A blocks 2, 5, 6 data added     


**2022-05-19** (9e37a3e) Removed excess error messages from Makefile     


**2022-05-19** (943c2b1) TMP setting now in blazar_properties. This
makes files saved in temp directory.     


**2022-05-13** (9c4123c) Make compression at end of mcmc cross-platform     


**2022-05-13** (575ba4f) Changed USE_BASE_PATH back to false to prevent
problems     


**2022-05-13** (56295ff) Add sed_calculations to gitignore     


**2022-05-13** (3c545cd) Descriptive names to results files, remove h5
and pickle files now stored on google drive to decrease size     


**2022-05-13** (f50913e) Remove circular dependency in blazar_model and
blazar_utils, change Hubble constant     

     blazar_model now contains the log_to_linear and the linear_to_log
     functions so that it does not need to reference blazar_utils.

     Hubble constant is now 69.6, not 71.0

**2022-05-13** (b7c739c) EIC results for 3C66A added from VHE     


**2022-05-07** (c3db53b) mcmc_config is now in .gitignore,
mcmc_config_template.txt is added as a guide     


**2022-05-07** (931f663) Fixed read_configs to allow for no
tau_variability provided     


**2022-05-07** (d3a3f4f) Add eic to configs file     


**2022-05-07** (ccd76d0) Merge branch 'main' of
github.com:sarah829/blazars-mcmc into main     

     branch from vhe

**2022-05-07** (f74bc2e) requirements.txt from pip environment created     


**2022-05-07** (37ae0b8) Fix error in plot_1sigma, add redshift and
use_variability to configs, add redshift as optional param to relevant
functions, optional parameters added     

     In get_params_1sigma_ranges, minima and maxima started
     as all copies of flat_samples[0], which was a problem because that
     is not necessarily within 1 sigma. It now starts with
     flat_samples[indices_within_1sigma[0]], which fixes the error. 
     Redshift is now in the parameter file and is passed into the
     log_probability function, so it was added to the chi_squared
     function, etc. Having a variability constraint is now optional and
     set by use_variability in the config file. Functions now have more
     of the optional parameters that the functions they call have,
     allowing for more flexibility. Code was re-formatted.

**2022-05-05** (fe4d13d) Typo fixed (added space before Todo)     


**2022-05-05** (017a3cd) Changed name stem used in make_model in
blazar_report, fixed typo in README     

     make_model in blazar_report now uses run<random number> for the
     file name stem, and the files are deleted after use. Additionally,
     it no longer uses a parameter file (now the MCMC code never uses
     one).

**2022-05-05** (a12213d) Add todo list to README, update name stems in
blazar_plots, add 3C66A data     

     Name stems in blazar_plots are now for_plotting<random number>
     instead of just for_plotting. This allows for running multiple runs
     of the program at the same time without problems. The functions now
     delete the created SED models after running.

**2022-04-27** (ce0b76f) readme updated to reflect changes in file
organization     


**2022-04-23** (ade3ddc) References to blazar_eic.py removed     


**2022-04-22** (81dc301) Create local_results and old directories and
add to .gitignore     

     local_results contains results files that don't need to be in the
     public repo, and old contains now defunct code and results. Files
     in .gitignore were dropped.

**2022-04-22** (7f65dd8) File reorganization starting to clear out old
files (results, defunct test scripts, etc.) from the repo     


**2022-04-22** (2ace9f0) Merged eic branch to main     


**2022-04-21** (bf80b56) Removed folder for gitignored backups     


**2022-04-21** (b48dafa) EIC code now integrated with regular code
instead of in blazar_eic.     

     Whether or not the model is an eic model is stored in a boolean eic
     for most methods. Additionally, in main of blazar_run_mcmc, the h5
     and pickle files are zipped.

**2022-04-03** (1664363) MCMC with EIC now works, BlazarProperties class
added.     

     blazar_properties changed to now have a BlazarProperties class.
     blazar_properties has a function modelProperties(is_eic) that will
     return a BlazarProperties object that either has SSC values or EIC
     values. Now to call for values, you can use
     modelProperties(eic).NUM_DIM, for example. This creates flexibility
     and makes it a lot easier to incorporate EIC into the main code.
     Updates were made throughout the code to make it consistent with
     this. MCMC was successfully run with EIC.

**2022-04-03** (e66092d) Code for adding more than just synchrotron and
inverse Compton data (add_data in blazar_model).     

     The code will line up and interpolate the existing data to add in
     new data for a given file suffix. Code for testing this is in the
     jupyter notebook blazars.ipynb and test_eic.py. It creates an
     interactive plot with different components added.

**2022-03-17** (8dcaeb9) Bug fixes in calling the subprocess and getting
the parameters     


**2022-03-16** (910a6ec) plots made for presentation     


**2022-03-09** (030f2cb) Integrating EIC     

     Made a copy of functions so that they could be used with EIC. There
     are now 13 parameters; the original 9 + disk black body
     temperature, luminosity of the blob, tau (fraction absorbed), and
     the distance of the the blob from the observer.

**2022-03-08** (bb8efec) Added EIC to cpp, fixed bug with cores in
config file     

     Now, code can be called to set EIC parameters using command line
     arguments. In the config file, there was #cores=`<num cores>` 
     instead of cores = `<num cores>`, which made the number of cores
     not be read and the default always be used.

**2022-02-26** (aa45c74) Folders for MCMC now set in
blazar_properties.py, SED creation w/out param file sped up, started
unit tests.     

     In blazar_properties.py, the folders for parameter files, results,
     and data are set, to make the program more flexible.

     In make_sed_no_file, previously, the list of 23 parameters was
     generated each time. Now, the list for before the custom parameters
     and the list for after the parameters are passed to the function.
     These are passed to log_probability, which passes to chi_squared,
     which passes to make_model, which calls make_sed_no_file. The
     original list is created in blazar_run_mcmc.py. It is now
     significantly faster.

     blazar_mcmc.py renamed to blazar_run_mcmc.py

     test_blazar_mcmc.py currently only has tests for blazar_clean.

     Documentation added to blazar_model.py but not finished.

**2022-02-25** (933111b) Create variables for number of dimensions and
parameter names, add Python scripts to initialize and clean before/after
mcmc.     

     Previously, the number of dimensions was hard-coded at 9 and the
     parameter names were hard-coded. They are now in variables in
     blazar_properties.py. This is to make future changes to the code
     easier. Most functions now function regardless the number of
     dimensions and the parameter names. The exceptions are the
     functions for creating parameter files and making models that rely
     on the exact parameters in the right places.

     blazar_initialize.py creates necessary folders and compiles code. 
     It creates folders for:
     - Data (by default, folder named data)
     - Mcmc results (by default, folder named mcmc_results)
     - If specified by a parameter that parameter files are being used
     (default is false), for parameter files (by default, folder named
     parameter_files) It will also compile the bjet code, removing
     previous object files and compiling.

     blazar_clean.py exists to delete all files in data folder and all
     files in parameter file folder except for params.txt if present.

     blazar_initialize.py and blazar_clean.py should be system
     independent because all commands are used with the python os module
     except for invoking `make` for the executable, which can handle
     different platforms.

**2022-02-25** (6fc1015) Integrated new C++ code, moved Python source
into module, switched to absolute file paths     

     File paths now use the absolute path by finding the parth to
     blazars-mcmc and using it as the bath path in script
     blazar_path.py.

     The c++ code, bjet_mcmc, has changed from bjet:
     - Option to not create prev files
     - Option to specify the data folder
     - Option to enter parameters via the command line (biggest change)
     - Code has been moved into functions. The main code is now in
     run_models(), and there are functions to set variables from a list
     while checking error bounds for all blob parameters. In main, the
     input from the command line is parsed and the proper variables are
     set and the proper way to load the variables (parameter file or
     command line) is called.

     The new C++ code is now being used in the Python code. The Python
     code now has the option to have prev files or not and to use a
     parameter file or not, both defaulting to True.

     There is now an option to add a custom message to the top of the
     info.txt doc in the results for the MCMC.

**2022-02-23** (b061f23) C++ code compiled     


**2022-02-23** (6f14c49) Large plot removed from jupyter notebook     


**2022-02-23** (e710474) fixed C++ executables     


**2022-02-23** (63f11ff) Dependencies added to README     


**2022-02-23** (8c3230c) Add calling C++ with a bunch of command
arguments for the parameters.     

     I added a new function in the C++ code and restructured it a little
     bit. Then, running it from Python is very easy. This is 18–19%
     faster than calling the C++ code with a parameter file while the
     previous method (pybind11) was only 4% faster or basically the
     same.

**2022-02-22** (571ba2f) Integrate C++ code using pybind11     

     Createdthe C++ bindings and the calls for linking and compiling the
     C++ code. Create test_run_cpp.py to test running making the model
     to make sure it works and compare it to the existing make_SED.
     Benchmarking started. Modified c++ code was added, now in the main
     folder (planning to change later). This code makes it so that it
     does not create prev files and it has the functions for inputting
     parameters with a function call.

**2022-02-18** (8857c2c) data files changed     


**2022-02-18** (8bba727) changelog updated     


**2022-02-16** (6593196) Documentation (module docstring and function
docstrings) completed for blazar_utils and blazar_mcmc     


**2022-02-16** (cd6e59e) docstrings explaining functions in the module
completed for blazar_mcmc, in progress for blazar_utils     


**2022-02-16** (65e94f8) added message in main docstring with order of
params     


**2022-02-13** (3ba1d68) increased documentation     


**2022-02-11** (3187c7f) images moved to presentation images     


**2022-02-11** (fe48250) new results added     


**2022-02-10** (5100981) new plots from run_2022-02-09-00:42:26     


**2022-02-09** (a28c1ff) Documentation additions in utils, added cores
to configs     


**2022-02-09** (aadbbbe) integrated config file changes, updated readme
with instructions     

     read_configs updated and configs now correctly referenced
     throughout program. README has basic instructions for use, though
     more are necessary.

**2022-02-08** (278889f) Fixed normalized chi squared calculation     

     Previously, I divided by the number of data points, not the number
     of data points minus the number of free parameters (9). This is now
     fixed.

**2022-02-08** (238f12e) Methods moved, mcmc method created, config file
being modified.     

     calculations.py no longer exists and blazar_reports.py is now here
     with most of the methods from blazar_mcmc.py besides the MCMC
     methods. All plotting is in blazar_plots.py. The Python file names
     are now blazar_mcmc.py, blazar_utils.py, blazar_plots.py,
     blazar_report.py, and blazar_model.py for consistency.

     The config file has had some removed true values and file name
     stem, alpha 2 limits are now all on one line. Added parallel
     processing, tau_variability value, and discard number These changes
     HAVE NOT been reflected in the read_configs method or in the
     functions.

     The program CURRENTLY DOES NOT WORK. (Yeah... in retrospect,
     definitely should've made a new branch.)

**2022-02-08** (8aeac92) Added text-based printout of best param values
and min/max 1 sigma converted into base 10     


**2022-02-06** (bcf1f7a) Updated corner plot, made image type flexible,
using svg for images     

     On the corner plots, there are now intersecting horizontal and
     vertical lines showing the best values in each of the plots. Above
     each of the histograms, the best value with the +/- within 1 sigma
     is shown. A title was added. Tick labels are now smaller,
     preventing axis labels from overlapping with tick labels.

     In plotting_and_images, there is a global variable image_type
     currently set to "svg" which controls what type all images are
     saved in. svg is now being used over jpeg because vectorized images
     allow for clearer plots. TODO: put the image type in the config
     file (along with discard number and whether parallel processing
     should be used + number of cores)

**2022-02-06** (2f9a11f) changed function names of error plots in
blazar_mcmc     


**2022-02-04** (c0646d3) Combined plotting random models within 1 sigma
w/ plotting the extreme models     

     Now, there is one method called plot_with_error_filled in
     calculations.py. It has an optional parameter setting extreme or
     not, which defaults to True. This determiens what kind of plot is
     generated.

     A function to find the min and max for each point value given the v
     values and a list of parameters was created, and this can be used
     in both cases.

**2022-02-04** (e4c1c48) Added tau variability constraint     

     The tau variability constraint is that tau_var >= (1 + redshift) /
     c * R / delta. This is added to the prior distribution function
     which has the other checks that a set of parameters is valid.

**2022-02-04** (461138d) Testing extrapolation now makes graphs in one
plot     

     Using an AxesSequence class from Joe Kington on StackOverflow in
     2012
     https://stackoverflow.com/questions/13443474/matplotlib-sequence-of-figures-in-the-same-window. 
     AxesSequence makes it so that several plots can be displayed in one
     window by using the arrow keys to navigate between them. This is
     used for the testing extrapolation, and now all of the plots with
     extrapolated points are in one plot.

**2022-02-04** (1bd4c91) Added progress bar to random plot (takes a long
time)     


**2022-02-04** (37d5428) Tests and convenient methods to run plot making
commands easily in testing.py     


**2022-02-04** (b966e0f) Added labels and legends to the plots     


**2022-02-04** (71f683b) Images switched from png to jpeg     


**2022-02-04** (823d98f) new data and results     


**2022-02-04** (29f72ad) Added random number to unique file names     

     Unique file names now are random number up to 2^30 + the first 50
     chars of the params This is because I saw that there were many
     duplicate parameter values in the end results of the MCMC, making
     naming the files with the parameters extremely dangerous.

**2022-02-03** (6f67ec0) Switched chi squared calculation to
interpolating with log values     

     Now, the interpolation from model frequency to model flux that real
     flux values are plugged into interpolates from log_v to log_vFv.

     In testing.py, there is testing under the extrapolation testing.
     The two chi squared values generated from the old method and the
     new method are similar enough to seem reasonable. I also added
     displaying the chi squared that would be gotten if data points were
     dropped, not extrapolated (must remove from program for other
     uses).

**2022-02-03** (f969738) Added plotting with models from extreme param
values     

     The parameters within 1 sigma that have the biggest and smallest
     values for each parameter are found, resulting in 2 arrays of
     dimension ndim * ndim. Models are created from these, and for each
     frequency value, the minimum and the maximum are found. The graph
     is made by filling in the space between the minimum and maximum for
     each frequency value. The best model and the actual data with error
     bars are plotted on top of this.

**2022-02-03** (54735e9) Changed error in plot_with_error, fixed corner
plot     

     The error is now calculated with the differences in logs, and upper
     limits were added. The problem with the corner plot was that in
     min_max_1sigma calculated the max by setting values as max if they
     were > min, not > max. It's fixed now and it works.

**2022-02-02** (56318ec) fixed weird problem, merged branches back     


**2022-02-02** (94c0c94) Turn code in calculations.py into functions,
add all plots in plotting and call in blazar_mcmc     

     Functions to plot chi squared, now can plot all, average, or best. 
     Function to plot real data with error bars and plot lines from
     model within 1 sigma

     Note: this is on a branch, rolled back a commit because something
     went wonky (log_prob vals were in the millions for no conceivable
     reason, never figured out the cause).

**2022-02-01** (331ff5a) remove old data     


**2022-02-01** (1698246) calculations.py moved into functions, all plots
added to blazar_mcmc save_plots_and_info     


**2022-01-27** (077758d) Change the order of parameters everywhere,
remove src directory     

     Parameters are now in the following order:
     [delta (linear), K (log), n1 (linear), n2 (linear), gamma_min
     (log), gamma_max (log), gamma_break (log), B (log), R (log)]

     They used to be in the following order:
     [R, B, delta, gamma_break, n1, K, gamma_min, gamma_max, n2]

     The new order is the same order that the parameter file uses for
     calling the c++ code, and all in all, this order makes a lot more
     sense. It has not been fully tested yet.

     ---

     src directory was causing problems because it wasn't added to path
     making it hard to import files from it, so all python code is back
     in the main directory.

**2022-01-27** (78231a1) removed extraneous files     


**2022-01-27** (a35a9dd) moved python files to src     


**2022-01-26** (af1ec84) Add plotting avg chi squared or best chi
squared by step     

     plot_chi_squared now has a "type" argument that can be 'avg',
     'best', or 'all', which controls what kind of plot of chi squared
     values are made.

**2022-01-26** (1f1919e) Integrate parallelization into main
blazar_mcmc; add setting min_freq and max_freq in make_model     

     blazar_mcmc_parallelization.py deleted and now parallel processing
     is a parameter in blazar_mcmc.

     min_freq and max_freq are now changed if the real data is outside
     the range of the defaults.

     position renamed to values in corner plotting (misleading since
     input is full chain)

**2022-01-26** (d93e5a9) Calculations for error added, added parameters
in plotting_and_images     

     calculations.py is currently a rough draft--it has no functions and
     is just code. It calculates chi squared values and calculates the
     median value of each parameter and the 16th and 84th percentiles.
     It also then calculates error by using all samples (set of
     parameters) such that chi^2 < chi^2_min + delta chi^2 for alpha =
     .32.

     It creates a plot with the lines from models of a random sample of
     the parameters meeting the above criterion with the best model line
     darker and in front and with the actual data. This image is saved
     to a png and pickled for later use.

     The generation of the plot is slow since it has to generate the
     model for every sample. I made plots initially with only 100 of the
     plots (out of over 30,000--over half of the plots met the
     criterion), and I am running a bigger sample overnight.

**2022-01-24** (192e8ab) Changed methods to be compatible with parallel
processing     

     In blazar_mcmc_utils, there is now an option to generate a file
     name for parameters and for data that comes from the parameters
     passed to it turned into a string. This is so that there are unique
     files for the parallel processing.

     When this option is in place, the files are deleted after each run,
     and this was added to the log_prob function in blazar_mcmc_utils.

     blazar_mcmc_parallelization was created, which is a copy of
     blazar_mcmc. Duplicating this file was not necessary as most
     changes for parallelization were made in blazar_mcmc_utils, but it
     has the code for running with parallel processing.

     The program was run with parallel processing first with 5,000 steps
     (in 2:45:..) and then 25,000 steps (in a little over 12 hours).

**2022-01-24** (63af187) add p0 to info.txt     


**2022-01-22** (9efcf63) Add functionality to scale plot to real data,
plot log_prob, renamed variables and added documentation     

     Function to get the min and max for scaling to real data added to
     plotting_and_images. All plotting functions now have a parameter
     for adjust_scale, along with an adjust_multiplier for what the
     scale should look like.

     plotting log_prob added to plotting_and_images. This also allows
     for scaling to remove outliers. (Otherwise the plot is very hard to
     see.)

     actual_* renamed to real_* for consistency, documentation for added
     parameters added, and code refactored to avoid too long lines and
     get rid of warnings

     additional runs of the mcmc and the random testing added more data
     files.

**2022-01-22** (5fe95d7) updated changelog     


**2022-01-22** (c9018fe) removed old practice data files     


**2022-01-22** (afea1d3) blazar_mcmc functionalities: automatically
saving plots, loading initial state from data     

     function that makes an info txt, a corner plot, a SED with the best
     model and the actual data on it, and a plot of the chain created

     p0 can now come from the last state in a previous run, effectively
     allowing it to continue.

     plotting functions were moved to plotting_and_images.py; they are
     just called from blazar_mcmc.

**2022-01-22** (61bf282) functions removed from random_testing and
make_model and put into plotting_and_images     


**2022-01-22** (03aedfe) plotting_and_images created, functions related
to plotting moved here     

     plotting functions from random_testing, make_model, and blazar_mcmc
     moved here

**2022-01-21** (ed2fdfb) changed chi_squared to extrapolation     

     Changed chi_squared function from discarding actual data that was
     out of range of the model to extrapolating from the data. In the
     future, I will try not to have data out of range, so this is a
     temporary fix.

     Switched likelihood function from -.5 * chi_squared / num_samples
     back to -.5 * chi_squared since now there are always the same
     number of data points again.

**2022-01-21** (c8215bf) Ran MCMC, added timing and info display to
mcmc, added git lfs     

     Ran MCMC overnight, so the results files are added, as are the
     images of the plots for the best results found created after.

     In blazar_mcmc, the function now times the MCMC run and records the
     time in the file with configs and info about the sampler.

     A function was added to blazars_mcmc to print salient information
     after MCMC run.

     The return type of the main MCMC function was switched to returning
     a tuple of the sampler itself and the time elapsed so that the main
     function can access properties of the sampler directly.

     Git lfs was enabled for *.h5 files (used to store sampler chains)

**2022-01-21** (fe3248e) blazar_mcmc functionality added: file saving,
plotting     

     data saved in backend in h5 files during run, and a function
     fetching data from the h5 file and writing information into a text
     file was created. This text file lists configurations, the best
     parameters and log_prob, all log_probs, autocorrrelation times, and
     the entire chain.

     a function plot_chain was added. this functionality was moved from
     the main run_mcmc method

     this is the first version of the mcmc post-unmessing up the
     likelihood function, so the code is starting to actually work.

**2022-01-21** (d3bd989) make print only when verbose in
blazar_mcmc_utils     


**2022-01-21** (db7467c) fixed log_prob in blazar_mcmc_utils     

     isfinite switched to 'not isfinite'. This... explains a lot of the
     problems seeing that previously it was only running the simulation
     when it shouldn't. Things work a lot better now.

**2022-01-20** (7ba236e) data from additional runs     


**2022-01-20** (0bfd77f) Added info about the order parameters are
listed in     


**2022-01-20** (7bb0d4c) Changed likelihood function     

     Trying averaging over the number of data points to compensate for
     removing data points that are outside of the interpolation range
     (see previous commit). Unsure if this is a good idea.

**2022-01-20** (91d9a45) Changed likelihood function and added saving to
file     

     Likelihood function is now using the avg over the number of data
     points. I don't know if this is a good idea. It is related to
     removing data points from the real data when they are outside of
     the range of the interpolation. I worry that by doing this, it will
     mess up the chi squared values, and averaging over the number of
     data points may help.

     This would not be necessary if there were another solution to
     dealing with points out of range, and I am not sure if just
     deleting them is the best idea.

**2022-01-19** (d38f9f6) typo in README     

     ##Formatting... to ## Formatting...

**2022-01-19** (1bceb3c) README updated, minor formatting changes     


**2022-01-19** (0c67223) old data files removed     


**2022-01-19** (a614d2b) fix refactoring error     

     forgot to change function references to blazar_mcmc_utils from
     mcmc_utils

**2022-01-19** (8137c99) documentation updating and refactoring     

     Improved documentation, split blazar_mcmc_9_params.py into
     blazar_mcmc_utils.py and blazar_mcmc.py

**2022-01-19** (b963d16) Documentation edited in make_model     


**2022-01-17** (3858f08) minor wording changes     


**2022-01-17** (caa0dd1) created to make changelog from git commits     


**2022-01-17** (2ecd9c9) exception handling for cd failing     


**2022-01-17** (cac0ac6) attempted mcmc initial implementation     

     currently not working, says the chain is too short for
     autocorrelation time

**2022-01-17** (6348763) process_model combining compton and synchrotron
data modified     

     load data so that it checks if all lines are the correct length,
     check if synchrotron is ever > max compton and the other way (had
     some issues with this)

**2022-01-12** (02cb1b5) jupyter notebook created     


**2022-01-12** (3324bfc) testing now plots successive best models and
saves as images, both for real data and fake model data     


**2022-01-11** (ea46d68) Testing uses random parameters to find the best     

     It plots it and saves the best ones as images and a gif

**2022-01-06** (2768b82) make process_model use interpolation     

     Previously, the code looked for synchrotron and compton frequency
     values that were exactly equal, which could cause all kinds of
     problems. Now, the program approximates compton vFv values using
     interpolation and does not require exactly equal frequencies.

**2021-12-28** (40aa282) Log prob functions and parameter validation
added to mcmc file     


**2021-12-18** (3b3294d) started work on mcmc code, read params and data     


**2021-12-18** (81ea148) plotting added to makemodel.py     


**2021-12-17** (5b17a62) makemodel.py mostly done; v & vFv data
processing     


**2021-12-17** (2700445) file reorganization, fixed param format so
create_params_file works     


**2021-12-16** (0b1c70a) Initial commit: Bjet_02 added     


**2021-12-16** (6b97964) Initial commit     


