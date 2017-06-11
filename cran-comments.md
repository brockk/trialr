
# Version 0.0.1

## Test environments
* local Windows install, R 3.3.3
* ubuntu 12.04 on travis-ci, devel fork
* local OS X install, R 3.3.2

## R CMD check results on my machine
0 errors | 0 warnings | 3 notes

### NOTEs
* checking installed package size ... NOTE
  installed size is  5.9Mb
  sub-directories of 1Mb or more:
    libs   5.5Mb
  
  The package contains compiled RStan models, hence the directory size.
  

* checking dependencies in R code ... NOTE
Packages in Depends field not imported from:
  'Rcpp' 'methods'
  These packages need to be imported from (in the NAMESPACE file)
  for when this namespace is loaded but not attached.
  
  This package extends rstan. The Rcpp and methods packages are in Depends 
  because this is the action taken by the rstantools package, intended to help
  developers write packages that interface with Stan.


* checking R code for possible problems ... NOTE
efftox_contour_plot: no visible binding for global variable 'dl'
efftox_utility_density_plot: no visible binding for global variable
  'Utility'
efftox_utility_density_plot: no visible binding for global variable
  'Dose'
Undefined global functions or variables:
  Dose Utility dl
  
  Dose, Utility and dl are items plotted by ggplot2. They are not global variables.
  

## Additional NOTEs observed in output at incoming_pretest at win-builder

** running examples for arch 'i386' ... [289s] NOTE
Examples with CPU or elapsed time > 10s
                  user system elapsed
efftox_dtps     107.08   0.32  108.47
peps2_run_sims   68.66   0.02   68.70
efftox_simulate  59.42   0.08   59.92
trialr-package   11.02   0.00   11.04

and

** running examples for arch 'x64' ... [248s] NOTE
Examples with CPU or elapsed time > 10s
                 user system elapsed
efftox_dtps     89.89   0.16   90.09
peps2_run_sims  62.46   0.03   62.53
efftox_simulate 44.02   0.06   44.12
trialr-package  10.88   0.01   10.89

There are some illustrative examples that perform posterior sampling in a small
number of iterations - these take 1-2 mins to run.

* checking compiled code ... NOTE
File 'trialr/libs/i386/trialr.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
File 'trialr/libs/x64/trialr.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
  
These are not functions of mine. There is some boilerplate code added by rstantools
that might explain their existence.

  
## Downstream dependencies
There are none.
