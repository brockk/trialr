
# Version 0.0.1

## Test environments
* local Windows install, R 3.3.3
* ubuntu 12.04 on travis-ci, devel fork
* local OS X install, R 3.3.2

## R CMD check results
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

  
## Downstream dependencies
There are none.
