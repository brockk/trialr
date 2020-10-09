
# Version 0.1.4

## Test environment - local Mac install, R 3.6.2

0 errors | 0 warnings | 2 notes

### NOTEs

❯ checking installed package size ... NOTE
    installed size is 12.6Mb
    sub-directories of 1Mb or more:
      doc    3.8Mb
      libs   8.0Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Both on these notes are OK.

## Test environment - WinBuilder R-release, 3.6.3
TODO

## Test environment - WinBuilder WinBuilder R-devel
TODO

## Test environment - TODO on travis-ci
TODO


# Version 0.1.4

## Test environment - local Mac install, R 3.6.2

0 errors | 0 warnings | 2 notes

### NOTEs

❯ checking installed package size ... NOTE
    installed size is 12.6Mb
    sub-directories of 1Mb or more:
      doc    3.8Mb
      libs   8.0Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Both on these notes are OK.


## Test environment - WinBuilder R-release, 3.6.3

0 errors | 0 warnings | 1 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kristian Brock <kristian.brock@gmail.com>'

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/2531628
    From: inst/doc/CRM-pathways.html
          inst/doc/CRM.html
          inst/doc/EffTox.html
    Status: 403
    Message: Forbidden

Found the following (possibly) invalid DOIs:
  DOI: 10.2307/2531628
    From: DESCRIPTION
    Status: Forbidden
    Message: 403

This URL and DOI are OK.


## Test environment - WinBuilder WinBuilder R-devel

1 errors | 0 warnings | 1 notes

* checking package dependencies ... ERROR
Package required but not available: 'rstan'

Why might rstan be unavailable under R-devel?
The latest version on production CRAN (2.19.3) exceeds the required version (2.18.2).
Perhaps rstan is not available for R v4.0.0?


* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kristian Brock <kristian.brock@gmail.com>'

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/2531628
    From: inst/doc/CRM-pathways.html
          inst/doc/CRM.html
          inst/doc/EffTox.html
          README.md
    Status: 403
    Message: Forbidden

Found the following (possibly) invalid DOIs:
  DOI: 10.2307/2531628
    From: DESCRIPTION
    Status: Forbidden
    Message: 403

This URL and DOI are OK.



# Version 0.1.3

## Test environment - local Mac install, R 3.5.2

0 errors | 0 warnings | 2 notes

### NOTEs

❯ checking installed package size ... NOTE
    installed size is 12.6Mb
    sub-directories of 1Mb or more:
      doc    4.1Mb
      libs   7.7Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Both on these notes are OK.

## Test environment - WinBuilder R-release, 3.6.2

0 errors | 0 warnings | 0 notes

## Test environment - WinBuilder WinBuilder R-devel

0 errors | 0 warnings | 0 notes



# Version 0.1.2

## Test environments
* local Mac install, R 3.5.2
* WinBuilder R-release, 3.6.0

0 errors | 0 warnings | 1 notes

## NOTEs

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kristian Brock <kristian.brock@gmail.com>'

Version contains large components (0.1.1.9002)

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/2531628
    From: inst/doc/CRM-pathways.html
          inst/doc/CRM.html
          inst/doc/EffTox.html
    Status: 403
    Message: Forbidden

This URL is OK.



# Version 0.1.1

## Test environments
* local Mac install, R 3.5.2
* WinBuilder R-release, 3.6.0

0 errors | 0 warnings | 1 notes

## NOTEs

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kristian Brock <kristian.brock@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  CRM (6:65)
  EffTox (6:74)

These words are OK.

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/2531628
    From: inst/doc/CRM-pathways.html
          inst/doc/CRM.html
          inst/doc/EffTox.html
    Status: 403
    Message: Forbidden

This URL is OK.



# Version 0.1.0

## Test environments
* local Mac install, R 3.5.2
* WinBuilder R-devel

0 errors | 0 warnings | 3 notes

## NOTEs

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Kristian Brock <kristian.brock@gmail.com>’

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/2531628
    From: inst/doc/CRM.html
          inst/doc/EffTox.html
    Status: 403
    Message: Forbidden
  URL: https://www.jstor.org/stable/2965714
    From: inst/doc/EffTox.html
    Status: 403
    Message: Forbidden

Both of these URLs are valid. 

* checking installed package size ... NOTE
  installed size is  x.xMb
  sub-directories of 1Mb or more:
    libs   x.xMb
  
  The package contains compiled RStan models, hence the directory size.
  This note was also present in previous CRAN versions.

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

This is true. 



# Version 0.0.7

## Test environments
* local Win install, R 3.5.1
* local Mac install, R 3.5.2

0 errors | 0 warnings | 3 notes

## NOTEs

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Kristian Brock <kristian.brock@gmail.com>’

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/2531628
    From: inst/doc/CRM.html
          inst/doc/EffTox.html
    Status: 403
    Message: Forbidden
  URL: https://www.jstor.org/stable/2965714
    From: inst/doc/EffTox.html
    Status: 403
    Message: Forbidden

Both of these URLs are valid. 

* checking installed package size ... NOTE
  installed size is  7.9Mb
  sub-directories of 1Mb or more:
    libs   6.5Mb
  
  The package contains compiled RStan models, hence the directory size.
  This note was also present in previous CRAN versions.

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

This is true. 

# Version 0.0.6

## Test environments
* local Win install, R 3.5.1
* local Mac install, R 3.5.0

0 errors | 0 warnings | 3 notes

## NOTEs

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Kristian Brock <kristian.brock@gmail.com>’

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.2307/2531628
    From: inst/doc/CRM.html
          inst/doc/EffTox.html
    Status: 403
    Message: Forbidden
  URL: https://www.jstor.org/stable/2965714
    From: inst/doc/EffTox.html
    Status: 403
    Message: Forbidden

Both of these URLs are valid. 

* checking installed package size ... NOTE
  installed size is  8.1Mb
  sub-directories of 1Mb or more:
    libs   6.5Mb
  
  The package contains compiled RStan models, hence the directory size.
  This note was also present in previous CRAN versions.

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

This is true. I understand that other packages in the Stan ecosystem 
(e.g. rstanarm) use the same.


# Version 0.0.5

## Test environments
* local Win install, R 3.5.1
* local Mac install, R 3.5.0

0 errors | 0 warnings | 2 notes

## NOTEs

* checking installed package size ... NOTE
  installed size is  8.1Mb
  sub-directories of 1Mb or more:
    libs   6.5Mb
  
  The package contains compiled RStan models, hence the directory size.
  This note was also present in previous CRAN versions.

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

This is true. I understand that other packages in the Stan ecosystem 
(e.g. rstanarm) use the same.


# Version 0.0.4

## Test environments
* local Win install, R 3.5.1
* local Mac install, R 3.5.0

0 errors | 0 warnings | 2 notes

### NOTEs

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Kristian Brock <kristian.brock@gmail.com>’

Hopefully this NOTE is safe to ignore.

* checking installed package size ... NOTE
  installed size is  8.1Mb
  sub-directories of 1Mb or more:
    libs   6.5Mb
  
  The package contains compiled RStan models, hence the directory size.
  This note was also present in previous CRAN versions.

Furthermore, I am aware of an issue with rstan 2.18 / C++14 (and thus, this package too) on Solaris.
Corresponding with the rstan maintainer, I gather a solution has been proposed but perhaps not yet accepted.




# Version 0.0.3

## Test environments
* ubuntu 14.04.5 on travis-ci, devel fork
* local OS X install, R 3.5.0

## R CMD check results on travis-ci
0 errors | 0 warnings | 1 notes

### NOTEs

* checking installed package size ... NOTE
  installed size is  7.7Mb
  sub-directories of 1Mb or more:
    libs   6.6Mb
  
  The package contains compiled RStan models, hence the directory size.
  This note was also present in previous CRAN versions.





# Version 0.0.2

## Test environments
* ubuntu 14.04.5 on travis-ci, devel fork
* local OS X install, R 3.5.0

## R CMD check results on travis-ci
0 errors | 0 warnings | 1 notes

### NOTEs

* checking installed package size ... NOTE
  installed size is  7.7Mb
  sub-directories of 1Mb or more:
    libs   6.6Mb
  
  The package contains compiled RStan models, hence the directory size.
  This note was also present in accepted version 0.0.1.
  
  
  


# Version 0.0.1

## Test environments
* local Windows install, R 3.4.2
* ubuntu 14.04.5 on travis-ci, devel fork
* local OS X install, R 3.4.2

## R CMD check results on my machine
0 errors | 0 warnings | 2 notes

### NOTEs


* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Kristian Brock <kristian.brock@gmail.com>'

  New submission

This is the first release to CRAN.



* checking installed package size ... NOTE
  installed size is  6.2Mb
  sub-directories of 1Mb or more:
    libs   5.9Mb
    
  The package contains compiled RStan models, hence the directory size.


  

## Additional NOTEs observed in output at incoming_pretest at win-builder
There are none.


  
## Downstream dependencies
There are none.
