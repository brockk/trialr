
# trialr 0.1.5
Added the Neuenschwander, Branson & Gsponer model for dose-escalation. Also
added a method for calculating EffTox parameter priors. Re-ordered vignettes
to be more intuitive. Suppressed all the MCMC log messages in the tests scripts.


# trialr 0.1.4

Updated dependencies to account for recent breaking changes in tibble v3.0.0.


# trialr 0.1.3

Adjusted uses of tidyr::unnest to suppress warnings when using v1.0 of tidyr.

# trialr 0.1.2

Added the AugBin model for phase II response assessment in cancer. Just the 
two-stage single-arm version added for now. The others will follow. Also added
a general routine for running simulation studies.

# trialr 0.1.1

Plumbed in support for tidybayes. Added pathways analysis for CRM and rewrote
the same for EffTox.

# trialr 0.1.0

Added TITE-CRM implmentation plus tests and vignette.
Rebased package to the format now advocated by rstantools.

# trialr 0.0.7

Making package work with staged installation, with help from Tomas Kalibera.

# trialr 0.0.6

Fixing some duplicated vignette titles

# trialr 0.0.5

Adding changes advised by rstan maintainer so that package may build on Solaris.

# trialr 0.0.4

Updated to use rstan 2.18.1, which in-turn has been updated to use C++14.

# trialr 0.0.3

This release updates some Stan code that was generating C++ code that would
compile on Linux, Mac & Windows, but not Solaris.

# trialr 0.0.2

This release adds the Continual Reassessment Method (CRM) for dose-finding.
Four model variants are currently given:
- empiric likelihood with normal prior on slope;
- logistic likelihood with constant intercept and normal prior on slope;
- logistic likelihood with constant intercept and gamma prior on slope;
- logistic likelihood with normal priors on the intercept and slope.

This release also adds general model-fitting functions stan_crm and stan_efftox.
It also adds S3 classes crm_fit and efftox_fit with the expected generic methods.

Unit tests have been added.

# trialr 0.0.1

First release with implementations of:
- Thall & Cook's EffTox dose-finding clinical trial design;
- Thall et al.'s hierarchical Bayesian phase II design for diseases with multiple subtypes
- Brock et al.'s BEBOP design for efficacy and toxicity outcomes in phase II where predictive information
