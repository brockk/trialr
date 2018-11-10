
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
