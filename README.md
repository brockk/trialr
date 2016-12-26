# trialr - clinical trial designs in R &amp; Stan

This R package implements clinical trial designs in R via RStan.
RStan is the R connector of the [Stan](http://mc-stan.org/) project.

Stan offers full Bayesian statistical inference so the designs implemented here tend to be Bayesian in nature.
Furthermore, there is a preponderance of early-phase clinical trial
designs because a) these are the designs that tend to use Bayesian methodology; and b) this is my research area.

Functions are provided to invoke a trial analysis on observed outcomes and perform simulations.

Trial designs currently implemented include:
- EffTox, by Thall et al. 
- Hierarchical Bayesian model for binary outcomes, by Thall et al. 

If there is a published Bayesian design you want implemented in Stan, get in touch.
Contact @brockk on GitHub
