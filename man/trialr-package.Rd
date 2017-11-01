\name{trialr-package}
\alias{trialr-package}
\alias{trialr}
\docType{package}
\title{
  \packageTitle{trialr}
}
\description{
  \packageDescription{trialr}
}
\details{
\code{trialr} is a collection of Bayesian clinical trial designs implemented
in \code{rstan} and R.
Bayesian designs are being used more commonly in clinical trials however
software implementation remains a frequent challenge.
\code{rstan} is the R implementation of Stan, a probabilistic programming
language that performs full Bayesian statistical inference using MCMC sampling.
It allows the statistician to focus on specifying the model (the good bit)
without having to worry about the prior-to-posterior analysis (the nuisance bit).
\code{trialr} seeks to be a showcase of published examples of Bayesian clinical
trial designs, implemented in a common language and style, and brought together
in a single R package.

Clinical trial designs currently implemented:

- EffTox, a design for seamless phase I-II dose-finding, published by Thall &
Cook 2004

- Hierarchical model for a phase II trial of a treatment in a disease with
multiple sub-types using binary response outcomes, published by Thall et al 2003

- BEBOP, a stratified medicine design for studying efficacy and toxicity in
phase II that incorporates predictive baseline information, as developed for the
PePS2 trial and submitted for publication by Brock, et al. 2017

\packageIndices{trialr}
}
\author{
  \packageAuthor{trialr}

  Maintainer: \packageMaintainer{trialr}
}
\references{
Thall, PF, and JD Cook. 2004. "Dose-Finding Based on Efficacy-Toxicity
Trade-Offs". Biometrics 60 (3): 684-93.

Thall, Herrick, Nguyen, Venier, and Norris. 2014. "Effective sample size for
computing prior hyperparameters in Bayesian phase I-II dose-finding".
Clinical Trials 11 (6): 657-66. doi:10.1177/1740774514547397.

Thall, Wathen, Bekele, Champlin, Baker, & Benjamin. 2003. "Hierarchical Bayesian
approaches to phase II trials in diseases with multiple subtypes".
Statistics in Medicine, 22(5), 763-780. https://doi.org/10.1002/sim.1399
}
\keyword{ package }
\keyword{ clinical trial }
\keyword{ Bayesian }
\keyword{ stan }
\keyword{ rstan }
\seealso{

\code{\link{model_EffTox}}

\code{\link{model_ThallHierarchicalBinary}}

\code{\link{model_BebopInPeps2}}

}
\examples{
library(trialr)

# EffTox usage
dat <- efftox_parameters_demo()
# Add some outcomes
dat$num_patients <- 3
dat$eff <- c(0, 1, 1)
dat$tox <- c(0, 0, 1)
dat$doses <- c(1, 2, 3)
# Use rstan to obtain samples from the posterior distributions
samp <- rstan::sampling(stanmodels$EffTox, data = dat)
# Perform inference on posterior samples. E.g. posterior mean drug efficacy
colMeans(rstan::extract(samp, 'prob_eff')[[1]])

## Not run
## Hierarchical model for responses in a disease with multiple subtypes
# dat <- thallhierarchicalbinary_parameters_demo()
# samp <- rstan::sampling(stanmodels$ThallHierarchicalBinary, data = dat)
# rstan::plot(samp, pars = 'p')  # Plot the modelled response rates in the ten subtypes

## BEBOP in PePS2 usage
# set.seed(123)
# dat <- peps2_get_data(num_patients = 60,
#                      prob_eff = c(0.167, 0.192, 0.5, 0.091, 0.156, 0.439),
#                      prob_tox = rep(0.1, 6),
#                      eff_tox_or = rep(1, 6))
# samp <- rstan::sampling(stanmodels$BebopInPeps2, data = dat)
# decision <- peps2_process(dat, samp)
# decision$Accept   # Accept in cohort 2, 3, 5, 6 but not 1 or 4
# decision$ProbEff  # The probability of efficacy is driving that decision
}
