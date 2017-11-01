\name{model_EffTox}
\alias{model_EffTox}
\alias{EffTox}
\title{Stan model for EffTox dose-finding design}
\description{
This is the Stan implementation of the EffTox design for dose-finding in
scenarios where dose decisions are guided by efficacy and toxicity outcomes.

The model is compiled when \code{trialr} is installed and is available at
run-time under \code{stanmodels$EffTox}.

The design was first published by Thall & Cook in 2004.
Advice on its use was materially updated in 2014.

See the EffTox vignette for a detailed description of the probability model and
a demonstration of the pertinent methods implemented in \code{trial}.
}
\references{
Thall & Cook. 2004, Dose-Finding Based on Efficacy-Toxicity Trade-Offs

Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample size for
computing prior hyperparameters in Bayesian phase I-II dose-finding
}
\examples{
# Get example parameters from Thall et al. 2014
dat <- efftox_parameters_demo()
# Add outcomes for three patients
dat$num_patients <- 3
dat$eff <- c(0, 1, 1)
dat$tox <- c(0, 0, 1)
dat$doses <- c(1, 2, 3)
# Fit the observed data to the model using rstan
fit <- rstan::sampling(stanmodels$EffTox, data = dat)
# The fit object is quite crude. Post-process to perform useful inference:
decision <- efftox_process(dat, fit)
decision$recommended_dose == 3
# Dose-level 3 is recommended for the next cohort of patients.
}
\seealso{
\code{\link{efftox_params}}

\code{\link{efftox_parameters_demo}}

\code{\link{efftox_process}}
}
