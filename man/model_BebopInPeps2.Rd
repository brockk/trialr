\name{model_BebopInPeps2}
\alias{model_BebopInPeps2}
\alias{BebopInPeps2}
\title{Stan model for BEBOP implementation in PePS2 clinical trial}
\description{
This is the Stan implementation of the BEBOP design for the PePS2 trial, that
incorporates baseline predictive information to study co-primary efficacy and
toxicity outcomes.

The model is compiled when \code{trialr} is installed and is available at
run-time under \code{stanmodels$BebopInPeps2}.

The design has been submitted for publication by Brock et al. in 2017.

See the BEBOP vignette for a detailed description of the probability model and
a demonstration of the pertinent methods implemented in \code{trial}.
}
\examples{
# Get model parameters as used in the PePS2 trial.
# This call randomly samples patient outcomes so set a seed
set.seed(123)
dat <- peps2_get_data(num_patients = 60,
                     prob_eff = c(0.167, 0.192, 0.5, 0.091, 0.156, 0.439),
                     prob_tox = rep(0.1, 6),
                     eff_tox_or = rep(1, 6))
# Fit the observed data to the model using rstan
samp <- rstan::sampling(stanmodels$BebopInPeps2, data = dat)
# The fit object is quite crude. Post-process to perform useful inference:
decision <- peps2_process(dat, samp)
decision$Accept   # Accept in cohort 2, 3, 5, 6 but not 1 or 4
decision$ProbEff  # The probability of efficacy is driving that decision
}
\seealso{
\code{\link{peps2_params}}

\code{\link{peps2_get_data}}

\code{\link{peps2_process}}
}
