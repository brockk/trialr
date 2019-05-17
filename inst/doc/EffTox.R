## ---- message=FALSE, results = "hide"------------------------------------
library(trialr)

outcomes <- '1NNE 2EEB'
fit <- stan_efftox_demo(outcomes, seed = 123)

## ---- caption = 'ProbEff and ProbTox are the probabilities of efficacy and toxicity at each dose. ProbAccEff is the probability that the efficacy rate exceeds the desired threshold and ProbAccTox the probability that the toxicity rate is less than the threshold. All probabilities are posterior means.'----
fit

## ---- results = "hide"---------------------------------------------------
fit <- stan_efftox(outcomes,
                   real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                   efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                   p_e = 0.1, p_t = 0.1,
                   eff0 = 0.5, tox1 = 0.65,
                   eff_star = 0.7, tox_star = 0.25,
                   alpha_mean = -7.9593, alpha_sd = 3.5487,
                   beta_mean = 1.5482, beta_sd = 3.5018,
                   gamma_mean = 0.7367, gamma_sd = 2.5423,
                   zeta_mean = 3.4181, zeta_sd = 2.4406,
                   eta_mean = 0, eta_sd = 0.2,
                   psi_mean = 0, psi_sd = 1, 
                   seed = 123)

## ------------------------------------------------------------------------
fit$recommended_dose

## ---- fig.width = 7, fig.height = 7, fig.cap = "Utility contours after observing outcomes 1NEN 2NBE."----
efftox_contour_plot(fit)
title('EffTox utility contours')

## ---- fig.width = 7, fig.height = 7, fig.cap = "Utility densities after observing outcomes 1NEN 2NBE."----
efftox_utility_density_plot(fit, doses = 1:3) +
  ggplot2::ggtitle("EffTox dose utility densities")

## ------------------------------------------------------------------------
knitr::kable(efftox_superiority(fit), digits = 2, row.names = TRUE)

## ------------------------------------------------------------------------
p <- efftox_solve_p(eff0 = 0.5, tox1 = 0.65, eff_star = 0.7, tox_star = 0.25)
dat <- list(
  num_doses = 5,
  real_doses = c(1, 2, 4, 6.6, 10),
  efficacy_hurdle = 0.5,
  toxicity_hurdle = 0.3,
  p_e = 0.1,
  p_t = 0.1,
  p = p,
  eff0 = 0.5,
  tox1 = 0.65,
  eff_star = 0.7,
  tox_star = 0.25,

  alpha_mean = -7.9593, alpha_sd = 3.5487,
  beta_mean = 1.5482, beta_sd = 3.5018,
  gamma_mean = 0.7367, gamma_sd = 2.5423,
  zeta_mean = 3.4181, zeta_sd = 2.4406,
  eta_mean = 0, eta_sd = 0.2,
  psi_mean = 0, psi_sd = 1,

  doses = c(),
  tox   = c(),
  eff   = c(),
  num_patients = 0
)

## ----run_sims, eval = FALSE, cache = TRUE, results = "hide"--------------
#  set.seed(123)
#  sims = efftox_simulate(dat, num_sims = 100, first_dose = 1,
#                         true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
#                         true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
#                         cohort_sizes = rep(3, 13))

## ---- eval = FALSE-------------------------------------------------------
#  table(sims$recommended_dose) / length(sims$recommended_dose)

## ---- eval = FALSE-------------------------------------------------------
#  table(unlist(sims$doses_given)) / length(unlist(sims$doses_given))

## ---- eval = FALSE-------------------------------------------------------
#  table(unlist(sims$doses_given)) / length(sims$recommended_dose)

