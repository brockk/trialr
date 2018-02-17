## ----message=FALSE-------------------------------------------------------
library(trialr)
dat <- efftox_parameters_demo()

## ------------------------------------------------------------------------
p <- efftox_solve_p(eff0 = 0.5, tox1 = 0.65, eff_star = 0.7, tox_star = 0.25)
p

## ------------------------------------------------------------------------
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
  
  eff = c(),
  tox = c(),
  doses = c(),
  num_patients = 0
)

## ------------------------------------------------------------------------
dat$doses = c(1, 1, 1, 2, 2, 2)
dat$tox   = c(0, 0, 0, 0, 1, 0)
dat$eff   = c(0, 1, 0, 0, 1, 1)
dat$num_patients = 6

## ---- results = "hide"---------------------------------------------------
samp <- rstan::sampling(stanmodels$EffTox, data = dat, seed = 123)

## ------------------------------------------------------------------------
x <- efftox_process(dat, samp)

## ------------------------------------------------------------------------
x$recommended_dose

## ----eval=FALSE----------------------------------------------------------
#  knitr::kable(efftox_analysis_to_df(x), digits = 3)

## ----echo=FALSE----------------------------------------------------------
knitr::kable(efftox_analysis_to_df(x), digits = 3, caption = 'ProbEff and ProbTox are the probabilities of efficacy and toxicity at each dose. ProbAccEff is the probability that the efficacy rate exceeds the desired threshold and ProbAccTox the probability that the toxicity rate is less than the threshold. All probabilities are posterior means.')

## ---- fig.width = 6, fig.height = 6, fig.cap = "Utility contours after observing outcomes 1NEN 2NBE."----
efftox_contour_plot(dat, prob_eff = x$prob_eff, prob_tox = x$prob_tox)
title('EffTox utility contours')

## ---- fig.width = 6, fig.height = 6, fig.cap = "Utility densities after observing outcomes 1NEN 2NBE."----
efftox_utility_density_plot(samp, doses = 1:3) +
  ggplot2::ggtitle("EffTox dose utility densities")

## ------------------------------------------------------------------------
knitr::kable(efftox_superiority(samp), digits = 2, row.names = TRUE)

## ------------------------------------------------------------------------
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

