
# library(trialr)

# EffTox -------
eff0 = 0.5
tox1 = 0.65
eff_star = 0.7
tox_star = 0.25
p = efftox_solve_p(eff0, tox1, eff_star, tox_star)
p_e = 0.10
p_t = 0.10

dat <- list(
  num_doses = 5,
  real_doses = c(1, 2, 4, 6.6, 10),
  efficacy_hurdle = 0.5,
  toxicity_hurdle = 0.3,
  p = p,
  eff0 = eff0,
  tox1 = tox1,

  alpha_mean = -7.9593, alpha_sd = 3.5487,
  beta_mean = 1.5482, beta_sd = 3.5018,
  gamma_mean = 0.7367, gamma_sd = 2.5423,
  zeta_mean = 3.4181, zeta_sd = 2.4406,
  eta_mean = 0, eta_sd = 0.2,
  psi_mean = 0, psi_sd = 1,

  num_patients = 3,
  eff = c(0,1,1),
  tox = c(0,1,0),
  doses = c(1,1,1)
)

samp = rstan::sampling(stanmodels$EffTox, data = dat)
class(samp)
plot(samp)

# Simulate
dat$eff = c()
dat$tox = c()
dat$doses = c()
dat$num_patients = 0
# TODO: provide function to strip patient outcomes? Doing just one of the above steps would be problematic.
set.seed(123)
sims = efftox_simulate(dat, num_sims = 2, first_dose = 1, p_e, p_t,
                       true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
                       true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
                       cohort_sizes = rep(3, 13),
                       chains = 2)
# Running 10 takes 2min30 = 150s => 150 / (4*13*10) = 150 / 520 = 0.3s per cycle to sample...1000?
table(sims$recommended_dose) / length(sims$recommended_dose)
table(unlist(sims$doses_given)) / length(unlist(sims$doses_given))
table(unlist(sims$doses_given)) / length(sims$recommended_dose)




# ThallHierarchicalBinary -------
dat <- list(
  m = 10,
  x = c(0, 0, 1, 3, 5, 0, 1, 2, 0, 0),  # Num responses, by cohort
  n = c(0, 2 ,1, 7, 5, 0, 2, 3, 1, 0),  # Num patients, by cohort
  target_resp = 0.3,
  mu_mean = -1.3863,
  mu_sd = sqrt(1 / 0.1),
  tau_alpha = 2,  # Thall et al. give a prior for the precision (tau) of rho,
  # as was the custom in BUGS.
  tau_beta = 20   # I use stdev (sigma) but retain his prior.
)
samp = rstan::sampling(stanmodels$ThallHierachicalBinary, data = dat)
plot(samp)


# BEBOP in PePS2 ------
set.seed(123)
dat <- peps2_get_data(num_patients = 60,
                      prob_eff = rep(0.3, 6),
                      prob_tox = 0.1,
                      eff_tox_or = 1.0)
samp = rstan::sampling(stanmodels$BebopInPeps2, data = dat)
plot(samp)
d <- peps2_trial_decision(dat, samp)
d$Decision
dat$cohort_eff

# Simulate
peps2_sc_1 <- function() peps2_get_data(num_patients = 60,
                                        prob_eff = rep(0.3, 6),
                                        prob_tox = 0.1,
                                        eff_tox_or = 1.0)
set.seed(123)
sims <- peps2_run_sims(num_sims = 10, sample_data_func = peps2_sc_1,
                       summarise_func = peps2_trial_decision)
apply(sapply(sims, function(x) x$Decision), 1, mean)


# Vignettes  -------
devtools::use_vignette("EffTox")
devtools::use_vignette("BEBOP")
devtools::use_vignette("HierarchicalBayesianResponse")



