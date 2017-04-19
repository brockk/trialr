
# R.version

# install.packages('rstan')
# install.packages('Rcpp')
# install.packages('rstantools')
# install.packages('gtools')
# install.packages('devtools')

# library(rstan)
# library(rstantools)
# library(Rcpp)
# library(trialr)
devtools::load_all()

# EffTox -------
# Get parameters for EffTox model
dat <- efftox_parameters_demo()
# Add outcomes for 3 patients: Neither event in patient 1; both efficacy and
# toxicity in patient 2; and just efficacy in patient 3
dat$num_patients <- 3
dat$eff <- c(0, 1, 1)
dat$tox <- c(0, 0, 1)
dat$doses <- c(1, 2, 3)

# Invoke RStan posterior sampling on model and data
samp <- rstan::sampling(stanmodels$EffTox, data = dat, seed = 123)
decision <- efftox_process(dat, samp, p_e = 0.10, p_t = 0.10)
decision
decision$recommended_dose  # 3
round(decision$utility, 2)
#  -0.64  0.04  0.25 -0.05 -0.20 on Win
#  -0.64  0.05  0.26 -0.04 -0.19 on Mac
plot(samp, par = 'utility') + ggtitle('Utility of doses after outcomes: 1NBE')

efftox_contour_plot(dat, prob_eff = decision$prob_eff, prob_tox = decision$prob_tox)
title('EffTox utility contours')

efftox_utility_density_plot(samp, doses = 1:3) +
  ggtitle("EffTox dose utility densities")

(sup_mat <- efftox_superiority(samp))
# Probability that utility of the dose in column i exceeds the utility of the
# dose in row j
# We propose that the least of these for each dose be used to infer the
# probability that that dose is superior to all others
dose_sup <- apply(sup_mat, 2, min, na.rm = TRUE)
round(dose_sup, 2)  # 0.05 0.37 0.63 0.23 0.20
# Dose 3 appears to be superior to all others, based on the limited data

# Simulate
dat <- efftox_parameters_demo()
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
devtools::build_vignettes()



