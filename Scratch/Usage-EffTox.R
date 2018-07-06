
library(trialr)
library(loo)

# Usage: EffTox -------
help("stan_efftox")

mod1 <- stan_efftox(outcome_str = '1N 2E 3B',
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
# OR
mod1 <- stan_efftox_demo('1N 2E 3B', seed = 123)
mod1

mod1$recommended_dose
mod1$utility

plot(mod1, pars = 'prob_eff') +
  ggtitle('Utility of doses after outcomes: 1NBE')

# Contour plot
efftox_contour_plot(mod1$dat, prob_eff = mod1$prob_eff, prob_tox = mod1$prob_tox)
title('EffTox utility contours')
# Or
efftox_contour_plot(mod1$dat, prob_eff = mod1$prob_eff, prob_tox = mod1$prob_tox,
                    use_ggplot = TRUE) + ggtitle('EffTox utility contours')
# Change to use mod1

# Utility density plot
library(tidyr)
library(dplyr)
library(ggplot2)

mod1 %>%
  as.data.frame(pars = 'utility') %>%
  tidyr::gather('Var', 'Utility') %>%
  mutate(Dose = factor(stringr::str_extract(Var, '\\d+'))) %>%
  filter(Dose %in% 1:3) %>%
  ggplot(aes(x = Utility, group = Dose)) +
  geom_density(aes(fill = Dose))

library(ggridges)
mod1 %>%
  as.data.frame(pars = 'utility') %>%
  tidyr::gather('Var', 'Utility') %>%
  mutate(Dose = factor(stringr::str_extract(Var, '\\d+'))) %>%
  filter(Dose %in% 1:3) %>%
  ggplot(aes(x = Utility, y = Dose, fill = Dose)) +
  geom_density_ridges()

# Superiority
sup_mat <- efftox_superiority(mod1$fit)

# DTPs
# Star
mod1 <- stan_efftox_demo('1N 2E 3B', seed = 123)
dtps1 <- efftox_dtps(mod1$dat, cohort_sizes = c(1, 1),
                     next_dose = mod1$recommended_dose)
dtps1

# dat <- efftox_params(real_doses = c(1, 2, 4, 6),
#                      efficacy_hurdle = 0.3, toxicity_hurdle = 0.3,
#                      p_e = 0.1, p_t = 0.1, eff0 = 0.2, tox1 = 0.8,
#                      eff_star = 0.4, tox_star = 0.2,
#                      alpha_mean = -4.0, alpha_sd = 2.8,
#                      beta_mean = 2.8, beta_sd = 2.8,
#                      gamma_mean = -1.9, gamma_sd = 2.5,
#                      zeta_mean = 3.2, zeta_sd = 2.5,
#                      eta_mean = 0, eta_sd = 0.2,
#                      psi_mean = 0, psi_sd = 1)
# dtps1 <- efftox_dtps(dat, cohort_sizes = c(3), next_dose = 1)
# dtps1

# Simulate
dat <- efftox_parameters_demo()
dat$doses = c()
dat$eff = c()
dat$tox = c()
# I think NULLs do not show up
dat$num_patients = 0
set.seed(123)
sims = efftox_simulate(dat, num_sims = 10, first_dose = 1,
                       true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
                       true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
                       cohort_sizes = rep(3, 13),
                       chains = 2)
sims

# Deprecate?
efftox_utility_density_plot(mod1$fit, doses = 1:3) +
  ggtitle("EffTox dose utility densities")


# The old way ----
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
decision <- efftox_process(dat, samp)
decision

decision$recommended_dose  # 3
round(decision$utility, 2)
#  -0.64  0.04  0.25 -0.05 -0.20 on Win
#  -0.64  0.05  0.26 -0.04 -0.19 on Mac
plot(samp, par = 'utility') + ggtitle('Utility of doses after outcomes: 1NBE')

# Contour plot
efftox_contour_plot(dat, prob_eff = decision$prob_eff, prob_tox = decision$prob_tox)
title('EffTox utility contours')
# Or
efftox_contour_plot(dat, prob_eff = decision$prob_eff, prob_tox = decision$prob_tox,
                    use_ggplot = TRUE) + ggtitle('EffTox utility contours')

# Utility density plot
efftox_utility_density_plot(samp, doses = 1:3) +
  ggtitle("EffTox dose utility densities")

# Superiority
(sup_mat <- efftox_superiority(samp))
# Probability that utility of the dose in column i exceeds the utility of the
# dose in row j
# We propose that the least of these for each dose be used to infer the
# probability that that dose is superior to all others
dose_sup <- apply(sup_mat, 2, min, na.rm = TRUE)
round(dose_sup, 2)  # 0.05 0.35 0.65 0.23 0.20
# Dose 3 appears to be superior to all others, based on the limited data

# DTPs

# cohort_sizes = c(2, 3)
# next_dose = 2

# dat = efftox_parameters_demo()
# dat$doses = c(1,2,2)
# dat$eff = c(0,1,1)
# dat$tox = c(0,0,1)
# dat$num_patients = 3

# dat = efftox_parameters_demo(); cohort_sizes = c(1,1,1);
# next_dose = 1

dat <- efftox_parameters_demo()
dtps1 <- efftox_dtps(dat = dat, cohort_sizes = c(3), next_dose = 1)
dtps1

dat <- efftox_parameters_demo()
dat$doses = array(c(1,1,1))
dat$eff = array(c(0,0,0))
dat$tox = array(c(1,1,1))
dat$num_patients = 3
samp <- rstan::sampling(stanmodels$EffTox, data = dat, seed = 123)
decision <- efftox_process(dat, samp)
decision
# This is different to MD Anderson. TODO

dat <- efftox_parameters_demo()
dat$doses = array(1)
dat$eff = array(1)
dat$tox = array(1)
dat$num_patients = 1
dtps2 <- efftox_dtps(dat = dat, cohort_sizes = c(1, 1, 1),
                     next_dose = 1)
dtps2


# Simulate
dat <- efftox_parameters_demo()
set.seed(123)
# sims = efftox_simulate(dat, num_sims = 2, first_dose = 1, p_e, p_t,
#                        true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
#                        true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
#                        cohort_sizes = rep(3, 13),
#                        chains = 2)
# Running 10 takes 2min30 = 150s => 150 / (4*13*10) = 150 / 520 = 0.3s per cycle to sample...1000?
# table(sims$recommended_dose) / length(sims$recommended_dose)
# table(unlist(sims$doses_given)) / length(unlist(sims$doses_given))
# table(unlist(sims$doses_given)) / length(sims$recommended_dose)




