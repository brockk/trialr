# R.version

# install.packages('rstan')
# install.packages('Rcpp')
# install.packages('rstantools')
# install.packages('gtools')
# install.packages('devtools')
# install.packages('curl')


# R.version
# library(Rcpp)
library(rstan)
# library(rstantools)
library(trialr)
# devtools::load_all()

# EffTox -------
help("efftox_analysis")
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

# Outcomes
efftox_parse_outcomes('1E 2E 3TTT') # OK
efftox_parse_outcomes('1E2E3TTT')  # Breaks! TODO
efftox_parse_outcomes(' 1ET 12EB ') # OK
efftox_parse_outcomes(' 1ET 12EB') # OK

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





# ThallHierarchicalBinary -------
dat <- thallhierarchicalbinary_parameters_demo()
samp = rstan::sampling(stanmodels$ThallHierarchicalBinary, data = dat)
plot(samp, pars = 'p')

plot(samp, pars = 'pg')
# Posterior Prob(Response)...
# In group 4
ggplot(data.frame(ProbResponse = extract(samp, 'p[4]')[[1]]),
       aes(x = ProbResponse)) + geom_density() + ggtitle('Prob(Response) in Sub-group 4')

ggplot(data.frame(ProbResponse = extract(samp, 'p[3]')[[1]]),
       aes(x = ProbResponse)) + geom_density() + ggtitle('Prob(Response) in Sub-group 3')


# BEBOP in PePS2 ------
set.seed(123)
dat <- peps2_get_data(num_patients = 60,
                      prob_eff = c(0.167, 0.192, 0.5, 0.091, 0.156, 0.439),
                      prob_tox = rep(0.1, 6),
                      #prob_tox = c(0.9, 0.15, 0.1, 0.9, 0.1, 0.1),
                      eff_tox_or = rep(1, 6))
samp = rstan::sampling(stanmodels$BebopInPeps2, data = dat)
colMeans(extract(samp, 'prob_eff')[[1]])
plot(samp, pars = 'prob_eff')

decision <- peps2_process(dat, samp)
decision$Accept
decision$ProbEff
decision$ProbAccEff
decision$ProbTox
decision$ProbAccTox

# Simulate
peps2_sc <- function() peps2_get_data(num_patients = 60,
                                      prob_eff = c(),
                                      prob_tox = 0.1,
                                      eff_tox_or = 1.0)
set.seed(123)
sims <- peps2_run_sims(num_sims = 10, sample_data_func = peps2_sc,
                       summarise_func = peps2_process)
apply(sapply(sims, function(x) x$Accept), 1, mean)


# CRM ----

# Example 1
library(dfcrm)
d <- c(1,1,1, 2,2,2, 2,2,2, 3,3,3)
tox <- c(0,0,0, 0,1,0, 0,0,0, 1,1,0)
target <- 0.25
skeleton <- c(0.1, 0.15, 0.25, 0.4, 0.65)
df_mod <- crm(prior = skeleton, target = target, tox = tox, level = d,
              model = 'empiric', method = 'bayes')
df_mod

crm_emp_dat <- list(beta_sd = sqrt(1.34),
                    num_doses = length(skeleton),
                    skeleton = skeleton,
                    num_patients = length(d),
                    tox = tox,
                    doses = d)
library(rstan)
crm_samp <- rstan::sampling(stanmodels$CrmEmpiricNormalPrior, data = crm_emp_dat, seed = 123)
# crm_emp_normal <- rstan::stan_model(file = 'exec_dev/CRM_Empirical_NormalPrior.stan')
# crm_samp <- rstan::sampling(object = crm_emp_normal, data = crm_emp_dat)
plot(crm_samp)
summary(crm_samp)$summary

# Example - p.21 Cheung (2011)
library(dfcrm)
target <- 0.25
skeleton <- c(0.05, 0.12, 0.25, 0.4, 0.55)
d <- c(3, 5, 5, 3, 4)
tox <- c(0, 0, 1, 0, 0)
crm(skeleton, target, tox, d, model = "logistic", intcpt = 3)

crm_logit1_dat <- list(a0 = 3,
                       beta_mean = 0,
                       beta_sd = sqrt(1.34),
                       num_doses = length(skeleton),
                       skeleton = skeleton,
                       num_patients = length(d),
                       tox = tox,
                       doses = d)

# crm_logit1_normal <- rstan::stan_model(file = 'exec/CRM_OneParamLogistic_NormalPrior.stan')
# crm_samp <- rstan::sampling(object = crm_logit1_normal, data = crm_logit1_dat)
crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticNormalPrior,
                            data = crm_logit1_dat, seed = 123)

summary(crm_samp)$summary


prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')

# Posterior mean
apply(prob_tox_samp, 2, mean)
# Posterior median
apply(prob_tox_samp, 2, median)
# Posterior probability that prob_tox exceeds target
apply(prob_tox_samp > target, 2, mean)
# Posterior probability that prob_tox exceeds target by 10%
apply(prob_tox_samp > target + 0.1, 2, mean)

# Dissonance over next dose
# Recommended dose using mean
which.min(abs(apply(prob_tox_samp, 2, mean) - target))
# Recommended dose using median
which.min(abs(apply(prob_tox_samp, 2, median) - target))
# Posterior MTD prob
table(apply(prob_tox_samp, 1, function(x) which.min(abs(x - target)))) / nrow(prob_tox_samp) # Note: apply by row, not col!!

pt_tall <- data.frame(ProbTox = as.vector(as.matrix(prob_tox_samp)),
                      Dose = as.factor(rep(1:ncol(prob_tox_samp), each = nrow(prob_tox_samp))))
ggplot(pt_tall, aes(x = ProbTox, group = Dose, col = Dose)) +
  geom_density()

# Levy example
# See Levy2006, Grieve2017
target <- 0.33
skeleton <- c(0.05, 0.10, 0.15, 0.33, 0.5)
d <- c(1, 1, 1)
tox <- c(0, 0, 0)
levy1_dat <- list(a0 = 3,
                  beta_shape = 1,
                  beta_inverse_scale = 1,
                  num_doses = length(skeleton),
                  skeleton = skeleton,
                  num_patients = length(d),
                  tox = tox,
                  doses = d)

crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = levy1_dat, seed = 123,
                            control = list(adapt_delta = 0.95))
summary(crm_samp)$summary
prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')

# Posterior mean
apply(prob_tox_samp, 2, mean)  # This matches bcrm but not Levy2006
# Posterior median
apply(prob_tox_samp, 2, median) # This is too high

levy1b_dat <- list(a0 = 3,
                   beta_mean = 0,
                   beta_sd = sqrt(1.34),
                   num_doses = length(skeleton),
                   skeleton = skeleton,
                   num_patients = length(d),
                   tox = tox,
                   doses = d)
crm_samp_norm <- rstan::sampling(stanmodels$CrmOneParamLogisticNormalPrior,
                            data = levy1b_dat, seed = 123)
beta_samp_1 <- exp(as.data.frame(crm_samp_norm, 'beta'))
summary(beta_samp_1)
beta_samp_2 <- as.data.frame(crm_samp, 'beta')
summary(beta_samp_2)
# They look kinda similar but the gamma prior yields much wider posterior, as
# might be expected when n=3

qgamma(p = c(0.05, 0.95), shape = 1, rate = 1) # 0.05 - 3.00
exp(qnorm(p = c(0.05, 0.95), mean = 0, sd = sqrt(1.34))) # 0.15 - 6.71
# The norm interval is much wider, and this is reflected in the posterior CI

# Posterior prob_tox
summary(crm_samp, 'prob_tox')$summary[, c('mean', 'sd', '50%')]
# The mean and SD match bcrm. But median does not
summary(crm_samp_norm, 'prob_tox')$summary[, c('mean', 'sd', '50%')]
# These all match bcrm


# Dissonance over next dose
# Recommended dose using mean
which.min(abs(apply(prob_tox_samp, 2, mean) - target)) # Dose 5
# Recommended dose using median
which.min(abs(apply(prob_tox_samp, 2, median) - target)) # Dose 5
# Posterior MTD prob
table(apply(prob_tox_samp, 1, function(x) which.min(abs(x - target)))) / nrow(prob_tox_samp) # Dose 5




sessionInfo()

# Vignettes  -------
# devtools::use_vignette("EffTox")
# devtools::use_vignette("BEBOP")
# devtools::use_vignette("HierarchicalBayesianResponse")
devtools::build_vignettes()
# Documentation
devtools::use_readme_rmd()
roxygen2::roxygenise()

# Help examples ----
help('trialr')
help("efftox_params")
help("efftox_parameters_demo")
help("efftox_analysis")
help("efftox_contour_plot")
help("efftox_utility")
help("efftox_process")

help("thallhierarchicalbinary_parameters_demo")

help("peps2_params")
help("peps2_get_data")
help("peps2_process")
help("peps2_run_sims")

devtools::build_vignettes()
# devtools::check()
devtools::check(manual = TRUE)
devtools::build()
devtools::release()

# RStan breaks? ----
# Time to do a rain dance:
# https://groups.google.com/forum/#!topic/stan-users/8e73htnTxro
# remove.packages("rstan")
# dir(system.file("libs", package = "rstan"))
# install.packages("rstan", dependencies = TRUE)
# Scrub the environ. Restart R session.
# Then Clean and Rebuild
# eight_schools <- stan_demo("eight_schools")
# eight_schools
