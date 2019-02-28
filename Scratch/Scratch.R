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
library(loo)


# ThallHierarchicalBinary -------
dat <- thallhierarchicalbinary_parameters_demo()
samp = rstan::sampling(stanmodels$ThallHierarchicalBinary, data = dat)
plot(samp, pars = 'p')

plot(samp, pars = 'pg')
# Posterior Prob(Response)...
# In group 4
ggplot(data.frame(ProbResponse = extract(samp, 'p[4]')[[1]]),
       aes(x = ProbResponse)) + geom_density() +
  ggtitle('Prob(Response) in Sub-group 4')

ggplot(data.frame(ProbResponse = extract(samp, 'p[3]')[[1]]),
       aes(x = ProbResponse)) + geom_density() +
  ggtitle('Prob(Response) in Sub-group 3')


# BEBOP in PePS2 ------
set.seed(123)
dat <- peps2_get_data(num_patients = 60,
                      prob_eff = c(0.167, 0.192, 0.5, 0.091, 0.156, 0.439),
                      prob_tox = rep(0.1, 6),
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

mod1 <- stan_crm(outcome_str = '1NNN 2NTN 2NNN 3TTN', skeleton = skeleton,
                 target = target, model = 'empiric', beta_sd = sqrt(1.34))
class(mod1)
mod1
head(as.data.frame(mod1, 'prob_tox'))
summary(mod1)
summary(mod1, 'prob_tox')
plot(mod1)
plot(mod1, pars = 'beta')

# Example - p.21 Cheung (2011)
library(dfcrm)
target <- 0.25
skeleton <- c(0.05, 0.12, 0.25, 0.4, 0.55)
d <- c(3, 5, 5, 3, 4)
tox <- c(0, 0, 1, 0, 0)
crm(skeleton, target, tox, d, model = "logistic", intcpt = 3)

mod2 <- stan_crm(outcome_str = '3N 5N 5T 3N 4N', skeleton = skeleton,
                 target = target, model = 'logistic',
                 a0 = 3, beta_mean = 0, beta_sd = sqrt(1.34))
mod2

library(magrittr)
library(ggplot2)
mod2 %>%
  gather_samples.crm_fit %>%
  ggplot(aes(x = DoseLevel, y = ProbTox, group = DoseLevel)) +
  geom_violin(fill = 'orange') + ylim(0, 1) +
  geom_hline(yintercept = target, col = 'red', linetype = 'dashed') +
  labs(title = 'Pr(DLT) after 18 patients in Levy, et al. (2006)')

mod2 %>%
  gather_samples.crm_fit %>%
  head %>%
  filter(as.name(".iteration") <= 1)


# TODO: Emulate Prob(too tox)
# E.g.
# apply(prob_tox_samp > target, 2, mean)
# apply(prob_tox_samp > target + 0.1, 2, mean)
# TODO: plots
# Post param summaries




# Control commands etc ----
sessionInfo()


# README & Vignettes  -------
# devtools::use_readme_rmd()
# devtools::use_vignette("trialr-overview")
# devtools::use_vignette("EffTox")
# devtools::use_vignette("BEBOP")
# devtools::use_vignette("HierarchicalBayesianResponse")
# devtools::use_vignette("CRM")
# devtools::use_vignette("CRM-visualisation")
# devtools::use_vignette("CRM-model-choice")
# devtools::use_build_ignore(c("Scratch", "exec_dev"))


# Documentation
# devtools::use_readme_rmd()





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

# Release etc
# devtools::check()
devtools::check(manual = TRUE)
devtools::build()
roxygen2::roxygenise()
devtools::build_vignettes()
pkgdown::build_site()
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

# Tests ----
devtools::use_testthat()
devtools::test()
