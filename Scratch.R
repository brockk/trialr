
# library(trialr)

# Eight schools
dat = eight_schools_get_data()
samp = rstan::sampling(stanmodels$EightSchools, data = dat)
plot(samp)


# EffTox
eff0 = 0.5
tox1 = 0.65
p = 2
p_e = 0.10
p_t = 0.10
max_patients = 39

dat <- list(
  num_doses = 5,
  real_doses = c(1, 2, 4, 6.6, 10),
  efficacy_hurdle = 0.5,
  toxicity_hurdle = 0.3,
  p = p,
  eff0 = eff0,
  tox1 = tox1,

  alpha_mean = -7.9593,
  alpha_sd = 3.5487,
  beta_mean = 1.5482,
  beta_sd = 3.5018,
  gamma_mean = 0.7367,
  gamma_sd = 2.5423,
  zeta_mean = 3.4181,
  zeta_sd = 2.4406,
  eta_mean = 0,
  eta_sd = 0.2,
  psi_mean = 0,
  psi_sd = 1,

  num_patients = 3,
  eff = c(0,1,1),
  tox = c(0,1,0),
  doses = c(1,1,1)
)

samp = rstan::sampling(stanmodels$EffTox, data = dat)
plot(samp)



# ThallHierarchicalBinary
dat = list(
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


# Vignettes
devtools::use_vignette("EffTox")


