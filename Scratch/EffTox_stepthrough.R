
# Using example trial in:
# Thall, P. F., Herrick, R. C., Nguyen, H. Q., Venier, J. J., & Norris, J. C. (2014).
# Effective sample size for computing prior hyperparameters in Bayesian phase I-II dose-finding.
# Clinical Trials, 11(6), 657–666. http://doi.org/10.1177/1740774514547397

eff0 = 0.5
tox1 = 0.65
p = efftox_solve_p(eff0 = eff0, tox1 = tox1, eff_star = 0.7, tox_star = 0.25)
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

  num_patients = 0,
  eff = c(),
  tox = c(),
  doses = c()
)


# Scenario 1
true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90)
true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40)

set.seed(seed)
dose = 1  # Starting dose
fit = NULL

# Cohort 1
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 0,0,0
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,0,0
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 2

# Cohort 2
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,0,0
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 1,0,1
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 3

# Cohort 3
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,0
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,0,0
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 4

# Cohort 4
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,0,0
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 5
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,0
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,0,1
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 6
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,1,1
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 7
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,0,1
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 8
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,1,0
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 9
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,1,0
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 10
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,1,0
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 11
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 0,1,0
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,1,0
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 12
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,0,1
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 13
dat$num_patients < max_patients
cohort_size = 3
prob_eff = true_eff[dose]
prob_tox = true_tox[dose]
# Simulate new efficacy events
new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)  # 1,1,1
new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)  # 0,0,1
# And append to trial data
dat$eff = c(dat$eff, new_eff)
dat$tox = c(dat$tox, new_tox)
# Also reflect doses delivered
dat$doses <- c(dat$doses, rep(dose, cohort_size))
dat$num_patients = dat$num_patients + cohort_size
fit <- efftox_fit(dat, fit = fit, seed = seed)
l <- efftox_process(dat, fit, p_e = p_e, p_t = p_e)
# Select a dose?
sum(l$admissible) > 0
# Select dose
dose = which.max(ifelse(l$admissible, l$utility, NA))  # 5

# Cohort 14?
dat$num_patients < max_patients
# No
