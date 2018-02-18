
# Levy, V., Zohar, S., Bardin, C., Vekhoff, A., Chaoui, D., Rio, B., … Marie, J. P. (2006).
# A phase I dose-finding and pharmacokinetic study of subcutaneous semisynthetic
# homoharringtonine (ssHHT) in patients with advanced acute myeloid leukaemia.
# British Journal of Cancer, 95(3), 253–259. https://doi.org/10.1038/sj.bjc.6603265

target <- 0.33
skeleton <- c(0.05, 0.10, 0.15, 0.33, 0.5)

# Cohort 1 ----
d <- c(1, 1, 1)
tox <- c(0, 0, 0)
levy_dat <- list(a0 = 3,
                  beta_shape = 1,
                  beta_inverse_scale = 1,
                  num_doses = length(skeleton),
                  skeleton = skeleton,
                  num_patients = length(d),
                  tox = tox,
                  doses = d)
crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = levy_dat, seed = 123,
                            control = list(adapt_delta = 0.95))
prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
apply(prob_tox_samp, 2, median) # 0.004260121 0.012230709 0.023349450 0.093317441 0.219839194
# Diverges from Table 1
which.min(abs(apply(prob_tox_samp, 2, median) - target)) # Dose 5

# Cohort 2 ----
d <- c(1, 1, 1, 3, 3, 3)
tox <- c(0, 0, 0, 0, 0, 1)
levy_dat <- list(a0 = 3,
                 beta_shape = 1,
                 beta_inverse_scale = 1,
                 num_doses = length(skeleton),
                 skeleton = skeleton,
                 num_patients = length(d),
                 tox = tox,
                 doses = d)
crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = levy_dat, seed = 123,
                            control = list(adapt_delta = 0.95))
prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
apply(prob_tox_samp, 2, median) # 0.08488401  0.15423275  0.21699597  0.41224564  0.57101296
# Getting closer to Table 1
which.min(abs(apply(prob_tox_samp, 2, median) - target)) # Dose 4

# Cohort 3 ----
d <- c(1, 1, 1, 3, 3, 3, 4, 4, 4)
tox <- c(0, 0, 0, 0, 0, 1, 0, 0, 1)
levy_dat <- list(a0 = 3,
                 beta_shape = 1,
                 beta_inverse_scale = 1,
                 num_doses = length(skeleton),
                 skeleton = skeleton,
                 num_patients = length(d),
                 tox = tox,
                 doses = d)
crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = levy_dat, seed = 123,
                            control = list(adapt_delta = 0.95))
prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
apply(prob_tox_samp, 2, median) # 0.07494966  0.13942849  0.19924889  0.39196263  0.55421598
# Close to Table 1
which.min(abs(apply(prob_tox_samp, 2, median) - target)) # Dose 4

# Cohort 4 ----
d = c(1, 1, 1, 3, 3, 3, 4, 4, 4, 4, 4, 4)
tox  = c(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0)
levy_dat <- list(a0 = 3,
                 beta_shape = 1,
                 beta_inverse_scale = 1,
                 num_doses = length(skeleton),
                 skeleton = skeleton,
                 num_patients = length(d),
                 tox = tox,
                 doses = d)
crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = levy_dat, seed = 123,
                            control = list(adapt_delta = 0.95))
prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
apply(prob_tox_samp, 2, median) #  0.03480635  0.07393290  0.11550117  0.28008492  0.45244122
# Close to Table 1
which.min(abs(apply(prob_tox_samp, 2, median) - target)) # Dose 4

# Cohort 5 ----
d = c(1, 1, 1, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4)
tox  = c(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0)
levy_dat <- list(a0 = 3,
                 beta_shape = 1,
                 beta_inverse_scale = 1,
                 num_doses = length(skeleton),
                 skeleton = skeleton,
                 num_patients = length(d),
                 tox = tox,
                 doses = d)
crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = levy_dat, seed = 123,
                            control = list(adapt_delta = 0.95))
prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
apply(prob_tox_samp, 2, median) # 0.03836539  0.08021418  0.12398132  0.29296965  0.46510317
# Close to Table 1
which.min(abs(apply(prob_tox_samp, 2, median) - target)) # Dose 4

# Cohort 6 ----
d = c(1, 1, 1, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4)
tox  = c(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1)
levy_dat <- list(a0 = 3,
                 beta_shape = 1,
                 beta_inverse_scale = 1,
                 num_doses = length(skeleton),
                 skeleton = skeleton,
                 num_patients = length(d),
                 tox = tox,
                 doses = d)
crm_samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = levy_dat, seed = 123,
                            control = list(adapt_delta = 0.95))
prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
apply(prob_tox_samp, 2, median) # 0.0645770   0.1234640   0.1796798   0.3684325   0.5341758
# Close to Table 1
which.min(abs(apply(prob_tox_samp, 2, median) - target)) # Dose 4

