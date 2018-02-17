
d <- c(1,1,1, 2,2,2, 2,2,2, 3,3,3)
tox <- c(0,0,0, 0,1,0, 0,0,0, 1,1,0)
target <- 0.25
skeleton <- c(0.1, 0.15, 0.25, 0.4, 0.65)

# Empiric model ----
crm_emp_normal <- rstan::stan_model(file = 'exec_dev/CRM_Empirical_NormalPrior.stan')
crm_emp_dat <- list(beta_sd = sqrt(1.34),
                    num_doses = length(skeleton),
                    skeleton = skeleton,
                    num_patients = length(d),
                    tox = tox,
                    doses = d)
crm_samp <- rstan::sampling(object = crm_emp_normal, data = crm_emp_dat, seed = 1)
plot(crm_samp, par = 'prob_tox')
summary(crm_samp, par = c('prob_tox', 'beta'))$summary

# Compare to dfcrm
library(dfcrm)
df_mod <- crm(prior = skeleton, target = target, tox = tox, level = d,
              model = 'empiric', method = 'bayes')
df_mod

library(loo)
log_lik_1 <- extract_log_lik(crm_samp)
loo_1 <- loo(log_lik_1)
print(loo_1)

# Logistic model ----
crm_logit1_normal <- rstan::stan_model(file = 'exec_dev/CRM_OneParamLogistic_NormalPrior.stan')
crm_logit1_dat <- list(a0 = 3,
                       beta_mean = 0,
                       beta_sd = sqrt(1.34),
                       num_doses = length(skeleton),
                       skeleton = skeleton,
                       num_patients = length(d),
                       tox = tox,
                       doses = d)
crm_samp2 <- rstan::sampling(object = crm_logit1_normal, data = crm_logit1_dat)
summary(crm_samp2)$summary

data.frame(emp = summary(crm_samp, 'prob_tox')$summary[, 'mean'],
     logit1 = summary(crm_samp2, 'prob_tox')$summary[, 'mean'])


log_lik_2 <- extract_log_lik(crm_samp2)
loo_2 <- loo(log_lik_2)

print(loo_1)
print(loo_2)
compare(loo_1, loo_2)
# Small difference.

# See:
# Vehtari, A., Gelman, A., and Gabry, J. (2016).
# Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.

library(tidyverse)

prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
prob_tox_samp_tall <- prob_tox_samp %>%
  gather(Label, ProbTox) %>%
  mutate(
    DoseLevel = rep(1:ncol(prob_tox_samp), each = nrow(prob_tox_samp)),
    Draw = rep(1:nrow(prob_tox_samp), times = ncol(prob_tox_samp))
  )
prob_tox_samp_tall %>% head()

ggplot(prob_tox_samp_tall, aes(x = DoseLevel, y = ProbTox, group = DoseLevel)) +
  geom_boxplot()
ggplot(prob_tox_samp_tall, aes(x = DoseLevel, y = ProbTox, group = DoseLevel)) +
  geom_violin()
ggplot(prob_tox_samp_tall, aes(x = DoseLevel, y = ProbTox, group = Draw)) +
  geom_line(alpha = 0.1, col = 'orange')
