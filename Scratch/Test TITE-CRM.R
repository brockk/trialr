
# remove.packages('trialr')
library(trialr)

# CRM
mod1 <- stan_crm('1N 2N 3T', skeleton = c(0.1, 0.2, 0.35, 0.6),
                 target = 0.2, model = 'empiric', beta_sd = sqrt(1.34),
                 seed = 123)
mod1
mod1$recommended_dose

mod2 <- stan_crm('1NNN 2NNN 3TTT', skeleton = c(0.1, 0.2, 0.35, 0.6),
                 target = 0.2, model = 'logistic', a0 = 3, beta_mean = 0,
                 beta_sd = sqrt(1.34), seed = 123)
mod2

# TITE-CRM
# Logistic, Cheung p.124
mod3a <-stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
                doses_given = c(3, 3, 3, 3),
                tox = c(0, 0, 0, 0),
                weights = c(73, 66, 35, 28) / 126,
                model = 'logistic', a0 = 3, beta_mean = 0, beta_sd = sqrt(1.34),
                seed = 123)
mod3a

mod3b <-stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
                 doses_given = c(3, 3, 3, 3),
                 tox = c(0, 0, 0, 0),
                 model = 'logistic', a0 = 3, beta_mean = 0, beta_sd = sqrt(1.34),
                 seed = 123)
mod3b

library(dplyr)
mod3a$fit %>% as.data.frame() %>% summarise(mean(beta))
mod3b$fit %>% as.data.frame() %>% summarise(mean(beta))
mod3a$fit %>% as.data.frame('prob_tox') %>% apply(2, median)
mod3b$fit %>% as.data.frame('prob_tox') %>% apply(2, median)



# Empiric
mod4a <-stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
                 doses_given = c(3, 3, 3, 3),
                 tox = c(0, 0, 0, 0),
                 weights = c(73, 66, 35, 28) / 126,
                 model = 'empiric', beta_sd = sqrt(1.34),
                 seed = 123)
mod4a
mod4a$fit

mod4b <-stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
                 doses_given = c(3, 3, 3, 3),
                 tox = c(0, 0, 0, 0),
                 model = 'empiric', beta_sd = sqrt(1.34),
                 seed = 123)
mod4b

mod4c <-stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
                 doses_given = c(3, 3, 3, 3),
                 tox = c(0, 0, 0, 0),
                 weights = c(0.01, 0.01, 0.01, 0.01),
                 model = 'empiric', beta_sd = sqrt(1.34),
                 seed = 123)
mod4c

mod4a$fit %>% as.data.frame() %>% head
mod4b$fit %>% as.data.frame() %>% head
mod4c$fit %>% as.data.frame() %>% head



cat(rstan::get_stancode(mod4b$fit))



mod4c$dat
