## ----message=FALSE-------------------------------------------------------
library(trialr)

## ---- results = "hide"---------------------------------------------------
fit <- stan_hierarchical_response_thall(
  group_responses = c(0, 0, 1, 3, 5, 0, 1, 2, 0, 0), 
  group_sizes = c(0, 2 ,1, 7, 5, 0, 2, 3, 1, 0), 
  mu_mean = -1.3863,
  mu_sd = sqrt(1 / 0.1),
  tau_alpha = 2,
  tau_beta = 20)

## ------------------------------------------------------------------------
fit

## ------------------------------------------------------------------------
knitr::kable(rstan::summary(fit, par = 'prob_response')$summary, digits = 3)

## ------------------------------------------------------------------------
colMeans(as.data.frame(fit, pars = 'prob_response') > 0.3)

## ---- message = FALSE, fig.width = 7, fig.height = 7, fig.cap = "Prob(Response | D) in subgroup 3"----
library(ggplot2)
library(rstan)
library(dplyr)

plot(fit, pars = 'prob_response') + 
  geom_vline(xintercept = 0.3, col = 'orange', linetype = 'dashed') +
  labs(title = 'Partially-pooled analysis of response rate in 10 sarcoma subtypes')

