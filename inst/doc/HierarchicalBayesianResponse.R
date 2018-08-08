## ----message=FALSE-------------------------------------------------------
library(trialr)
help(thallhierarchicalbinary_parameters)

## ------------------------------------------------------------------------
dat <- thallhierarchicalbinary_parameters_demo()
dat

## ---- results = "hide"---------------------------------------------------
samp <- rstan::sampling(stanmodels$ThallHierarchicalBinary, data = dat, 
                        seed = 123)

## ------------------------------------------------------------------------
knitr::kable(rstan::summary(samp, par = 'pg')$summary, digits = 3)

## ---- fig.width = 6, fig.height = 6, fig.cap = "Prob(Response | D) in subgroup 3"----
library(ggplot2)
library(rstan)
library(dplyr)

as.data.frame(samp, 'p[3]') %>% 
  mutate(ProbResponse = `p[3]`) %>% 
  ggplot(aes(x = ProbResponse)) + 
  geom_density() + 
  ggtitle('Prob(Response | D) in Sub-group 3')

## ---- fig.width = 6, fig.height = 6, fig.cap = "Prob(Response | D) in subgroup 4"----
as.data.frame(samp, 'p[4]') %>% 
  mutate(ProbResponse = `p[4]`) %>% 
  ggplot(aes(x = ProbResponse)) + 
  geom_density() + 
  ggtitle('Prob(Response | D) in Sub-group 4')

