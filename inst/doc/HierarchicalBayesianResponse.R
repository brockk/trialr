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
ggplot(data.frame(ProbResponse = extract(samp, 'p[3]')[[1]]),
       aes(x = ProbResponse)) + geom_density() + ggtitle('Prob(Response | D) in Sub-group 3')

## ---- fig.width = 6, fig.height = 6, fig.cap = "Prob(Response | D) in subgroup 4"----
ggplot(data.frame(ProbResponse = extract(samp, 'p[4]')[[1]]),
       aes(x = ProbResponse)) + geom_density() + ggtitle('Prob(Response | D) in Sub-group 4')

