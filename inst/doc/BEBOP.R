## ---- echo=FALSE---------------------------------------------------------
knitr::kable(data.frame(i = 1:6, 
                        Pretreated = c(F,F,F,T,T,T),
                        PDL1 = c('Low', 'Medium', 'High', 'Low', 'Medium', 'High'),
                        x1 = rep(c(0,1), each=3),
                        x2 = c(1,0,0,1,0,0),
                        x3 = c(0,1,0,0,1,0)
                        )
             )

## ----message=FALSE-------------------------------------------------------
library(trialr)
peps2_sc <- function() peps2_get_data(num_patients = 60,
                                      prob_eff = c(0.167, 0.192, 0.5, 0.091, 0.156, 0.439),
                                      prob_tox = rep(0.1, 6),
                                      eff_tox_or = rep(1, 6))

set.seed(123)
dat <- peps2_sc()

## ------------------------------------------------------------------------
c(dat$alpha_mean, dat$alpha_sd)

## ------------------------------------------------------------------------
knitr::kable(
  head(with(dat, data.frame(eff, tox, x1, x2, x3)), 10)
)

## ---- results = "hide"---------------------------------------------------
samp <- rstan::sampling(stanmodels$BebopInPeps2, data = dat)

## ---- fig.width = 6, fig.height = 6, fig.cap = "Posterior Prob(Efficacy) in the six PePS2 cohorts", eval=FALSE----
#  rstan::plot(samp, pars = 'prob_eff')

## ---- eval=T-------------------------------------------------------------
decision <- peps2_process(dat, samp)
knitr::kable(
  with(decision, data.frame(ProbEff, ProbAccEff, ProbTox, ProbAccTox, Accept)), 
  digits = 3
)

## ---- eval=FALSE---------------------------------------------------------
#  set.seed(123)
#  sims <- peps2_run_sims(num_sims = 10, sample_data_func = peps2_sc,
#                         summarise_func = peps2_process)

## ---- eval=FALSE---------------------------------------------------------
#  apply(sapply(sims, function(x) x$Accept), 1, mean)

