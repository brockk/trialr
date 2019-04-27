
#' @title Process RStan samples from an EffTox model
#'
#' @description Internal function to process rstan samples from an EffTox model
#' to make inferences about dose-acceptability, dose-utility and which dose
#' should be recommended next.
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param fit An instance of \code{rstan::stanmodel}, derived by fitting the
#' trialr EffTox model.
#' @return An instance of \code{\link{efftox_fit}}.
#'
efftox_process <- function(dat, fit) {
  # Posterior mean estimates
  prob_eff <- colMeans(rstan::extract(fit, 'prob_eff')[[1]])
  prob_acc_eff <- colMeans(rstan::extract(fit, 'prob_acc_eff')[[1]])
  prob_tox <- colMeans(rstan::extract(fit, 'prob_tox')[[1]])
  prob_acc_tox <- colMeans(rstan::extract(fit, 'prob_acc_tox')[[1]])
  post_utility <- colMeans(rstan::extract(fit, 'utility')[[1]])
  # Derived quantities
  utility = efftox_utility(dat$p, dat$eff0, dat$tox1,
                           prob_eff, prob_tox)
  # Dose admissibility
  dose_indices <- 1:dat$num_doses
  lowest <- min(dat$doses)
  highest <- max(dat$doses)
  in_range <- sapply(dose_indices,
                     function(x) (x >= lowest - 1) & (x <= highest + 1))
  acceptable <- (prob_acc_eff > dat$p_e) & (prob_acc_tox > dat$p_t) & in_range
  if(sum(acceptable) > 0) {
    recommended_dose <- which.max(ifelse(acceptable, utility, NA))  # 2
  } else {
    recommended_dose <- NA
  }

  x <- efftox_fit(dose_indices, recommended_dose, prob_eff, prob_tox,
                  prob_acc_eff, prob_acc_tox, utility, post_utility,
                  acceptable, dat, fit)
  return(x)
}
