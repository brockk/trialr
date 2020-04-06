
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

  dose_indices <- 1:dat$num_doses

  # Posterior estimates
  prob_eff_samp <- rstan::extract(fit, 'prob_eff')[[1]]
  prob_eff <- colMeans(prob_eff_samp)
  median_prob_eff <- apply(prob_eff_samp, 2, stats::median)
  prob_acc_eff <- colMeans(prob_eff_samp > dat$efficacy_hurdle)
  prob_tox_samp <- rstan::extract(fit, 'prob_tox')[[1]]
  prob_tox <- colMeans(prob_tox_samp)
  median_prob_tox <- apply(prob_tox_samp, 2, stats::median)
  prob_acc_tox <- colMeans(prob_tox_samp < dat$toxicity_hurdle)
  post_utility_samp <- rstan::extract(fit, 'utility')[[1]]
  post_utility <- colMeans(post_utility_samp)
  obd_candidate <- apply(post_utility_samp, 1, function(x) which.max(x))
  prob_obd <- sapply(dose_indices, function(x) mean(obd_candidate == x))

  # Derived quantities
  # Utility estimate from plugging in prob_eff and prob_tox, a la Thall et al
  utility = efftox_utility(dat$p, dat$eff0, dat$tox1, prob_eff, prob_tox)
  # post_utility in contrast is full Bayesian posterior mean

  # Dose admissibility and recommended dose
  if(dat$num_patients > 0) {
    lowest <- min(dat$doses)
    highest <- max(dat$doses)
    in_range <- sapply(dose_indices,
                       function(x) (x >= lowest - 1) & (x <= highest + 1))
    acceptable <- (prob_acc_eff > dat$p_e) & (prob_acc_tox > dat$p_t) & in_range
    if(sum(acceptable) > 0) {
      recommended_dose <- which.max(ifelse(acceptable, utility, NA))
    } else {
      recommended_dose <- NA
    }
  } else {
    acceptable <- rep(NA, dat$num_doses)
    recommended_dose <- NA
  }

  x <- efftox_fit(dose_indices = dose_indices,
                  num_patients = dat$num_patients,
                  doses = dat$doses,
                  tox = dat$tox,
                  eff = dat$eff,
                  prob_tox = prob_tox,
                  prob_eff = prob_eff,
                  median_prob_tox = median_prob_tox,
                  median_prob_eff = median_prob_eff,
                  prob_acc_tox = prob_acc_tox,
                  prob_acc_eff = prob_acc_eff,
                  utility = utility,
                  post_utility = post_utility,
                  prob_obd = prob_obd,
                  acceptable = acceptable,
                  recommended_dose = recommended_dose,
                  dat = dat,
                  fit = fit)
  return(x)
}
