#' @title Process RStan samples from a CRM model.
#'
#' @description Internal function to process rstan samples from a CRM model to
#' make inferences about dose-toxicity and which dose should be recommended next.
#' Typically, this function is not required to be called explicitly by the user
#' because \code{\link{stan_crm}} will call it implicitly.
#'
#' @param dat An instance of \code{\link{crm_params}}, a list of CRM
#' parameters.
#' @param fit An instance of \code{rstan::stanmodel}, derived by fitting one of
#' the trialr CRM models.
#' @return An instance of \code{\link{crm_fit}}.
crm_process <- function(dat, fit) {
  dose_indices <- seq(from = 1, to = dat$num_doses, by = 1)

  # Posterior estimates
  prob_tox_samp <- rstan::extract(fit, 'prob_tox')[[1]]
  prob_tox <- colMeans(prob_tox_samp)
  median_prob_tox <- apply(prob_tox_samp, 2, stats::median)
  recommended_dose <- which.min(abs(prob_tox - dat$target))
  model_dose <- which.min(abs(prob_tox - dat$target))

  # Implied MTD
  implied_mtd <- apply(prob_tox_samp, 1, function(x) which.min(abs(x - dat$target)))
  prob_mtd <- sapply(dose_indices, function(x) mean(implied_mtd == x))

  x <- crm_fit(dose_indices = dose_indices,
               num_patients = dat$num_patients,
               doses = dat$doses,
               tox = dat$tox,
               weights = dat$weights,
               prob_tox = prob_tox,
               median_prob_tox = median_prob_tox,
               prob_mtd = prob_mtd,
               recommended_dose = recommended_dose,
               dat = dat,
               fit = fit)
  return(x)
}
