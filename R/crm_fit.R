#' Class of model fit by \pkg{trialr} using the CRM dose-finding design.
#'
#' @name crm_fit-class
#' @aliases crm_fit
#' @docType class
#'
#' @details
#' See \code{methods(class = "crm_fit")} for an overview of available
#' methods.
#'
#' @param dose_indices A vector of integers representing the dose-levels under
#' consideration.
#' @param num_patients Integer, the number of patients analysed.
#' @param doses vector of integers representing the dose given to the patients.
#' @param tox vector of integers representing the toxicity status of the
#' patients.
#' @param weights Vector of numeric weights for the observations for patients
#' 1:num_patients, thus facilitating the TITE-CRM design.
#' @param prob_tox The posterior mean probabilities of toxicity at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param median_prob_tox The posterior median probabilities of toxicity at doses
#' 1:n; a vector of numbers between 0 and 1.
#' @param prob_mtd The posterior probability that each dose is the MTD, by the
#' chosen model; a vector of numbers between 0 and 1. This probability reflects
#' the uncertainty remaining in the parameter distributions, whereas
#' \code{prob_tox} and \code{median_prob_tox} do not.
#' @param recommended_dose An integer representing the dose-level that is
#' recommended for the next patient or cohort. Contrast to
#' \code{modal_mtd_candidate}.
#' @param dat Object \code{\link{crm_params}} containing data passed to
#' \code{\link[rstan:sampling]{sampling}}.
#' @param fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
#' @param samples An optional \code{data.frame} like object of samples.
#'
#' @export
#'
#' @seealso
#' \code{\link{stan_crm}}
crm_fit <- function(dose_indices,
                    num_patients,
                    doses,
                    tox,
                    weights,
                    prob_tox,
                    median_prob_tox,
                    prob_mtd,
                    recommended_dose,
                    dat,
                    fit,
                    samples = NULL) {

  # Elements in base class
  x <- dose_finding_fit(dose_indices = dose_indices,
                        num_patients = num_patients,
                        doses = doses,
                        tox = tox,
                        prob_tox = prob_tox,
                        median_prob_tox = median_prob_tox,
                        recommended_dose = recommended_dose,
                        dat = dat,
                        fit = fit)

  # Elements in this class
  x$weights <- weights
  x$prob_mtd <- prob_mtd
  x$modal_mtd_candidate <- which.max(prob_mtd)
  x$entropy <- .entropy(prob_mtd)

  class(x) <- c("crm_fit", "dose_finding_fit", "list")
  x
}
