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
#' @param modal_mtd_candidate An integer representing the dose-level most likely
#' to be the MTD, i.e. the dose-level that maximises \code{prob_mtd}.
#' @param entropy The value of the entropy implied by prob_mtd. If the values of
#' prob_mtd are equal, entropy is maximised, taking value log(n) where n is
#' the number of doses, and we could not be more ignorant about the identity of
#' the maximum tolerable dose. If prob_mtd is a vector of zeroes except for a
#' single value of 1, then entropy is minimisied, taking value 0, and we are
#' sure of the identity of the maximum tolerable dose.
#' @param dat Object \code{\link{crm_params}} containing data passed to
#' \code{\link[rstan:sampling]{sampling}}.
#' @param fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
#' @param An optional \code{data.frame} like object of samples.
#'
#' @export
#'
#' @seealso
#' \code{\link{stan_crm}}
crm_fit <- function(dose_indices,
                    prob_tox,
                    median_prob_tox,
                    prob_mtd,
                    recommended_dose,
                    # modal_mtd_candidate,
                    dat,
                    fit,
                    samples = NULL) {
  # crm_fit class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )
  x <- list(dose_indices = dose_indices,
            prob_tox = prob_tox,
            median_prob_tox = median_prob_tox,
            prob_mtd = prob_mtd,
            recommended_dose = recommended_dose,
            # modal_mtd_candidate = modal_mtd_candidate,
            modal_mtd_candidate = which.max(prob_mtd),
            entropy = .entropy(prob_mtd),
            dat = dat,
            fit = fit,
            samples = samples,
            version = version)
  class(x) <- "crm_fit"
  x
}
