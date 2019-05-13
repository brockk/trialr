#' Class of dose-finding model fit by \pkg{trialr} using Stan.
#'
#' @name dose_finding_fit-class
#' @aliases dose_finding_fit
#' @docType class
#'
#' @param dose_indices A vector of integers representing the dose-levels under
#' consideration.
#' @param num_patients Integer, the number of patients analysed.
#' @param doses vector of integers representing the dose given to the patients.
#' @param tox vector of integers representing the toxicity status of the
#' patients.
#' @param prob_tox The posterior mean probabilities of toxicity at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param median_prob_tox The posterior median probabilities of toxicity at
#' doses 1:n; a vector of numbers between 0 and 1.
#' @param recommended_dose An integer representing the dose-level that is
#' recommended for the next patient or cohort.
#' @param dat Object \code{\link{crm_params}} containing data passed to
#' \code{\link[rstan:sampling]{sampling}}.
#' @param fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
#'
#' @export
#'
#' @seealso
#' \code{\link{crm_fit}},
#' \code{\link{efftox_fit}}
dose_finding_fit <- function(dose_indices,
                             num_patients,
                             doses,
                             tox,
                             prob_tox,
                             median_prob_tox,
                             recommended_dose,
                             dat,
                             fit) {
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )

  x <- list(dose_indices = dose_indices,
            num_patients = num_patients,
            doses = doses,
            tox = tox,
            prob_tox = prob_tox,
            median_prob_tox = median_prob_tox,
            recommended_dose = recommended_dose,
            dat = dat,
            fit = fit,
            version = version)

  class(x) <- c("dose_finding_fit", "list")
  x
}
