
#' Class of model fit by \pkg{trialr} using the EffTox dose-finding design.
#'
#' Phase I/II dose-finding trials, i.e. those that search for a dose my efficacy
#' and toxicity outcomes search for the optimal biological dose (OBD), rather
#' than the maximum tolerated dose (MTD) that is typically sought be traditional
#' toxicity-only dose-finding.
#'
#' @name efftox_fit-class
#' @aliases efftox_fit
#' @docType class
#'
#' @details
#' See \code{methods(class = "efftox_fit")} for an overview of available
#' methods.
#'
#' @param dose_indices A vector of integers representing the dose-levels under
#' consideration.
#' @param num_patients Integer, the number of patients analysed.
#' @param doses vector of integers representing the dose given to the patients.
#' @param tox vector of integers representing the toxicity status of the
#' patients.
#' @param eff vector of integers representing the efficacy status of the
#' patients.
#' @param recommended_dose An integer representing the dose-level recommended
#' for the next patient or cohort; or \code{NA} if stopping is recommended.
#' @param prob_tox The posterior mean probabilities of toxicity at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param prob_eff The posterior mean probabilities of efficacy at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param median_prob_tox The posterior median probabilities of toxicity at
#' doses 1:n; a vector of numbers between 0 and 1.
#' @param median_prob_eff The posterior mean probabilities of efficacy at doses
#' 1:n; a vector of numbers between 0 and 1.
#' @param prob_acc_tox The posterior mean probabilities that toxicity at the
#' doses is acceptable, i.e. that it is less than the maximum toxicity
#' threshold; a vector of numbers between 0 and 1.
#' @param prob_acc_eff The posterior mean probabilities that efficacy at the
#' doses is acceptable, i.e. that it exceeds the minimum acceptable efficacy
#' threshold; a vector of numbers between 0 and 1.
#' @param utility The utilities of doses 1:n, calculated by plugging the
#' posterior mean probabilities of efficacy and toxicity into the utility
#' formula, as advocated by Thall & Cook. Contrast to \code{post_utility};
#' a vector of numbers.
#' @param prob_obd The posterior probability that each dose is the optimal
#' biological dose (OBD); a vector of numbers between 0 and 1. This probability
#' reflects the uncertainty remaining in the parameter distributions, whereas
#' \code{prob_tox} and \code{prob_eff} (etc) do not.
#' @param post_utility The posterior mean utilities of doses 1:n, calculated
#' from the posterior distributions of the utilities. This is in contrast to
#' \code{utility}, which uses plug-in posterior means of efficacy and toxicity,
#' as advocated by Thall & Cook; a vector of numbers.
#' @param acceptable A vector of logical values to indicate whether doses 1:n
#' are acceptable, according to the rules for acceptable efficacy & toxicity,
#' and rules on not skipping untested doses.
#' @param fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
#' @param dat Object \code{\link{efftox_params}} containing data passed to
#' \code{\link[rstan:sampling]{sampling}}.
#'
#' @export
#'
#' @seealso
#' \code{\link{stan_efftox}}
#' \code{\link{stan_efftox_demo}}
efftox_fit <- function(dose_indices,
                       num_patients,
                       doses,
                       tox,
                       eff,
                       prob_tox,
                       prob_eff,
                       median_prob_tox,
                       median_prob_eff,
                       prob_acc_tox,
                       prob_acc_eff,
                       utility,
                       post_utility,
                       prob_obd,
                       acceptable,
                       recommended_dose,
                       dat,
                       fit) {

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
  x$eff <- eff
  x$prob_eff <- prob_eff
  x$median_prob_eff <- median_prob_eff
  x$prob_acc_tox <- prob_acc_tox
  x$prob_acc_eff <- prob_acc_eff
  x$utility <- utility
  x$post_utility <- post_utility
  x$prob_obd <- prob_obd
  x$modal_obd_candidate <- which.max(prob_obd)
  x$entropy <- .entropy(prob_obd)
  x$acceptable <- acceptable

  class(x) <- c("efftox_fit", "dose_finding_fit", "list")
  x
}
