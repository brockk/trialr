
#' Class of model fit by \pkg{trialr} using the EffTox dose-finding design.
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
#' @param recommended_dose An integer representing the dose-level recommended
#' for the next patient or cohort; or \code{NA} if stopping is recommended.
#' @param prob_eff The posterior mean probabilities of efficacy at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param prob_tox The posterior mean probabilities of toxicity at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param prob_acc_eff The posterior mean probabilities that efficacy at the
#' doses is acceptable, i.e. that it exceeds the minimum acceptable efficacy
#' threshold; a vector of numbers between 0 and 1.
#' @param prob_acc_tox The posterior mean probabilities that toxicity at the
#' doses is acceptable, i.e. that it is less than the maximum toxicity
#' threshold; a vector of numbers between 0 and 1.
#' @param utility The utilities of doses 1:n, calculated by plugging the
#' posterior mean probabilities of efficacy and toxicity into the utility
#' formula, as advocated by Thall & Cook. Contrast to \code{post_utility};
#' a vector of numbers.
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
efftox_fit <- function(dose_indices, recommended_dose, prob_eff, prob_tox,
                       prob_acc_eff, prob_acc_tox, utility, post_utility,
                       acceptable, dat, fit) {
  # efftox_fit class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )
  x <- list(dose_indices = dose_indices, recommended_dose = recommended_dose,
            prob_eff = prob_eff, prob_tox = prob_tox,
            prob_acc_eff = prob_acc_eff, prob_acc_tox = prob_acc_tox,
            utility = utility, post_utility = post_utility,
            acceptable = acceptable,
            dat = dat, fit = fit, version = version)
  class(x) <- "efftox_fit"
  x
}
