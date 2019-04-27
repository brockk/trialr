#' @title Run EffTox simulations
#'
#' @description Run EffTox simulations for assumed true efficacy and toxicity
#' curves.
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param num_sims integer, number of simulated iterations
#' @param first_dose integer, the dose-level to give to patient 1, e.g. 1 for
#' the lowest dose.
#' @param true_eff the true probabilities of efficacy at the doses under
#' investigation; a vector of numbers between 0 and 1.
#' @param true_tox the true probabilities of toxicity at the doses under
#' investigation; a vector of numbers between 0 and 1.
#' @param cohort_sizes a vector of integer cohort sizes. A dose decision is made
#' when each cohort is completed and the next cohort is treated at the
#' recommended dose. To conduct a trial using at most 20 patients, where dose is
#' re-evaluated after every second patient, use \code{rep(2, 10)}. To conduct a
#' trial of 8 patients where dose is re-evaluated after each single patient, use
#' \code{rep(1, 8)}. Cohort size need not be uniform. E.g.
#' \code{c(rep(1, 5), rep(3, 10))} represents a trial where the dose is
#' re-evaluated after each patient for the first 5 patients, and then after
#' every third patient for a further 30 patients.
#' @param ... Extra parameters provided via the ellipsis are passed to
#' \code{stan::sampling}
#'
#' @return A list with named elements \code{recommended_dose},
#' \code{efficacies}, \code{toxicities}, and \code{doses_given}.
#'
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' set.seed(123)
#' # Let's say we want to use only 2 chains. Extra args are passed to stan
#' \dontrun{
#' sims <- efftox_simulate(dat, num_sims = 10, first_dose = 1,
#'                         true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
#'                         true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
#'                         cohort_sizes = rep(3, 13),
#'                         chains = 2)
#' table(sims$recommended_dose) / length(sims$recommended_dose)
#' table(unlist(sims$doses_given)) / length(unlist(sims$doses_given))
#' table(unlist(sims$doses_given)) / length(sims$recommended_dose)
#' }
#' # In real life, we would run thousands of iterations, not 10.
#' # This is an example.
efftox_simulate <- function(dat, num_sims, first_dose, true_eff, true_tox,
                            cohort_sizes, ...) {

  recommended_dose <- integer(length = num_sims)
  efficacies <- list()
  toxicities <- list()
  doses_given <- list()

  for(i in 1:num_sims) {
    print(paste('Starting iteration', i))
    this_dat <- dat
    dose <- first_dose
    for(cohort_size in cohort_sizes) {
      prob_eff <- true_eff[dose]
      prob_tox <- true_tox[dose]
      # Simulate new efficacy events
      new_eff <- stats::rbinom(n = cohort_size, size = 1, prob = prob_eff)
      new_tox <- stats::rbinom(n = cohort_size, size = 1, prob = prob_tox)
      # And append to trial data
      this_dat$eff <- c(this_dat$eff, new_eff)
      this_dat$tox <- c(this_dat$tox, new_tox)
      # Also reflect doses delivered
      this_dat$doses <- c(this_dat$doses, rep(dose, cohort_size))
      this_dat$num_patients <- this_dat$num_patients + cohort_size
      samp <- rstan::sampling(stanmodels$EffTox, data = this_dat, ...)
      l <- efftox_process(this_dat, samp)
      # Select a dose?
      if(sum(l$acceptable) > 0) {
        # Select dose
        dose = which.max(ifelse(l$acceptable, l$utility, NA))
      } else {
        dose <- NA
        break()
      }
    }
    recommended_dose[i] = dose
    efficacies[[i]] = this_dat$eff
    toxicities[[i]] = this_dat$tox
    doses_given[[i]] = this_dat$doses
  }

  return(list(recommended_dose = recommended_dose,
              efficacies = efficacies,
              toxicities = toxicities,
              doses_given = doses_given))
}
