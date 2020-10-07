#' Container class for parameters to fit the EffTox model in trialr.
#'
#' @name efftox_params-class
#' @aliases efftox_params
#' @docType class
#'
#' @param real_doses a vector of numbers.The doses under investigation.
#' They should be ordered from lowest to highest and be in consistent units.
#' E.g., to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
#' c(10, 20, 50).
#' @param efficacy_hurdle Minimum acceptable efficacy probability.
#' A number between 0 and 1.
#' @param toxicity_hurdle Maximum acceptable toxicity probability.
#' A number between 0 and 1.
#' @param p_e Certainty required to infer a dose is acceptable with regards to
#' being probably efficacious; a number between 0 and 1.
#' @param p_t Certainty required to infer a dose is acceptable with regards to
#' being probably tolerable; a number between 0 and 1.
#' @param eff0 Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1 (see Details).
#' @param tox1 Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1 (see Details).
#' @param eff_star Efficacy probability of an equi-utility third point (see
#' Details).
#' @param tox_star Toxicity probability of an equi-utility third point (see
#' Details).
#' @param priors instance of class \code{\link{efftox_priors}}, the
#' hyperparameters for normal priors on the six model parameters.
#'
#' @export
#'
#' @seealso
#' \code{\link{efftox_priors}}
#' \code{\link{get_efftox_priors}}
#' \code{\link{stan_efftox}}
#' \code{\link{stan_efftox_demo}}
efftox_params <- function(real_doses, efficacy_hurdle, toxicity_hurdle,
                          p_e, p_t, eff0, tox1, eff_star, tox_star,
                          priors) {

  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )

  p <- efftox_solve_p(eff0, tox1, eff_star, tox_star)
  x <- list(real_doses = real_doses, num_doses = length(real_doses),
            efficacy_hurdle = efficacy_hurdle, toxicity_hurdle = toxicity_hurdle,
            p_e = p_e, p_t = p_t, p = p, eff0 = eff0, tox1 = tox1,
            eff_star = eff_star, tox_star = tox_star,
            alpha_mean = priors$alpha_mean, alpha_sd = priors$alpha_sd,
            beta_mean = priors$beta_mean, beta_sd = priors$beta_sd,
            gamma_mean = priors$gamma_mean, gamma_sd = priors$gamma_sd,
            zeta_mean = priors$zeta_mean, zeta_sd = priors$zeta_sd,
            eta_mean = priors$eta_mean, eta_sd = priors$eta_sd,
            psi_mean = priors$psi_mean, psi_sd = priors$psi_sd,
            version = version)

  # Initialise with no patients observed
  x$doses = c()
  x$eff = c()
  x$tox = c()
  x$num_patients = 0

  # Set type. This is, at heart, just a list.
  class(x) <- c("efftox_params", "list")
  return(x)
}
