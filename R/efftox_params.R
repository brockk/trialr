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
#' @param alpha_mean The prior normal mean of the intercept term in the toxicity
#' logit model. A number.
#' @param alpha_sd The prior normal standard deviation of the intercept term in
#' the toxicity logit model. A number.
#' @param beta_mean The prior normal mean of the slope term in the toxicity
#' logit model. A number.
#' @param beta_sd The prior normal standard deviation of the slope term in the
#' toxicity logit model. A number.
#' @param gamma_mean The prior normal mean of the intercept term in the efficacy
#' logit model. A number.
#' @param gamma_sd The prior normal standard deviation of the intercept term in
#' the efficacy logit model. A number.
#' @param zeta_mean The prior normal mean of the slope term in the efficacy
#' logit model. A number.
#' @param zeta_sd The prior normal standard deviation of the slope term in the
#' efficacy logit model. A number.
#' @param eta_mean The prior normal mean of the squared term coefficient in the
#' efficacy logit model. A number.
#' @param eta_sd The prior normal standard deviation of the squared term
#' coefficient in the efficacy logit model. A number.
#' @param psi_mean The prior normal mean of the association term in the combined
#' efficacy-toxicity model. A number.
#' @param psi_sd The prior normal standard deviation of the association term in
#' the combined efficacy-toxicity model. A number.
#'
#' @export
#'
#' @seealso
#' \code{\link{stan_efftox}}
#' \code{\link{stan_efftox_demo}}
efftox_params <- function(real_doses, efficacy_hurdle, toxicity_hurdle,
                          p_e, p_t, eff0, tox1, eff_star, tox_star,
                          alpha_mean, alpha_sd, beta_mean, beta_sd,
                          gamma_mean, gamma_sd, zeta_mean, zeta_sd,
                          eta_mean, eta_sd, psi_mean, psi_sd) {
  # efftox_params class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )

  p <- efftox_solve_p(eff0, tox1, eff_star, tox_star)
  x <- list(real_doses = real_doses, num_doses = length(real_doses),
            efficacy_hurdle = efficacy_hurdle, toxicity_hurdle = toxicity_hurdle,
            p_e = p_e, p_t = p_t, p = p, eff0 = eff0, tox1 = tox1,
            eff_star = eff_star, tox_star = tox_star,
            alpha_mean = alpha_mean, alpha_sd = alpha_sd,
            beta_mean = beta_mean, beta_sd = beta_sd,
            gamma_mean = gamma_mean, gamma_sd = gamma_sd,
            zeta_mean = zeta_mean, zeta_sd = zeta_sd,
            eta_mean = eta_mean, eta_sd = eta_sd,
            psi_mean = psi_mean, psi_sd = psi_sd,
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
