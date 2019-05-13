
#' Fit an EffTox model
#'
#' Fit an EffTox model using Stan for full Bayesian inference.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{efftox_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given}, \code{eff} and
#' \code{tox} parameters. See Details.
#' @param real_doses A vector of numbers.The doses under investigation. They
#' should be ordered from lowest to highest and be in consistent units.
#' E.g., #' to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
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
#' @param zeta_mean The prior normal mean of the slope term in the efficacy logit
#' model. A number.
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
#' @param doses_given A optional vector of dose-levels given to patients
#' 1:num_patients, where 1=lowest dose, 2=second dose, etc. Only required when
#' \code{outcome_str} is not provided.
#' @param eff An optional vector of efficacy outcomes for patients
#' 1:num_patients, where 1=efficacy and 0=no efficacy. Only required when
#' \code{outcome_str} is not provided.
#' @param tox An optional vector of toxicity outcomes for patients
#' 1:num_patients, where 1=toxicity and 0=no toxicity. Only required when
#' \code{outcome_str} is not provided.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. \code{\link[rstan:sampling]{sampling}}.
#'
#' @details
#' The quickest and easiest way to fit an EffTox model to some observed outcomes
#' is to describe the outcomes using \pkg{trialr}'s syntax for efficacy-toxicity
#' dose-finding outcomes. See \code{\link{efftox_parse_outcomes}} for full
#' details and examples.
#'
#' Utility or attractivess scores are calculated in EffTox using L^p norms.
#' Imagine the first quadrant of a scatter plot with prob_eff along the x-axis
#' and prob_tox along the y-axis.
#' The point (1, 0) (i.e. guaranteed efficacy & no toxicity) is the holy grail.
#' The neutral contour intersects the points (eff0, 0), (1, tox1) and
#' (eff_star, tox_star). A unique curve intersects these three points and
#' identifies a value for p, the exponent in the L^p norm. On this neutral-
#' utility contour, scores are equal to zero. A family of curves with different
#' utility scores is defined that are "parallel" to this neutral curve.
#' Points with probabilities of efficacy and toxicity that are nearer to (1, 0)
#' will yield greater scores, and vice-versa.
#'
#' @return An object of class \code{\link{efftox_fit}}
#'
#' @author Kristian Brock \email{kristian.brock@@gmail.com}
#'
#' @references
#'   Thall, P., & Cook, J. (2004). Dose-Finding Based on Efficacy-Toxicity
#'     Trade-Offs. Biometrics, 60(3), 684-693.
#'
#'   Thall, P., Herrick, R., Nguyen, H., Venier, J., & Norris, J. (2014).
#'     Effective sample size for computing prior hyperparameters in Bayesian
#'     phase I-II dose-finding. Clinical Trials, 11(6), 657-666.
#'     https://doi.org/10.1177/1740774514547397
#'
#'   Brock, K., Billingham, L., Copland, M., Siddique, S., Sirovica, M., &
#'     Yap, C. (2017). Implementing the EffTox dose-finding design in the
#'     Matchpoint trial. BMC Medical Research Methodology, 17(1), 112.
#'     https://doi.org/10.1186/s12874-017-0381-x
#'
#' @seealso
#'   \code{\link{efftox_fit}}
#'   \code{\link{stan_efftox_demo}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This model is presented in Thall et al. (2014)
#' mod1 <- stan_efftox('1N 2E 3B',
#'                      real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'                      efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
#'                      p_e = 0.1, p_t = 0.1,
#'                      eff0 = 0.5, tox1 = 0.65,
#'                      eff_star = 0.7, tox_star = 0.25,
#'                      alpha_mean = -7.9593, alpha_sd = 3.5487,
#'                      beta_mean = 1.5482, beta_sd = 3.5018,
#'                      gamma_mean = 0.7367, gamma_sd = 2.5423,
#'                      zeta_mean = 3.4181, zeta_sd = 2.4406,
#'                      eta_mean = 0, eta_sd = 0.2,
#'                      psi_mean = 0, psi_sd = 1, seed = 123)
#'
#' # Shorthand for the above is:
#' mod2 <- stan_efftox_demo('1N 2E 3B', seed = 123)
#'
#' # the seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
#' }
stan_efftox <- function(outcome_str = NULL,
                        real_doses, efficacy_hurdle, toxicity_hurdle, p_e, p_t,
                        eff0, tox1, eff_star, tox_star,
                        alpha_mean, alpha_sd, beta_mean, beta_sd,
                        gamma_mean, gamma_sd, zeta_mean, zeta_sd,
                        eta_mean, eta_sd, psi_mean, psi_sd,
                        doses_given = NULL,
                        eff = NULL,
                        tox = NULL,
                        ...) {

  # Create parameters object to pass to Stan
  dat <- efftox_params(real_doses, efficacy_hurdle, toxicity_hurdle,
                       p_e, p_t, eff0, tox1, eff_star, tox_star,
                       alpha_mean, alpha_sd, beta_mean, beta_sd,
                       gamma_mean, gamma_sd, zeta_mean, zeta_sd,
                       eta_mean, eta_sd, psi_mean, psi_sd)

  # Add outcomes
  if(is.null(outcome_str)) {
    if(length(doses_given) != length(eff))
      stop('doses_given and eff vectors should have same length')
    if(length(tox) != length(eff))
      stop('tox and eff vectors should have same length')
    dat$doses <- array(doses_given)
    dat$eff <- array(eff)
    dat$tox <- array(tox)
    dat$num_patients <- length(doses_given)
  } else {
    outcomes_df <- efftox_parse_outcomes(outcome_str, as.list = TRUE)
    dat$num_patients <- outcomes_df$num_patients
    dat$doses <- array(outcomes_df$doses)
    dat$eff <- array(outcomes_df$eff)
    dat$tox <- array(outcomes_df$tox)
  }

  # Fit data to model using Stan
  samp <- rstan::sampling(stanmodels$EffTox, data = dat, ...)
  # Create useful output from posterior samples
  decision <- efftox_process(dat, samp)

  return(decision)
}
