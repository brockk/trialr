
#' Fit an EffTox model
#'
#' Fit an EffTox model using Stan for full Bayesian inference.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{efftox_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given}, \code{eff} and
#' \code{tox} parameters. See Details.
#' @param real_doses A vector of numbers, the doses under investigation. They
#' should be ordered from lowest to highest and be in consistent units.
#' E.g. to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
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
#' @param alpha_mean Optional, the prior normal mean of the intercept term in
#' the toxicity logit model. A number. You should prioritise specifying this
#' value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param alpha_sd Optional, the prior normal standard deviation of the
#' intercept term in the toxicity logit model. A number.You should prioritise
#' specifying this value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param beta_mean Optional, the prior normal mean of the slope term in the
#' toxicity logit model. A number. You should prioritise specifying this
#' value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param beta_sd Optional, the prior normal standard deviation of the slope
#' term in the toxicity logit model. A number. You should prioritise specifying
#' this value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param gamma_mean Optional, The prior normal mean of the intercept term in
#' the efficacy logit model. A number. You should prioritise specifying this
#' value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param gamma_sd Optional, the prior normal standard deviation of the
#' intercept term in the efficacy logit model. A number. You should prioritise
#' specifying this value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param zeta_mean Optional, the prior normal mean of the slope term in the
#' efficacy logit model. A number. You should prioritise specifying this value
#' via  \code{priors} but this option is provided for backwards-compatibility.
#' @param zeta_sd Optional, the prior normal standard deviation of the slope
#' term in the efficacy logit model. A number. You should prioritise specifying
#' this value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param eta_mean Optional, the prior normal mean of the squared term
#' coefficient in the efficacy logit model. A number. You should prioritise
#' specifying this value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param eta_sd Optional, the prior normal standard deviation of the squared
#' term coefficient in the efficacy logit model. A number. You should prioritise
#' specifying this value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param psi_mean Optional, the prior normal mean of the association term in
#' the combined efficacy-toxicity model. A number. You should prioritise
#' specifying this value via \code{priors} but this option is provided for
#' backwards-compatibility.
#' @param psi_sd Optional, the prior normal standard deviation of the
#' association term in the combined efficacy-toxicity model. A number. You
#' should prioritise specifying this value via \code{priors} but this option is
#' provided for backwards-compatibility.
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
#' @author Kristian Brock \email{kristian.brock@gmail.com}
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
#' \code{\link{efftox_priors}}
#' \code{\link{get_efftox_priors}}
#' \code{\link{efftox_fit}}
#' \code{\link{stan_efftox_demo}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This model is presented in Thall et al. (2014).
#'p <- efftox_priors(alpha_mean = -7.9593, alpha_sd = 3.5487,
#'                   beta_mean = 1.5482, beta_sd = 3.5018,
#'                   gamma_mean = 0.7367, gamma_sd = 2.5423,
#'                   zeta_mean = 3.4181, zeta_sd = 2.4406,
#'                   eta_mean = 0, eta_sd = 0.2,
#'                   psi_mean = 0, psi_sd = 1)
#' mod1 <- stan_efftox('1N 2E 3B',
#'                     real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'                     efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
#'                     p_e = 0.1, p_t = 0.1,
#'                     eff0 = 0.5, tox1 = 0.65,
#'                     eff_star = 0.7, tox_star = 0.25,
#'                     priors = p,
#'                     seed = 123)
#'
#' # The above is a longhad version of:
#' mod2 <- stan_efftox_demo('1N 2E 3B', seed = 123)
#'
#' # the seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
#' }
stan_efftox <- function(outcome_str = NULL,
                        real_doses, efficacy_hurdle, toxicity_hurdle, p_e, p_t,
                        eff0, tox1, eff_star, tox_star,
                        priors = NULL,
                        alpha_mean = NULL, alpha_sd = NULL,
                        beta_mean = NULL, beta_sd = NULL,
                        gamma_mean = NULL, gamma_sd = NULL,
                        zeta_mean = NULL, zeta_sd = NULL,
                        eta_mean = NULL, eta_sd = NULL,
                        psi_mean = NULL, psi_sd = NULL,
                        doses_given = NULL,
                        eff = NULL,
                        tox = NULL,
                        ...) {

  if(is.null(priors)) {
    priors <- efftox_priors(alpha_mean = alpha_mean, alpha_sd = alpha_sd,
                            beta_mean = beta_mean, beta_sd = beta_sd,
                            gamma_mean = gamma_mean, gamma_sd = gamma_sd,
                            zeta_mean = zeta_mean, zeta_sd = zeta_sd,
                            eta_mean = eta_mean, eta_sd = eta_sd,
                            psi_mean = psi_mean, psi_sd = psi_sd)
  }

  # Create parameters object to pass to Stan
  dat <- efftox_params(real_doses, efficacy_hurdle, toxicity_hurdle,
                       p_e, p_t, eff0, tox1, eff_star, tox_star, priors)

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
