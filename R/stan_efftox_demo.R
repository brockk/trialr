
#' Fit the EffTox model presented in Thall et al. (2014)
#'
#' Fit the EffTox model presented in Thall et al. (2014) using Stan for full
#' Bayesian inference.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{efftox_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given}, \code{eff} and
#' \code{tox} parameters. See Details.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. \code{\link[rstan:sampling]{sampling}}.
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
#'   \code{\link{efftox_fit}}
#'   \code{\link{stan_efftox}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This model is presented in Thall et al. (2014)
#' mod2 <- stan_efftox_demo('1N 2E 3B', seed = 123)
#'
#' # The seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
#' }
stan_efftox_demo <- function(outcome_str, ...) {
  p <- efftox_priors(alpha_mean = -7.9593, alpha_sd = 3.5487,
                     beta_mean = 1.5482, beta_sd = 3.5018,
                     gamma_mean = 0.7367, gamma_sd = 2.5423,
                     zeta_mean = 3.4181, zeta_sd = 2.4406,
                     eta_mean = 0, eta_sd = 0.2,
                     psi_mean = 0, psi_sd = 1)
  stan_efftox(outcome_str,
              real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
              efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
              p_e = 0.1, p_t = 0.1,
              eff0 = 0.5, tox1 = 0.65,
              eff_star = 0.7, tox_star = 0.25,
              priors = p,
              ...)
}
