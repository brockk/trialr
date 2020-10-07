#' @title Get parameters to run the EffTox demo
#'
#' @description Get parameters to run the EffTox demo. These match those used
#' to demonstrate EffTox in Thall et al. 2014.
#'
#' @return a \code{list} of parameters, described in \code{efftox_params}
#'
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' names(dat)
#' dat$real_doses == c(1, 2, 4, 6.6, 10)
#'
#' @seealso
#' \code{\link{efftox_params}}
#'
#' @references Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample
#' size for computing prior hyperparameters in Bayesian phase I-II dose-finding
efftox_parameters_demo <- function() {
  # Demonstration from 'Effective sample size for computing prior
  # hyperparameters in Bayesian phase I-II dose-finding', Thall et al., 2014

  p <- efftox_priors(alpha_mean = -7.9593, alpha_sd = 3.5487,
                     beta_mean = 1.5482, beta_sd = 3.5018,
                     gamma_mean = 0.7367, gamma_sd = 2.5423,
                     zeta_mean = 3.4181, zeta_sd = 2.4406,
                     eta_mean = 0, eta_sd = 0.2,
                     psi_mean = 0, psi_sd = 1)

  x <- efftox_params(real_doses = c(1, 2, 4, 6.6, 10),
                     efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                     p_e = 0.1, p_t = 0.1, eff0 = 0.5, tox1 = 0.65,
                     eff_star = 0.7, tox_star = 0.25,
                     priors = p)

  return(x)
}
