
#' Simple class to hold prior hyperparameters for the EffTox model.
#'
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
#'
#' @return list-like, instance of class \code{efftox_priors}.
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
#' @export
#'
#' @examples
#' # The priors used in Thall et al. (2014)
#' p <- efftox_priors(alpha_mean = -7.9593, alpha_sd = 3.5487,
#'                    beta_mean = 1.5482, beta_sd = 3.5018,
#'                    gamma_mean = 0.7367, gamma_sd = 2.5423,
#'                    zeta_mean = 3.4181, zeta_sd = 2.4406,
#'                    eta_mean = 0, eta_sd = 0.2,
#'                    psi_mean = 0, psi_sd = 1)
#' # The class exists simply to hold these twelve values.
efftox_priors <- function(alpha_mean, alpha_sd,
                         beta_mean, beta_sd,
                         gamma_mean, gamma_sd,
                         zeta_mean, zeta_sd,
                         eta_mean, eta_sd,
                         psi_mean, psi_sd) {

  l <- list(alpha_mean = alpha_mean, alpha_sd = alpha_sd,
            beta_mean = beta_mean, beta_sd = beta_sd,
            gamma_mean = gamma_mean, gamma_sd = gamma_sd,
            zeta_mean = zeta_mean, zeta_sd = zeta_sd,
            eta_mean = eta_mean, eta_sd = eta_sd,
            psi_mean = psi_mean, psi_sd = psi_sd)
  class(l) <- c('efftox_priors')
  l
}
