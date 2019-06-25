
#' @title Get the prior beliefs for a CRM trial scenario.
#'
#' @description Infer the prior beliefs consistent with the parameters and model
#' form for a CRM dose-finding trial. This function could be interpreted as
#' fitting the model to no data, thus examining the beliefs on dose-toxicity
#' that are suggested by the parameter priors alone. This function provides the
#' task analagous to \code{\link{stan_crm}} before any data has been collected.
#'
#' @param skeleton a vector of the prior guesses of toxicity at doses.
#' This should be a monotonically-increasing vector of numbers between 0 and 1.
#' @param target the target toxicity probability, a number between 0 and 1.
#' This value would normally be one of the values in \code{skeleton}, but that
#' is not a requirement.
#' @param model Character string to denote desired model. One of \code{empiric},
#' \code{logistic}, \code{logistic_gamma}, or \code{logistic2}.
#' The choice of model determines which parameters are required. See Details.
#' @param a0 Value of fixed intercept parameter.
#' Only required for certain models. See Details.
#' @param alpha_mean Prior mean of intercept variable for normal prior.
#' Only required for certain models. See Details.
#' @param alpha_sd Prior standard deviation of intercept variable for normal prior.
#' Only required for certain models. See Details.
#' @param beta_mean Prior mean of gradient variable for normal prior.
#' Only required for certain models. See Details.
#' @param beta_sd Prior standard deviation of slope variable for normal prior.
#' Only required for certain models. See Details.
#' @param beta_shape Prior shape parameter of slope variable for gamma prior.
#' Only required for certain models. See Details.
#' @param beta_inverse_scale Prior inverse scale parameter of slope variable for
#' gamma prior. Only required for certain models. See Details.
#' @param ... extra parameters passed to \code{\link{stan_crm}}.
#'
#' @details
#' Different model choices require that different parameters are
#' provided. See below.
#'
#' @section Parameter requirements of \code{empiric} model:
#' \itemize{
#'   \item \code{beta_sd}
#' }
#'
#' @section Parameter requirements of \code{logistic} model:
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_mean}
#'   \item \code{beta_sd}
#' }
#'
#' @section Parameter requirements of \code{logistic_gamma} model:
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_shape}
#'   \item \code{beta_inverse_scale}
#' }
#'
#' @section Parameter requirements of \code{logistic2} model:
#' \itemize{
#'   \item \code{alpha_mean}
#'   \item \code{alpha_sd}
#'   \item \code{beta_mean}
#'   \item \code{beta_sd}
#' }
#'
#' @return An object of class \code{\link{crm_fit}}
#'
#' @author Kristian Brock
#'
#' @references
#'   O'Quigley, J., Pepe, M., & Fisher, L. (1990).
#'   Continual reassessment method: a practical design for phase 1 clinical
#'   trials in cancer.
#'   Biometrics, 46(1), 33-48. https://www.jstor.org/stable/2531628
#'
#'   Cheung, Y.K. (2011).
#'   Dose Finding by the Continual Reassessment Method.
#'   CRC Press. ISBN 9781420091519
#'
#' @seealso
#'   \code{\link{stan_crm}}
#'   \code{\link{crm_fit}}
#'
#' @export
#'
#' @examples
#' skeleton <- c(0.05, 0.1, 0.15, 0.33, 0.5)
#' target <- 0.33
#'
#' prior_fit1 <- crm_prior_beliefs(skeleton, target, model = 'empiric',
#'                                 beta_sd = sqrt(1.34))
#' prior_fit2 <- crm_prior_beliefs(skeleton, target, model = 'logistic_gamma',
#'                                 a0 = 3, beta_shape = 1,
#'                                 beta_inverse_scale = 2)
crm_prior_beliefs <- function(skeleton, target,
                              model = c('empiric', 'logistic', 'logistic_gamma',
                                        'logistic2'),
                              a0 = NULL,
                              alpha_mean = NULL, alpha_sd = NULL,
                              beta_mean = NULL, beta_sd = NULL,
                              beta_shape = NULL, beta_inverse_scale = NULL,
                              ...) {

  stan_crm(outcome_str = '', skeleton = skeleton, target = target,
           model = model, a0 = a0, alpha_mean = alpha_mean, alpha_sd = alpha_sd,
           beta_mean = beta_mean, beta_sd = beta_sd,
           beta_shape = beta_shape, beta_inverse_scale = beta_inverse_scale,
           ...)
}
