#' Container class for parameters to fit the CRM models in trialr.
#'
#' @name crm_params-class
#' @aliases crm_params
#' @docType class
#'
#' @param skeleton a vector of the prior guesses of toxicity at doses.
#' This should be a monotonically-increasing vector of numbers between 0 and 1.
#' @param target the target toxicity probability, a number between 0 and 1.
#' This value would normally be one of the values in \code{skeleton}, but that
#' is not a requirement.
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
#'
#' @details
#' Different model parameterisations require that difference parameter values
#' are specified.
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
#' @export
#'
#' @seealso
#' \code{\link{stan_crm}}
crm_params <- function(skeleton, target, a0 = NULL,
                       alpha_mean = NULL, alpha_sd = NULL,
                       beta_mean = NULL, beta_sd = NULL,
                       beta_shape = NULL, beta_inverse_scale = NULL) {

  # crm_params class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )

  x <- list(num_doses = length(skeleton), skeleton = skeleton, target = target,
            a0 = a0, alpha_mean = alpha_mean, alpha_sd = alpha_sd,
            beta_mean = beta_mean, beta_sd = beta_sd,
            beta_shape = beta_shape, beta_inverse_scale = beta_inverse_scale,
            version = version)

  # Initialise with no patients observed
  x$num_patients = 0
  x$doses = integer(length = 0)
  x$tox = integer(length = 0)
  x$weights = numeric(length = 0)

  # Set type. This is, at heart, just a list.
  class(x) <- c("crm_params", "list")
  return(x)
}
