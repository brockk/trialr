#' Fit a CRM model
#'
#' Fit a continual reassessment method (CRM) model for dose-finding using Stan
#' for full Bayesian inference. There are several likelihood and prior
#' combinations supported. See model-specific sections below.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{df_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given} and \code{tox}
#' parameters. See Details.
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
#' @param doses_given A optional vector of dose-levels given to patients
#' 1:num_patients, where 1=lowest dose, 2=second dose, etc. Only required when
#' \code{outcome_str} is not provided.
#' @param tox An optional vector of toxicity outcomes for patients
#' 1:num_patients, where 1=toxicity and 0=no toxicity. Only required when
#' \code{outcome_str} is not provided.
#' @param weights An optional vector of numeric weights for the observations
#' for patients 1:num_patients, thus facilitating the TITE-CRM design.
#' Can be used with \code{outcome_str}, or with \code{doses_given} and
#' \code{tox}. It is generally tidier to specify \code{doses_given},
#' \code{tox} and \code{weights} when a TITE-CRM analysis is desired.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores}, and
#' \code{control}.
#'
#' @details
#' The quickest and easiest way to fit a CRM model to some observed outcomes
#' is to describe the outcomes using \pkg{trialr}'s syntax for dose-finding
#' outcomes. See \code{\link{df_parse_outcomes}} for full details and examples.
#'
#' Different model choices require that different parameters are
#' provided. See sections below.
#'
#' @section The \code{empiric} model:
#' The model form is:
#'
#' \eqn{F(x_{i}, \beta) = x_{i}^{\exp{\beta}}}
#'
#' and the required parameters are:
#'
#' \itemize{
#'   \item \code{beta_sd}
#' }
#'
#' @section The \code{logistic} model:
#' The model form is:
#'
#' \eqn{F(x_{i}, \beta) = 1 / (1 + \exp{(-a_{0} - \exp{(\beta)} x_{i}})) }
#'
#' and the required parameters are:
#'
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_mean}
#'   \item \code{beta_sd}
#' }
#'
#' @section The \code{logistic_gamma} model:
#' The model form is:
#'
#' \eqn{F(x_{i}, \beta) = 1 / (1 + \exp{(-a_{0} - \exp{(\beta)} x_{i}})) }
#'
#' and the required parameters are:
#'
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_shape}
#'   \item \code{beta_inverse_scale}
#' }
#'
#' @section The \code{logistic2} model:
#' The model form is:
#'
#' \eqn{F(x_{i}, alpha, \beta) = 1 / (1 + \exp{(-\alpha - \exp{(\beta)} x_i)}) }
#'
#' and the required parameters are:
#'
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
#'   \code{\link{crm_fit}}
#'   \code{\link[rstan:sampling]{sampling}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # CRM example
#' fit1 <- stan_crm('1N 2N 3T', skeleton = c(0.1, 0.2, 0.35, 0.6),
#'                  target = 0.2, model = 'empiric', beta_sd = sqrt(1.34),
#'                  seed = 123)
#'
#' fit2 <- stan_crm('1NNN 2NNN 3TTT', skeleton = c(0.1, 0.2, 0.35, 0.6),
#'                  target = 0.2, model = 'logistic', a0 = 3, beta_mean = 0,
#'                  beta_sd = sqrt(1.34), seed = 123)
#'
#' # The seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
#'
#' # TITE-CRM example, p.124 of Dose Finding by the CRM, Cheung (2010)
#' fit3 <-stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
#'                 doses_given = c(3, 3, 3, 3),
#'                 tox = c(0, 0, 0, 0),
#'                 weights = c(73, 66, 35, 28) / 126,
#'                 model = 'empiric', beta_sd = sqrt(1.34), seed = 123)
#' fit3$recommended_dose
#' }
stan_crm <- function(outcome_str = NULL, skeleton, target,
                     model = c('empiric', 'logistic', 'logistic_gamma',
                               'logistic2'),
                     a0 = NULL,
                     alpha_mean = NULL, alpha_sd = NULL,
                     beta_mean = NULL, beta_sd = NULL,
                     beta_shape = NULL, beta_inverse_scale = NULL,
                     doses_given = NULL,
                     tox = NULL,
                     weights = NULL,
                     ...) {

  model <- match.arg(model)

  # CRM assumes a strictly monotonically increasing skeleton.
  # Stop if this is not the case.
  if(length(skeleton) > 1) {
    if(any(skeleton[-1] <= skeleton[-length(skeleton)])) {
      stop('Skeleton must be strictly monotonically-increasing.')
    }
  }

  # Create parameters object to pass to Stan
  dat <- crm_params(skeleton = skeleton,
                    target = target,
                    a0 = a0,
                    alpha_mean = alpha_mean,
                    alpha_sd = alpha_sd,
                    beta_mean = beta_mean,
                    beta_sd = beta_sd,
                    beta_shape = beta_shape,
                    beta_inverse_scale = beta_inverse_scale)

  # Add outcomes
  if(is.null(outcome_str)) {
    if(length(doses_given) != length(tox))
      stop('doses_given and tox vectors should have same length')

    dat$doses <- array(doses_given)
    dat$tox <- array(tox)
    dat$num_patients <- length(doses_given)
  } else {
    outcomes_df <- df_parse_outcomes(outcome_str, as.list = TRUE)
    dat$num_patients <- outcomes_df$num_patients
    dat$doses <- array(outcomes_df$doses)
    dat$tox <- array(outcomes_df$tox)
  }
  # Add weights if specified; infer all to be 1 if not.
  if(is.null(weights))
    dat$weights <- array(rep(1, dat$num_patients))
  else
    dat$weights <- array(weights)

  # Fit data to model using Stan, after performing model-specific checks.
  if(model == 'empiric') {
    # Check parameters
    if(is.null(beta_sd)) stop('beta_sd parameter must be specified.')
    if(beta_sd <= 0) stop('beta_sd parameter must be strictly positive.')
    # Fit
    samp <- rstan::sampling(stanmodels$CrmEmpiricNormalPrior,
                            data = dat, ...)
  } else if(model == 'logistic') {
    # Check parameters
    if(is.null(a0)) stop('a0 parameter must be specified.')
    if(is.null(beta_mean)) stop('beta_mean parameter must be specified.')
    if(is.null(beta_sd)) stop('beta_sd parameter must be specified.')
    if(beta_sd <= 0) stop('beta_sd parameter must be strictly positive.')
    # Fit
    samp <- rstan::sampling(stanmodels$CrmOneParamLogisticNormalPrior,
                            data = dat, ...)
  } else if(model == 'logistic_gamma') {
    # Check parameters
    if(is.null(a0)) stop('a0 parameter must be specified.')
    if(is.null(beta_shape)) stop('beta_shape parameter must be specified.')
    if(is.null(beta_inverse_scale)) {
      stop('beta_inverse_scale parameter must be specified.')
    }
    if(beta_shape <= 0) stop('beta_shape parameter must be strictly positive.')
    if(beta_inverse_scale <= 0) {
      stop('beta_inverse_scale parameter must be strictly positive.')
    }
    # Fit
    samp <- rstan::sampling(stanmodels$CrmOneParamLogisticGammaPrior,
                            data = dat, ...)
  } else if(model == 'logistic2') {
    # Check parameters
    if(is.null(alpha_mean)) stop('alpha_mean parameter must be specified.')
    if(is.null(alpha_sd)) stop('alpha_sd parameter must be specified.')
    if(alpha_sd <= 0) stop('alpha_sd parameter must be strictly positive.')
    if(is.null(beta_mean)) stop('beta_mean parameter must be specified.')
    if(is.null(beta_sd)) stop('beta_sd parameter must be specified.')
    if(beta_sd <= 0) stop('beta_sd parameter must be strictly positive.')
    # Fit
    samp <- rstan::sampling(stanmodels$CrmTwoParamLogisticNormalPrior,
                            data = dat, ...)
  } else {
    stop(paste0("Model type '", model, "' is not recognised. See ?stan_crm"))
  }

  # Create useful output from posterior samples
  decision <- crm_process(dat, samp)

  return(decision)
}
