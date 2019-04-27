
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
#' @param prior_samples Integer, number of prior samples to draw. 4000 by
#' default to match rstan.
#'
#' @details
#' Different model choices require that different parameters are
#' provided. See below.
#'
#' @section Requirements of \code{empiric} model:
#' \itemize{
#'   \item \code{beta_sd}
#' }
#'
#' @section Requirements of \code{logistic} model:
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_mean}
#'   \item \code{beta_sd}
#' }
#'
#' @section Requirements of \code{logistic_gamma} model:
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_shape}
#'   \item \code{beta_inverse_scale}
#' }
#'
#' @section Requirements of \code{logistic2} model:
#' \itemize{
#'   \item \code{a0}
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
#'   Biometrics, 46(1), 33-48. https://doi.org/10.2307/2531628
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
                              prior_samples = 4000) {

  model <- match.arg(model)
  dose_indices <- 1:length(skeleton)

  closest_to_target <- function(prob_tox, target) {
    which.min(abs(prob_tox - target))
  }
  select_dose_func = closest_to_target

  # Create parameters object to pass to Stan
  dat <- crm_params(skeleton = skeleton,
                    target = target,
                    a0 = a0,
                    alpha_mean = alpha_mean,
                    alpha_sd = alpha_sd,
                    beta_mean = beta_mean,
                    beta_sd = beta_sd,
                    beta_shape = beta_shape,
                    beta_inverse_scale = beta_inverse_scale
  )

  # CRM assumes a strictly monotonically increasing skeleton.
  # Stop if this is not the case.
  if(length(skeleton) > 1) {
    if(any(skeleton[-1] <= skeleton[-length(skeleton)])) {
      stop('Skeleton must be strictly monotonically-increasing.')
    }
  }

  if(model == 'empiric') {
    # Check parameters
    if(is.null(beta_sd)) stop('beta_sd parameter must be specified.')
    if(beta_sd <= 0) stop('beta_sd parameter must be strictly positive.')
    # Sample
    coded_doses <- skeleton
    beta_sample <- rnorm(n = prior_samples, mean = 0, sd = beta_sd)
    prior_curves <- t(sapply(beta_sample,
                             function(beta) coded_doses^(exp(beta))))
    samples <- cbind(beta_sample, prior_curves)
    colnames(samples) <- c('beta', paste0('prob_tox[', dose_indices, ']'))
  } else if(model == 'logistic') {
    # Check parameters
    if(is.null(a0)) stop('a0 parameter must be specified.')
    if(is.null(beta_mean)) stop('beta_mean parameter must be specified.')
    if(is.null(beta_sd)) stop('beta_sd parameter must be specified.')
    if(beta_sd <= 0) stop('beta_sd parameter must be strictly positive.')
    # Sample
    coded_doses <- crm_codified_dose_logistic(skeleton, a0, exp(beta_mean))
    beta_sample <- rnorm(n = prior_samples, mean = beta_mean, sd = beta_sd)
    prior_curves <- t(
      sapply(beta_sample,
             function(beta) gtools::inv.logit(a0 + exp(beta) * coded_doses))
    )
    samples <- cbind(beta_sample, prior_curves)
    colnames(samples) <- c('beta', paste0('prob_tox[', dose_indices, ']'))
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
    # Sample
    coded_doses <- crm_codified_dose_logistic(skeleton, a0,
                                              beta_shape / beta_inverse_scale)
    beta_sample <- rgamma(n = prior_samples, shape = beta_shape,
                          rate = beta_inverse_scale)
    prior_curves <- t(
      sapply(beta_sample,
             function(beta) gtools::inv.logit(a0 + beta * coded_doses))
    )
    samples <- cbind(beta_sample, prior_curves)
    colnames(samples) <- c('beta', paste0('prob_tox[', dose_indices, ']'))
  } else if(model == 'logistic2') {
    # Check parameters
    if(is.null(alpha_mean)) stop('alpha_mean parameter must be specified.')
    if(is.null(alpha_sd)) stop('alpha_sd parameter must be specified.')
    if(alpha_sd <= 0) stop('alpha_sd parameter must be strictly positive.')
    if(is.null(beta_mean)) stop('beta_mean parameter must be specified.')
    if(is.null(beta_sd)) stop('beta_sd parameter must be specified.')
    if(beta_sd <= 0) stop('beta_sd parameter must be strictly positive.')
    # Sample
    coded_doses <- crm_codified_dose_logistic(skeleton, alpha_mean, beta_mean)
    alpha_sample <- rnorm(n = prior_samples, mean = alpha_mean, sd = beta_sd)
    beta_sample <- rnorm(n = prior_samples, mean = beta_mean, sd = beta_sd)
    prior_curves <- t(
      apply(data.frame(alpha_sample, beta_sample), 1,
            function(x) gtools::inv.logit(x[1] + exp(x[2]) * coded_doses))
    )
    samples <- cbind(alpha_sample, beta_sample, prior_curves)
    colnames(samples) <- c('alpha', 'beta',
                           paste0('prob_tox[', dose_indices, ']'))
  } else {
    stop(paste0(
      "Model type '", model, "' is not recognised. See ?crm_prior_beliefs"))
  }

  prob_tox <- colMeans(prior_curves)
  median_prob_tox = apply(prior_curves, 2, median)
  prior_mtd <- apply(prior_curves, 1,
                     function(x) select_dose_func(prob_tox = x, target = target))
  prob_mtd <- sapply(dose_indices, function(i) mean(prior_mtd == i))
  recommended_dose <- select_dose_func(prob_tox = skeleton, target = target)

  xyz <- crm_fit(dose_indices = dose_indices,
                 prob_tox = prob_tox,
                 median_prob_tox = median_prob_tox,
                 prob_mtd = prob_mtd,
                 recommended_dose = recommended_dose,
                 # model_dose = which.min(abs(prob_tox - target)),
                 # modal_mtd_candidate = which.max(prob_mtd),
                 dat = dat,
                 fit = NULL,
                 samples = samples)
  xyz
}
