
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
#' @section Requirements of \code{empiric} model:
#' * beta_sd
#'
#' @section Requirements of \code{logistic} model:
#' * a0
#' * beta_mean
#' * beta_sd
#'
#' @section Requirements of \code{logistic_gamma} model:
#' * a0
#' * beta_shape
#' * beta_inverse_scale
#'
#' @section Requirements of \code{logistics} model:
#' * a0
#' * alpha_mean
#' * alpha_sd
#' * beta_mean
#' * beta_sd
#'
#' @seealso
#' \code{\link{stan_crm}}
#' \code{\link{crm_process}}
crm_params <- function(skeleton, target, a0 = NULL,
                       alpha_mean = NULL, alpha_sd = NULL,
                       beta_mean = NULL, beta_sd = NULL,
                       beta_shape = NULL, beta_inverse_scale = NULL) {

  # crm_params class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )

  x <- loo::nlist(num_doses = length(skeleton), skeleton, target,
                  a0, alpha_mean, alpha_sd,
                  beta_mean, beta_sd, beta_shape, beta_inverse_scale)

  # Initialise with no patients observed
  x$doses = c()
  x$tox = c()
  x$num_patients = 0

  # Set type. This is, at heart, just a list.
  class(x) <- c("crm_params", "list")
  return(x)
}


#' Fit a CRM model
#'
#' Fit a CRM model using Stan for full Bayesian inference.
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
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. \code{\link[rstan:sampling]{sampling}}.
#'
#' @details
#' The quickest and easiest way to fit a CRM model to some observed outcomes
#' is to describe the outcomes using \pkg{trialr}'s syntax for dose-finding
#' outcomes. See \code{\link{df_parse_outcomes}} for full details and examples.
#'
#' Different model parameterisations require that difference parameter values
#' are specified.
#'
#' @section Requirements of \code{empiric} model:
#' * beta_sd
#'
#' @section Requirements of \code{logistic} model:
#' * a0
#' * beta_mean
#' * beta_sd
#'
#' @section Requirements of \code{logistic_gamma} model:
#' * a0
#' * beta_shape
#' * beta_inverse_scale
#'
#' @section Requirements of \code{logistics} model:
#' * a0
#' * alpha_mean
#' * alpha_sd
#' * beta_mean
#' * beta_sd
#'
#' @return An object of class \code{\link{crm_fit}}
#'
#' @author Kristian Brock \email{kristian.brock@@gmail.com}
#'
#' @references
#'   O'Quigley, J., Pepe, M., & Fisher, L. (1990).
#'   Continual reassessment method: a practical design for phase 1 clinical
#'   trials in cancer.
#'   Biometrics, 46(1), 33-48. https://doi.org/10.2307/2531628
#'
#' @seealso
#'   \code{\link{crm_fit}}
#'   \code{\link{crm_process}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This model is presented in Thall et al. (2014)
#' mod1 <- stan_crm('1N 2N 3T', skeleton = c(0.1, 0.2, 0.35, 0.6),
#'                  target = 0.2, model = 'empiric', beta_sd = sqrt(1.34),
#'                  seed = 123)
#'
#' # Shorthand for the above is:
#' mod2 <- stan_efftox_demo('1N 2E 3B', seed = 123)
#'
#' # the seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
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
                    beta_inverse_scale = beta_inverse_scale
  )

  # Add outcomes
  if(is.null(outcome_str)) {
    if(length(doses_given) != length(tox))
      stop('doses_given and tox vectors should have same length')
    dat$doses <- doses_given
    dat$tox <- tox
    dat$num_patients <- length(doses_given)
  } else {
    outcomes_df <- df_parse_outcomes(outcome_str, as.list = TRUE)
    dat$num_patients <- outcomes_df$num_patients
    dat$doses <- outcomes_df$doses
    dat$tox <- outcomes_df$tox
  }

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

#' Class of model fit by \pkg{trialr} using the CRM dose-finding design.
#'
#' @name crm_fit-class
#' @aliases crm_fit
#' @docType class
#'
#' @details
#' See \code{methods(class = "crm_fit")} for an overview of available
#' methods.
#'
#' @param dose_indices A vector of integers representing the dose-levels under
#' consideration.
#' @param recommended_dose An integer representing the dose-level recommended
#' for the next patient or cohort; or \code{NA} if stopping is recommended.
#' The recommended dose typically has associated probability of DLT closest to
#' the target toxicity rate. Contrast to \code{modal_mtd_candidate}.
#' @param prob_tox The posterior mean probabilities of toxicity at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param median_prob_tox The posterior median probabilities of toxicity at doses
#' 1:n; a vector of numbers between 0 and 1.
#' @param modal_mtd_candidate An integer representing the dose-level most likely
#' to be the MTD, i.e. the dose-level that maximises \code{prob_mtd}.
#' @param prob_mtd The posterior probability that each dose is the MTD, by the
#' chosen model; a vector of numbers between 0 and 1.
#' @param dat Object \code{\link{crm_params}} containing data passed to
#' \code{\link[rstan:sampling]{sampling}}.
#' @param fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
#'
#' @seealso
#' \code{\link{stan_crm}}
#' \code{\link{crm_process}}
crm_fit <- function(dose_indices, recommended_dose, prob_tox, median_prob_tox,
                    modal_mtd_candidate, prob_mtd, dat, fit) {
  # crm_fit class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )
  x <- loo::nlist(dose_indices, recommended_dose, prob_tox, median_prob_tox,
                  modal_mtd_candidate, prob_mtd, dat, fit, version)
  class(x) <- "crm_fit"
  x
}

#' @title Process RStan samples from a CRM model
#'
#' @description Process RStan samples from a CRM model to make inferences
#' about dose-toxicity and which dose should be recommended next.
#' Typically, this function is not required to be called explicitly by the user
#' because \code{\link{stan_crm}} will call it implicitly.
#'
#' @param dat An instance of \code{\link{crm_params}}, a list of CRM
#' parameters.
#' @param fit An instance of \code{rstan::stanmodel}, derived by fitting one of
#' the trialr CRM models.
#' @return An instance of \code{\link{crm_fit}}.
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- list(
#'   num_doses = 5,
#'   skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55),
#'   target = 0.25,
#'   beta_sd = sqrt(1.34),
#'   num_patients = 3,
#'   doses = c(1, 2, 3),
#'   tox = c(0, 0, 1)
#' )
#' samp <- rstan::sampling(stanmodels$CrmEmpiricNormalPrior,
#'                         data = dat, seed = 123)
#' decision <- crm_process(dat, samp)
#' }
#'
#' @seealso
#' \code{\link{stan_crm}}
#' \code{\link{crm_params}}
crm_process <- function(dat, fit) {
  dose_indices <- seq(from = 1, to = dat$num_doses, by = 1)
  # Posterior estimates
  prob_tox_samp <- rstan::extract(fit, 'prob_tox')[[1]]
  prob_tox <- colMeans(prob_tox_samp)
  recommended_dose <- which.min(abs(prob_tox - dat$target))
  median_prob_tox <- apply(prob_tox_samp, 2, stats::median)
  # Implied MTD
  implied_mtd <- apply(prob_tox_samp, 1, function(x) which.min(abs(x - dat$target)))
  prob_mtd <- sapply(dose_indices, function(x) mean(implied_mtd == x))
  modal_mtd_candidate <- which.max(prob_mtd)

  x <- crm_fit(dose_indices, recommended_dose, prob_tox, median_prob_tox,
               modal_mtd_candidate, prob_mtd, dat, fit)
  return(x)
}


# Generics ----
#' Print crm_fit object.
#'
#' @param x \code{\link{crm_fit}} object to convert.
#' @param ... Extra parameters, passed onwards.
#' @sdname print
#' @method print crm_fit
#' @S3method print crm_fit
print.crm_fit <- function(x, ...) {
  # Patient-level data
  treated <- data.frame(
    Patient = 1:length(x$dat$doses),
    Dose = x$dat$doses,
    Toxicity = x$dat$tox
  )
  print(treated)
  cat('\n')

  # Dose-level data
  df <- data.frame(
    DoseLevel = factor(x$dose_indices),
    Skeleton = x$dat$skeleton,
    N = sapply(1:x$dat$num_doses, function(i) sum(x$dat$doses == i)),
    Tox = sapply(1:x$dat$num_doses, function(i) sum(x$dat$tox[x$dat$doses == i])),
    ProbTox = x$prob_tox,
    ProbMTD = x$prob_mtd
  )
  print(df)
  cat('\n')

  # Extras
  cat(paste0('The model targets a toxicity level of ',
             x$dat$target, '.'))
  cat('\n')
  cat(paste0('The dose with estimated toxicity probability closest to target is ',
               x$recommended_dose, '.'))
  cat('\n')
  cat(paste0('The dose most likely to be the MTD is ',
             x$modal_mtd_candidate, '.'))
}

#' Convert crm_fit object to \code{data.frame}.
#'
#' @param x \code{\link{crm_fit}} object to convert.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A \code{data.frame}
#' @sdname as.data.frame
#' @method as.data.frame crm_fit
#' @S3method as.data.frame crm_fit
as.data.frame.crm_fit <- function(x, ...) {
  as.data.frame(x$fit, ...)
}

#' Plot an crm_fit
#'
#' @param x \code{\link{crm_fit}} object to plot.
#' @param pars Parameters to plot. Plots utility scores by default.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A plot
#' @sdname plot
#' @method plot crm_fit
#' @S3method plot crm_fit
plot.crm_fit <- function(x, pars = 'prob_tox', ...) {
  rstan::plot(x$fit, pars = pars, ...)
}

#' Obtain summary of an crm_fit
#'
#' @param object \code{\link{crm_fit}} object to summarise.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A summary object.
#' @sdname summary
#' @method summary crm_fit
#' @S3method summary crm_fit
summary.crm_fit <- function(object, ...) {
  rstan::summary(object$fit, ...)
}


# Not generic yet.... ----
# tidybayes is not yet on CRAN but once it is, add it as an import and
# implement as_sample_tibble.crm_fit(x).
# However, for now:

#' Extract tall data.frame of posterior prob_tox samples.
#'
#' @param x \code{\link{crm_fit}} object.
#'
#' @return data.frame
#' @export
gather_samples.crm_fit <- function(x) {
  df <- as.data.frame(x, 'prob_tox')
  Label <- ProbTox <- NULL
  df_tall <- df %>%
    tidyr::gather(Label, ProbTox) %>%
    dplyr::mutate(
      DoseLevel = rep(1:ncol(df), each = nrow(df)),
      Draw = rep(1:nrow(df), times = ncol(df))
    )
  df_tall
}
