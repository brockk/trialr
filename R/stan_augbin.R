#' Fit Wason & Seaman's Augmented Binary model for tumour response.
#'
#' Phase II clinical trials in oncology commonly assess response as a key outcome
#' measure. Patients achieve a RECIST response if their tumour size post-baseline
#' has changed in size by some threshold amount and they do not experience
#' non-shrinkage failure. An example of non-shrinkage failure is the appearance
#' of new lesions. As a dichtotomisation of the underlying continuous tumour size
#' measurement, RECIST response is inefficient. Wason & Seaman introduced the
#' Augmented Binary method to incorporate mechanisms for non-shrinkage failure
#' whilst modelling the probability of response based on the continuous tumour
#' size measurements. See model-specific sections below, and the references.
#'
#' @param tumour_size matrix-like object containing tumour size measures, with
#' rows representing patients and columns representing chronological
#' standardised assessment points. Column one is baseline.
#' @param non_shrinkage_failure matrix-like object containing logical indicators
#' of non-shrinkage failure, with rows representing patients and columns
#' representing chronological standardised assessment points.
#' @param arm optional vector of integers representing the allocated treatment
#' arms for patients, assumed in the same order as \code{tumour_size} and
#' \code{non_shrinkage_failure}. NULL to fit the augbin variant for single-arm
#' trials. NULL is the default.
#' @param model Character string to denote the desired model. Currently, only
#' \code{2t-1a} is supported, representing the model variant with two
#' post-baseline assessments in a single arm trial. Multi-period and multi-arm
#' versions will be added in future releases. The model choice determines the
#' prior parameters that must be provided. See sections below.
#' @param prior_params list of prior parameters. These are combined with the
#' data and passed to \code{rstan::sampling}. The parameters required depend on
#' the model form being fit. See sections below.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. See \code{\link[rstan:sampling]{sampling}}.
#'
#' @section Single-arm model with two post-baseline assessments:
#' The complete model form is:
#' \deqn{(y_{1i}, y_{2i})^T \sim N( (\mu_{1i}, \mu_{2i})^T, \Sigma) }
#' \deqn{ \mu_{1i} = \alpha + \gamma z_{0i} }
#' \deqn{ \mu_{2i} = \beta + \gamma z_{0i} }
#' \deqn{ logit(Pr(D_{1i} = 1 | Z_{0i})) = \alpha_{D1} + \gamma_{D1} z_{0i} }
#' \deqn{ logit(Pr(D_{2i} = 1 | D_{1i} = 0, Z_{0i}, Z_{1i})) = \alpha_{D2} + \gamma_{D2} z_{1i} }
#' where \eqn{z_{0i}, z_{1i}, z_{2i}} are tumour sizes at baseline, period 1,
#' and period 2, for patient i; \eqn{y_{1i}, y_{2i}} are the log-tumour-size
#' ratios with respect to baseline; \eqn{D_{1i}, D_{2i}} are indicators of
#' non-shrinkage failure; and \eqn{\Sigma} is assumed to be unstructured
#' covariance matrix, with associated correlation matrix having an LKJ prior.
#'
#' The following prior parameters are required:
#' \itemize{
#'   \item \code{alpha_mean} & \code{alpha_sd} for normal prior on \eqn{\alpha}.
#'   \item \code{beta_mean} & \code{beta_sd} for normal prior on \eqn{\beta}.
#'   \item \code{gamma_mean} & \code{gamma_sd} for normal prior on \eqn{\gamma}.
#'   \item \code{sigma_mean} & \code{sigma_sd} for normal priors on diagonal elements of \eqn{\Sigma};
#'   \item \code{omega_lkj_eta} for a LKJ prior on the two-period correlation matrix associated with Sigma. omega_lkj_eta = 1 is uniform, analogous to a Beta(1,1) prior on a binary probability.
#'   \item \code{alpha_d1_mean} & \code{alpha_d1_sd} for normal prior on \eqn{\alpha_{D1}}.
#'   \item \code{gamma_d1_mean} & \code{gamma_d1_sd} for normal prior on \eqn{\gamma_{D1}}.
#'   \item \code{alpha_d2_mean} & \code{alpha_d2_sd} for normal prior on \eqn{\alpha_{D2}}.
#'   \item \code{gamma_d2_mean} & \code{gamma_d2_sd} for normal prior on \eqn{\gamma_{D2}}.
#' }
#'
#' @return an instance or subclass of type \code{\link{augbin_fit}}.
#'
#' @author Kristian Brock
#'
#' @seealso
#'   \code{\link{augbin_fit}}
#'   \code{\link{prior_predictive_augbin_2t_1a}}
#'   \code{\link[rstan:sampling]{sampling}}
#'
#' @export
#'
#' @references
#' Wason JMS, Seaman SR. Using continuous data on tumour measurements to improve
#' inference in phase II cancer studies.
#' Statistics in Medicine. 2013;32(26):4639-4650. doi:10.1002/sim.5867
#'
#' Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation
#' criteria in solid tumours: Revised RECIST guideline (version 1.1).
#' European Journal of Cancer. 2009;45(2):228-247. doi:10.1016/j.ejca.2008.10.026
#'
#' @examples
#' priors <- list(alpha_mean = 0, alpha_sd = 1,
#'                beta_mean = 0, beta_sd = 1,
#'                gamma_mean = 0, gamma_sd = 1,
#'                sigma_mean = 0, sigma_sd = 1,
#'                omega_lkj_eta = 1,
#'                alpha_d1_mean = 0, alpha_d1_sd = 1,
#'                gamma_d1_mean = 0, gamma_d1_sd = 1,
#'                alpha_d2_mean = 0, alpha_d2_sd = 1,
#'                gamma_d2_mean = 0, gamma_d2_sd = 1)
#' # Scenario 1 of Table 1 in Wason & Seaman (2013)
#' N <- 50
#' sigma <- 1
#' delta1 <- -0.356
#' mu <- c(0.5 * delta1, delta1)
#' Sigma = matrix(c(0.5 * sigma^2, 0.5 * sigma^2, 0.5 * sigma^2, sigma^2),
#'                ncol = 2)
#' alphaD <- -1.5
#' gammaD <- 0
#' set.seed(123456)
#' y <- MASS::mvrnorm(n = N, mu, Sigma)
#' z0 <- runif(N, min = 5, max = 10)
#' z1 <- exp(y[, 1]) * z0
#' z2 <- exp(y[, 2]) * z0
#' d1 <- rbinom(N, size = 1, prob = gtools::inv.logit(alphaD + gammaD * z0))
#' d2 <- rbinom(N, size = 1, prob = gtools::inv.logit(alphaD + gammaD * z1))
#' tumour_size <- data.frame(z0, z1, z2) # Sizes in cm
#' non_shrinkage_failure <- data.frame(d1, d2)
#' # Fit
#' \dontrun{
#' fit <- stan_augbin(tumour_size, non_shrinkage_failure,
#'                    prior_params = priors, model = '2t-1a', seed = 123)
#' }
stan_augbin <- function(tumour_size, non_shrinkage_failure,
                        arm = NULL,
                        model = c('2t-1a'),
                        prior_params = list(),
                        ...) {

  model <- match.arg(model)
  if(model != '2t-1a') {
    msg <- paste0('Only the two-period, single-arm AugBin model is currently ',
                  'supported via "2t-1a". Multi-period and multi-arm variants ',
                  'will be added in future releases.')
    stop(msg)
  }

  num_times <- ncol(tumour_size)
  if(num_times != 3)
    stop('tumour_size should have three columns.')
  if(ncol(non_shrinkage_failure) != 2)
    stop('non_shrinkage_failure should have two columns.')
  if(nrow(tumour_size) != nrow(non_shrinkage_failure))
    stop('tumour_size and non_shrinkage_failure should have same number of rows.')
  if(!is.null(arm)) {
    stop('The multi-arm version of AugBin is not yet supported. Set arm = NULL')
    if(nrow(tumour_size) != length(arm))
      stop('The number of items in arm should match nrow(tumour_size).')
  }
  # Check priors
  if(model == '2t-1a') {
    reqd_params <- c("alpha_mean", "alpha_sd", "beta_mean", "beta_sd",
                     "gamma_mean", "gamma_sd", "sigma_mean", "sigma_sd",
                     "omega_lkj_eta", "alpha_d1_mean", "alpha_d1_sd",
                     "gamma_d1_mean", "gamma_d1_sd", "alpha_d2_mean",
                     "alpha_d2_sd", "gamma_d2_mean", "gamma_d2_sd")
    if(!all(reqd_params %in% names(prior_params))) {
      missing_params <- reqd_params[!(reqd_params %in% names(prior_params))]
      stop(paste0('\nParameter ', missing_params, ' is required in prior_params',
                  ' for model of type ', model))
    }
  }

  data <- list(N = nrow(tumour_size),
               z0 = tumour_size[, 1],
               z1 = tumour_size[, 2],
               z2 = tumour_size[, 3],
               d1 = non_shrinkage_failure[, 1],
               d2 = non_shrinkage_failure[, 2]
  )

  # Append prior parameters
  data <- append(data, prior_params)

  # Fit in Stan
  samp <- rstan::sampling(stanmodels$AugBin2T1A, data = data, ...)
  fit <- augbin_2t_1a_fit(num_patients = nrow(tumour_size),
                          tumour_size = tumour_size,
                          non_shrinkage_failure = non_shrinkage_failure,
                          fit = samp)
  fit
}
