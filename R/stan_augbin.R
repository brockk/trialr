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
#' size measurements. See References.
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
#' versions will be added in future releases. The model choice determiens the
#' prior parameters that must be provided. See sections below.
#' @param prior_params list of prior parameters. These are combined with the
#' data and passed to \code{rstan::sampling}. The parameters required depend on
#' the model form being fit. See sections below.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. See \code{\link[rstan:sampling]{sampling}}.
#'
#' @section Parameters for single-arm model with two post-baseline assessments
#' \itemize{
#'   \item \code{alpha_mean} & \code{alpha_sd} for normal prior on alpha.
#'   \item \code{beta_mean} & \code{beta_sd} for normal prior on beta.
#'   \item \code{gamma_mean} & \code{gamma_sd} for normal prior on gamma.
#'   \item \code{sigma_mean} & \code{sigma_sd} for normal priors on diagonal
#'   elements of Sigma, the two-period covariance matrix of the log-tumour-size
#'   ratios.
#'   \item \code{omega_lkj_eta} for a LKJ prior on the two-period correlation
#'   matrix of the log-tumour-size ratios. omega_lkj_eta = 1 is uniform,
#'   analogous to a Beta(1,1) prior on a binary probability.
#'   \item \code{alpha_d1_mean} & \code{alpha_d1_sd} for normal prior on alpha_D1.
#'   \item \code{gamma_d1_mean} & \code{gamma_d1_sd} for normal prior on gamma_D1.
#'   \item \code{alpha_d2_mean} & \code{alpha_d2_sd} for normal prior on alpha_D2.
#'   \item \code{gamma_d2_mean} & \code{gamma_d2_sd} for normal prior on gamma_D2.
#' }
#'
#' @return an instance or subclass of type \code{\link{augbin_fit}}.
#'
#' @author Kristian Brock
#'
#' @seealso
#'   \code{\link{crm_fit}}
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
#' \dontrun{
#'   # TODO
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
