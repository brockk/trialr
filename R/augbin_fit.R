#' Class used by \pkg{trialr} to fit Wason & Seaman's Augmented Binary method.
#'
#' @param num_patients Integer, the number of patients analysed.
#' @param tumour_size matrix-like object containing tumour size measures, with
#' rows representing patients and columns representing chronological
#' standardised assessment points. Column one is baseline.
#' @param non_shrinkage_failure matrix-like object containing logical indicators
#' of non-shrinkage failure, with rows representing patients and columns
#' representing chronological standardised assessment points.
#' @param fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
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
#' @seealso
#' \code{\link{stan_augbin}}
augbin_fit <- function(num_patients,
                       tumour_size,
                       non_shrinkage_failure,
                       fit) {

  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )

  x <- list(num_patients = num_patients,
            tumour_size = tumour_size,
            non_shrinkage_failure = non_shrinkage_failure,
            fit = fit,
            version = version)

  class(x) <- c("augbin_fit", "list")
  x
}

#' Print augbin_fit object.
#'
#' @param x \code{\link{augbin_fit}} object to print.
#' @param pars parameters in model to summarise.
#' @param ... Extra parameters, passed onwards.
#' @method print augbin_fit
#' @export
print.augbin_fit <- function(x,
                             pars = c('alpha', 'beta', 'gamma', 'Omega',
                                      'sigma', 'alphaD1', 'gammaD1', 'alphaD2',
                                      'gammaD2'),
                             ...) {
  print(x$fit, pars = pars, ...)
}
