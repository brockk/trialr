
#' Class used by \pkg{trialr} to fit Wason & Seaman's Augmented Binary method in
#' single arm trials with two post-baseline tumour assessments.
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
#' \code{\link{augbin_fit}}
#' \code{\link{stan_augbin}}
augbin_2t_1a_fit <- function(num_patients,
                             tumour_size,
                             non_shrinkage_failure,
                             fit) {
  x <- augbin_fit(
    num_patients = num_patients,
    tumour_size = tumour_size,
    non_shrinkage_failure = non_shrinkage_failure,
    fit = fit)

  class(x) <- c("augbin_2t_1a_fit", class(x))
  x
}

#' Cast \code{augbin_2t_1a_fit} object to \code{\link[tibble]{tibble}}.
#'
#' @param x Object of class \code{augbin_2t_1a_fit}.
#' @param ... Extra args passed onwards.
#' @return Object of class \code{\link[tibble]{tibble}}
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#' @export
as_tibble.augbin_2t_1a_fit <- function(fit, ...) {
  data.frame(fit$tumour_size, fit$non_shrinkage_failure) %>%
    as_tibble %>%
    mutate(y1 = log(.data$z1 / .data$z0),
           y2 = log(.data$z2 / .data$z0))
}
