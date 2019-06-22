
#' Class used by \pkg{trialr} to fit Wason & Seaman's Augmented Binary method in
#' single arm trials with two post-baseline tumour assessments.
#'
#' @param num_patients Integer, the number of patients analysed.
#' @param tumour_size matrix-like object containing tumour size measures, with
#' rows representing patients and columns representing chronological
#' assessment points. Column one is baseline.
#' @param non_shrinkage_failure matrix-like object containing logical indicators
#' of non-shrinkage failure, with rows representing patients and columns
#' representing chronological assessment points.
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
#' @importFrom rlang .data
#' @export
as_tibble.augbin_2t_1a_fit <- function(x, ...) {
  data.frame(x$tumour_size, x$non_shrinkage_failure) %>%
    as_tibble %>%
    mutate(y1 = log(.data$z1 / .data$z0),
           y2 = log(.data$z2 / .data$z0))
}

#' Predict probability of success for given tumour size measurements.
#'
#' This method simply forwards to \code{\link{prob_success}}.
#'
#' @param object Object of class \code{augbin_2t_1a_fit}.
#' @param y1_lower numeric, minimum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 1 to baseline.
#' Defaults to negative infinity.
#' @param y1_upper numeric, maximum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 1 to baseline.
#' Defaults to positive infinity.
#' @param y2_lower numeric, minimum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 2 to baseline.
#' @param y2_upper numeric, maximum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 2 to baseline.
#' Defaults to log(0.7).
#' @param probs pair of probabilities to use to calculate the credible interval
#' for the probability of success.
#' @param newdata data for which to infer the probability of success.
#' A dataframe-like object with baseline tumour sizes in first column, and first
#' and second post-baseline tumour sizes in columns 2 and 3. Omitted by default.
#' When omitted, newdata is set to be the \code{object$tumour_size}.
#' @param ... Extra args passed onwards.
#'
#' @return Object of class \code{\link[tibble]{tibble}}
#'
#' @export
predict.augbin_2t_1a_fit <- function(object,
                                     y1_lower = -Inf, y1_upper = Inf,
                                     y2_lower = -Inf, y2_upper = log(0.7),
                                     probs = c(0.025, 0.975),
                                     newdata = NULL,
                                     ...) {
  prob_success(object, y1_lower = y1_lower, y1_upper = y1_upper,
               y2_lower = y2_lower, y2_upper = y2_upper,
               probs = probs, newdata = newdata, ...)
}
