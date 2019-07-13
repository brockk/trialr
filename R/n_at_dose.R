
#' Get the number of patients treated at the doses under investigation.
#'
#' @param x An R object of class \code{"dose_finding_fit"}
#' @param dose Optional integer, at which dose-level? Omit to get data on all doses.
#' @param ... arguments passed to other methods
#'
#' @return integer vector
#' @export
#'
#' @examples
#' \dontrun{
#' # CRM example
#' target <- 0.2
#' fit <- stan_crm('1N 2N 3T', skeleton = c(0.1, 0.2, 0.35, 0.6),
#'                  target = target, model = 'empiric', beta_sd = sqrt(1.34),
#'                  seed = 123)
#' n_at_dose(fit)            # c(1, 1, 1, 0)
#' n_at_dose(fit, dose = 3)  # 1
#' }
n_at_dose <- function(x, dose, ...) {
  UseMethod('n_at_dose')
}

#' @rdname n_at_dose
#' @export
n_at_dose.dose_finding_fit <- function(x, dose = NULL, ...) {
  if(is.null(dose))
    sapply(x$dose_indices, function(i) sum(x$doses == i))
  else
    sum(x$doses == dose)
}
