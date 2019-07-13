
#' Get the number of efficacy events seen at the doses under investigation.
#'
#' @param x An R object of class \code{"dose_finding_fit"}
#' @param dose Optional integer, at which dose-level? Omit to get data on all doses.
#' @param ... arguments passed to other methods
#'
#' @return integer vector
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # EffTox example
#' x <- stan_efftox_demo(outcome_str = '1N 2E')
#' eff_at_dose(fit)            # c(0, 1, 0, 0)
#' eff_at_dose(fit, dose = 2)  # 1
#' eff_at_dose(fit, dose = 3)  # 0
#' }
eff_at_dose <- function(x, dose, ...) {
  UseMethod('eff_at_dose')
}

#' @rdname eff_at_dose
#' @export
eff_at_dose.efftox_fit <- function(x, dose = NULL, ...) {
  if(is.null(dose))
    sapply(x$dose_indices, function(i) sum(x$eff[x$doses == i]))
  else
    sum(x$eff[x$doses == dose])
}
