
#' Get the total weight of patient outcomes at the doses under investigation.
#'
#' @param x An R object of class \code{"dose_finding_fit"}
#' @param dose Optional integer, at which dose-level? Omit to get data on all doses.
#' @param ... arguments passed to other methods
#'
#' @return numerical vector
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # CRM example
#' fit <- stan_crm(skeleton = c(0.1, 0.2, 0.35, 0.6), target = 0.2,
#'                 model = 'empiric', beta_sd = sqrt(1.34), seed = 123,
#'                 doses = c(1, 1, 2, 2, 2),
#'                 tox   = c(0, 0, 0, 0, 0),
#'                 weights = c(1, 1, 0.9, 0.1, 0.1))
#'
#' total_weight_at_dose(fit)            # c(2, 1.1, 0, 0)
#' total_weight_at_dose(fit, dose = 2)  # 1.1
#' }
total_weight_at_dose <- function(x, dose, ...) {
  UseMethod("total_weight_at_dose")
}

#' @rdname total_weight_at_dose
#' @export
total_weight_at_dose.default <- function(x, dose = NULL, ...) {
  if(is.null(dose)) {
    weights <- weights_at_dose(x, dose = dose)
    map_dbl(weights, sum)
  } else {
    sum(x$dat$weights[x$doses == dose])
  }
}
