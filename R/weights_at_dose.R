
#' Get the weights of patient outcomes at the doses under investigation.
#'
#' @param x An R object of class \code{"dose_finding_fit"}
#' @param dose Optional integer, at which dose-level? Omit to get data on all doses.
#' @param ... arguments passed to other methods
#'
#' @return list if \code{dose} omitted, numerical vector if \code{dose} provided.
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
#' l <- weights_at_dose(fit)
#'
#' length(l)  # 4
#' l[[1]]  # c(1, 1)
#' l[[2]]  # c(0.9, 0.1, 0.1)
#' l[[3]]  # c()
#'
#' weights_at_dose(fit, dose = 2)  # c(0.9, 0.1, 0.1)
#' }
weights_at_dose <- function(x, dose, ...) {
  UseMethod("weights_at_dose")
}

#' @rdname weights_at_dose
#' @export
weights_at_dose.default <- function(x, dose = NULL, ...) {
  if(is.null(dose)) {
    nn <- n_at_dose(x, dose = dose, ...)
    lapply(nn, function(i) rep(1, i))
  } else {
    rep(1, n_at_dose(x, dose = dose))
  }
}

#' @rdname weights_at_dose
#' @export
weights_at_dose.crm_fit <- function(x, dose = NULL, ...) {
  if(is.null(dose))
    lapply(x$dose_indices, function(i) x$weights[x$doses == i])
  else
    x$weights[x$doses == dose]
}
