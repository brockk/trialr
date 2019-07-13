
#' @title Calculate codified CRM doses.
#'
#' @description Calculate the codified CRM doses that map to probability of
#' toxicity \code{prob_tox} in a logistic model with expected values for
#' intercept and gradient. I.e. find \eqn{x[i]} such that
#' \eqn{logit(p[i]) = \alpha + \beta x[i]}, were \eqn{p} is
#' \code{prob_tox}.
#'
#' @param prob_tox Numeric vector, seek codified doses that yield these
#' probabilities of toxicity.
#' @param alpha_mean Numeric, expected value of intercept.
#' @param beta_mean Numeric, expected value of gradient with respect to dose.
#'
#' @return Numeric vector of codified doses.
#' @export
#'
#' @examples
#' skeleton <- c(0.05, 0.1, 0.2, 0.5)
#' crm_codified_dose_logistic(skeleton, 1, 0)
#' crm_codified_dose_logistic(skeleton, 3, 0.5)
crm_codified_dose_logistic <- function(prob_tox, alpha_mean, beta_mean) {
  (gtools::logit(prob_tox) - alpha_mean) / beta_mean
}
