#' @title Calculate the p-index for EffTox utility contours
#'
#' @description
#' Calculate the p-index for EffTox utility contours so that the neutral utility
#' contour intersects the following points in the
#' Prob(Efficacy) - Prob(Toxicity) plane:
#' (\code{eff0}, 0), (1, \code{tox1}) and (\code{eff_star}, \code{tox_star})
#'
#' @param eff0 Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1
#' @param tox1 Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1
#' @param eff_star Efficacy probability of an equi-utility third point
#' @param tox_star Toxicity probability of an equi-utility third point
#'
#' @return The p-index
#'
#' @export
#'
#' @examples
#' efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' @references Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample
#' size for computing prior hyperparameters in Bayesian phase I-II dose-finding
efftox_solve_p <- function(eff0, tox1, eff_star, tox_star) {
  # Calculate p for the efficacy/toxicity contours that will intersect points
  # (eff0, 0), (eff.star, tox.star), and (1, tox1)

  if(any(eff0 <= 0, tox1 <= 0, eff_star <= 0, tox_star <= 0,
         eff0 >= 1, tox1 >= 1, eff_star >= 1, tox_star >= 1)) {
    stop('eff0, tox1, eff_star and tox_star must all be between 0 and 1.')
  }

  .objective = function(p, eff0, tox1, eff_star, tox_star) {
    a <- ((1 - eff_star) / (1 - eff0))
    b <- tox_star / tox1
    return(a^p + b^p - 1)
  }

  rt <- stats::uniroot(.objective, interval = c(0, 100),
                       eff0 = eff0, tox1 = tox1, eff_star = eff_star,
                       tox_star = tox_star)
  return(rt$root)
}
