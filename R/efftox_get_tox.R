
#' @title Get the Prob(Tox) for Prob(Eff) and utility pairs
#'
#' @description Get the probability of toxicity for probability-of-efficacy and
#' utility pairs
#'
#' @param eff Probability of efficacy; number between 0 and 1
#' @param util Utility score; number
#' @param p p-index of EffTox utility contours. Use \code{efftox_solve_p}
#' @param eff0 Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1
#' @param tox1 Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1
#'
#' @return Probability(s) of toxicity
#'
#' @export
#'
#' @examples
#' p <- efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' prob_tox <- efftox_get_tox(0.7, 0, p, eff0 = 0.5, tox1 = 0.65)
#' round(prob_tox, 2) == 0.25
#'
#' prob_tox <- efftox_get_tox(0.7, seq(-0.5, 0.25, by = 0.25), p, eff0 = 0.5,
#'                            tox1 = 0.65)
#' round(prob_tox, 2) == c(0.57, 0.41, 0.25, 0.09)
#'
#' prob_tox <- efftox_get_tox(c(0.5, 0.7, 0.8), 0.25, p, eff0 = 0.5, tox1 = 0.65)
#' round(prob_tox, 2) == c(NaN, 0.09, 0.22)
#'
#' prob_tox <- efftox_get_tox(c(0.5, 0.7, 0.8), c(-1, 0, 1), p, eff0 = 0.5,
#'                            tox1 = 0.65)
#' round(prob_tox, 2) == c(0.63, 0.25, NaN)
#'
#' @note Various ways of vectorising the function are demonstrated in the
#' examples
#'
#' @seealso \code{\link{efftox_solve_p}}
efftox_get_tox <- function(eff, util, p, eff0, tox1) {

  a = ((1 - eff) / (1 - eff0))
  return(tox1 * ((1 - util)^p - a^p)^(1 / p))
}
