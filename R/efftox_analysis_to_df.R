
#' @title EffTox analysis to data.frame
#'
#' @description Convenient function to turn an \code{\link{efftox_fit}}
#' into a \code{data.frame}.
#'
#' @param x An instance of \code{\link{efftox_fit}}
#'
#' @return a \code{data.frame}
#'
#' @export
#'
#' @examples
#' fit <- stan_efftox_demo(outcome_str = '1N 2E 3B')
#' df <- efftox_analysis_to_df(fit)
#' df
#'
#' @seealso
#' \code{\link{stan_efftox}}
efftox_analysis_to_df <- function(x) {
  df <- data.frame(
    Dose = factor(x$dose_indices),
    N = sapply(1:length(x$dose_indices), function(i) sum(x$doses == i)),
    ProbEff = x$prob_eff,
    ProbTox = x$prob_tox,
    ProbAccEff = x$prob_acc_eff,
    ProbAccTox = x$prob_acc_tox,
    Utility = x$utility, Acceptable = x$acceptable,
    ProbOBD = x$prob_obd
  )

  return(df)
}
