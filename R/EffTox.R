
# Simple helpers ----
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
    N = sapply(1:x$dat$num_doses, function(i) sum(x$dat$doses == i)),
    ProbEff = x$prob_eff,
    ProbTox = x$prob_tox,
    ProbAccEff = x$prob_acc_eff,
    ProbAccTox = x$prob_acc_tox,
    Utility = x$utility, Acceptable = x$acceptable,
    ProbOBD = x$prob_obd
  )

  return(df)
}

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




# Generics ----
#' Print efftox_fit object.
#'
#' @param x \code{\link{efftox_fit}} object to convert.
#' @param ... Extra parameters, passed onwards.
#' @method print efftox_fit
#' @export
print.efftox_fit <- function(x, ...) {
  # Patient-level data
  if(x$dat$num_patients > 0) {
    treated <- data.frame(
      Patient = 1:length(x$dat$doses),
      Dose = x$dat$doses,
      Toxicity = x$dat$tox,
      Efficacy = x$dat$eff
    )
    print(treated)
  } else {
    cat('No patients have been treated.\n')
  }
  cat('\n')

  # Dose-level data
  df <- efftox_analysis_to_df(x)
  print(df, digits = 3)
  cat('\n')

  # Extras
  if(x$num_patients > 0 & sum(x$acceptable) == 0) {
    cat('The model advocates stopping.')
  } else {
    if(!is.na(x$recommended_dose)) {
      cat(paste0('The model recommends selecting dose-level ',
                 x$recommended_dose, '.'))
      cat('\n')
    }

    cat(paste0('The dose most likely to be the OBD is ',
               x$modal_obd_candidate, '.'))
    cat('\n')
    cat(paste0('Model entropy: ', format(round(x$entropy, 2), nsmall = 2)))
  }
}

#' Convert efftox_fit object to \code{data.frame}.
#'
#' @param x \code{\link{efftox_fit}} object to convert.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A \code{data.frame}
#' @method as.data.frame efftox_fit
#' @export
as.data.frame.efftox_fit <- function(x, ...) {
  as.data.frame(x$fit, ...)
}

#' Plot an efftox_fit
#'
#' @param x \code{\link{efftox_fit}} object to plot.
#' @param pars Parameters to plot. Plots utility scores by default.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A plot
#' @method plot efftox_fit
#' @export
plot.efftox_fit <- function(x,  pars = 'utility', ...) {
  rstan::plot(x$fit, pars = pars, ...)
}

#' Obtain summary of an efftox_fit
#'
#' @param object \code{\link{efftox_fit}} object to summarise.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A summary object.
#' @method summary efftox_fit
#' @export
summary.efftox_fit <- function(object, ...) {
  rstan::summary(object$fit, ...)
}

#' @title Convert \code{\link{efftox_fit}} to instance of
#' \code{\link[coda]{mcmc.list}}
#'
#' @description This function allows trialr to use tidybayes functions.
#'
#' @param efftox_fit Object of class \code{\link{efftox_fit}}
#' @param ... Extra variables that are passed onwards.
#'
#' @return Object of class \code{\link[coda]{mcmc.list}}
#' @method as.mcmc.list efftox_fit
#'
#' @importFrom coda as.mcmc.list
#' @importFrom rstan As.mcmc.list
#' @export
as.mcmc.list.efftox_fit <- function(efftox_fit, ...) {
  As.mcmc.list(efftox_fit$fit, ...)
}

