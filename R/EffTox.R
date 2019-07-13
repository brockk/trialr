
# Generics ----
#' Print efftox_fit object.
#'
#' @param x \code{\link{efftox_fit}} object to convert.
#' @param ... Extra parameters, passed onwards.
#' @method print efftox_fit
#' @export
print.efftox_fit <- function(x, ...) {
  # Patient-level data
  if(x$num_patients > 0) {
    treated <- data.frame(
      Patient = 1:length(x$doses),
      Dose = x$doses,
      Toxicity = x$tox,
      Efficacy = x$eff
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

