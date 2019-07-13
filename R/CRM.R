


# Generics ----
#' Print crm_fit object.
#'
#' @param x \code{\link{crm_fit}} object to print.
#' @param ... Extra parameters, passed onwards.
#' @method print crm_fit
#' @export
print.crm_fit <- function(x, ...) {
  # Patient-level data
  if(x$num_patients > 0) {
    treated <- data.frame(
      Patient = 1:length(x$doses),
      Dose = x$doses,
      Toxicity = x$tox,
      Weight = x$weights
    )
    print(treated)
  } else {
    cat('No patients have been treated.\n')
  }
  cat('\n')

  # Dose-level data
  df <- data.frame(
    Dose = factor(x$dose_indices),
    Skeleton = x$dat$skeleton,
    N = sapply(1:length(x$dose_indices), function(i) sum(x$doses == i)),
    Tox = sapply(1:length(x$dose_indices), function(i) sum(x$tox[x$doses == i])),
    ProbTox = x$prob_tox,
    MedianProbTox = x$median_prob_tox,
    ProbMTD = x$prob_mtd
  )
  print(df, digits = 3)
  cat('\n')

  # Extras
  cat(paste0('The model targets a toxicity level of ',
             x$dat$target, '.'))
  cat('\n')
  cat(paste0('The dose with estimated toxicity probability closest to target is ',
             x$recommended_dose, '.'))
  cat('\n')
  cat(paste0('The dose most likely to be the MTD is ',
             x$modal_mtd_candidate, '.'))
  cat('\n')
  cat(paste0('Model entropy: ', format(round(x$entropy, 2), nsmall = 2)))
}

#' Convert crm_fit object to \code{data.frame}.
#'
#' @param x \code{\link{crm_fit}} object to convert.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A \code{data.frame}
#' @method as.data.frame crm_fit
#' @export
as.data.frame.crm_fit <- function(x, ...) {
  if(!is.null(x$fit))
    as.data.frame(x$fit, ...)
  else
    as.data.frame(x$samples)
}

#' Plot an crm_fit
#'
#' @param x \code{\link{crm_fit}} object to plot.
#' @param pars Parameters to plot. Plots utility scores by default.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A plot
#' @method plot crm_fit
#' @export
plot.crm_fit <- function(x, pars = 'prob_tox', ...) {
  if(!is.null(x$fit))
    rstan::plot(x$fit, pars = pars, ...)
}

#' Obtain summary of an crm_fit
#'
#' @param object \code{\link{crm_fit}} object to summarise.
#' @param ... Extra parameters, passed onwards.
#'
#' @return A summary object.
#' @method summary crm_fit
#' @export
#' @seealso
#' \code{\link{stan_crm}}
summary.crm_fit <- function(object, ...) {
  if(!is.null(object$fit))
    rstan::summary(object$fit, ...)
}

#' @title Convert \code{\link{crm_fit}} to instance of
#' \code{\link[coda]{mcmc.list}}
#'
#' @description This function allows trialr to use tidybayes functions.
#'
#' @param crm_fit Object of class \code{\link{crm_fit}}
#' @param ... Extra variables that are passed onwards.
#'
#' @return Object of class \code{\link[coda]{mcmc.list}}
#' @method as.mcmc.list crm_fit
#'
#' @importFrom coda as.mcmc.list
#' @importFrom coda mcmc
#' @importFrom rstan As.mcmc.list
#' @export
as.mcmc.list.crm_fit <- function(crm_fit, ...) {
  if(is.null(crm_fit$fit)) {
    as.mcmc.list(mcmc(crm_fit$samples))
  } else {
    As.mcmc.list(crm_fit$fit, ...)
  }
}
