
# Note: nbg_fit does not have its own constructor, it just piggy-backs crm_fit
# because nbg_fit requires no further logic.

# Generics ----
#' Print nbg_fit object.
#'
#' @param x \code{nbg_fit} object to print.
#' @param ... Extra parameters, passed onwards.
#' @method print nbg_fit
#' @export
print.nbg_fit <- function(x, ...) {
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
