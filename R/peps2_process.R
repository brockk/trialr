
#' @title Process RStan samples from a BEBOP model fit to PePS2 data
#'
#' @description Process RStan samples from a BEBOP model fit to PePS2 data. This
#' step lets us make inferences about whether the modelled efficacy and toxicity
#' probabilities suggest the treatment is acceptable in each of the cohorts
#' under study.
#' The parameters have default values to match those used in the PePS2 trial.
#' See the accompanying vignette for a full description.
#'
#' @param fit An instance of \code{rstan::stanmodel}, derived by fitting data to
#' the BEBOP in PePS2 model.
#' Use \code{stan_peps2}.
#' @param min_eff The lower efficacy probability threshold; a number between 0
#' and 1.
#' @param max_tox The upper toxicity probability threshold; a number between 0
#' and 1.
#' @param eff_cert Certainty required to infer the treatment is acceptable with
#' regards to being probably efficacious; a number between 0 and 1.
#' @param tox_cert Certainty required to infer the treatment is acceptable with
#' regards to being probably tolerable; a number between 0 and 1.
#'
#' @return a list with the following items:
#' \itemize{
#' \item \code{ProbEff}, the posterior mean probability of efficacy in the 6
#' cohorts.
#' \item \code{ProbAccEff}, the posterior mean probability that the probability
#' of efficacy exceeds \code{min_eff}, in the 6 cohorts.
#' \item \code{ProbTox}, the posterior mean probability of toxicity in the 6
#' cohorts.
#' \item \code{ProbAccTox}, the posterior mean probability that the probability
#' of toxicity is less than \code{max_tox}, in the 6 cohorts.
#' \item \code{Accept}, a vector of logical values to show whether treatment
#' should be accepted in the 6 cohorts. Treatment is acceptable when it is
#' probably efficacious and probably not toxic, with respect to the described
#' rules.
#' \item \code{alpha}, the posterior mean estimate of alpha.
#' \item \code{beta}, the posterior mean estimate of beta.
#' \item \code{gamma}, the posterior mean estimate of gamma.
#' \item \code{zeta}, the posterior mean estimate of zeta.
#' \item \code{lambda}, the posterior mean estimate of lambda.
#' \item \code{psi}, the posterior mean estimate of psi.
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' fit <- stan_peps2(
#'   eff = c(0, 1, 0, 1, 0, 0),
#'   tox = c(0, 0, 1, 1, 0, 0),
#'   cohorts = c(3, 1, 1, 4, 5, 6)
#' )
#' decision <- peps2_process(fit)
#' decision$Accept
#' decision$ProbEff
#' decision$ProbAccEff
#'
#' @seealso
#' \code{\link{peps2_get_data}}
peps2_process <- function(fit, min_eff = 0.1, max_tox = 0.3,
                          eff_cert = 0.7, tox_cert = 0.9) {
  #fit = samp
  acc_eff <- rstan::extract(fit, par = 'prob_eff')[[1]] > min_eff
  acc_tox <- rstan::extract(fit, par = 'prob_tox')[[1]] < max_tox
  accept <- (apply(acc_eff, 2, mean) > eff_cert) &
    (apply(acc_tox, 2, mean) > tox_cert)
  prob_eff <- colMeans(rstan::extract(fit, 'prob_eff')[[1]])
  prob_tox <- colMeans(rstan::extract(fit, 'prob_tox')[[1]])
  l <- list(ProbEff = prob_eff, ProbAccEff = apply(acc_eff, 2, mean),
            ProbTox = prob_tox, ProbAccTox = apply(acc_tox, 2, mean),
            Accept = accept)
  # Append posterior parameter means
  l <- append(l, lapply(rstan::extract(fit, pars=c('alpha', 'beta', 'gamma',
                                                   'zeta', 'lambda')), mean))
  # Add psi posterior mean
  if('psi' %in% names(fit))
    l <- append(l, lapply(rstan::extract(fit, pars=c('psi')), mean))
  return(l)
}
