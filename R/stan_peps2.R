#' Fit the P2TNE model developed for the PePS2 trial to some outcomes.
#'
#' @description The PePS2 trial
#' investigates pembrolizumab in non-small-cell lung cancer. Patients may be
#' previously treated (PT) or treatment naive (TN). Response rates in
#' lung cancer have been shown to increase with PD-L1 tumour proportion score.
#' PD-L1 score is measured at baseline. Each patient belongs to one of the
#' categories <1%, 1-49% or >=50%. PT vs TN status and PD-L1 category jointly
#' stratify the patient
#' population and are used as predictive variables to stratify the analysis.
#' The BEBOP model studies co-primary efficacy and toxicity outcomes in the
#' presence of predictive data. Thus, PePS2 studies efficacy and toxicity in
#' 6 distinct cohorts:
#' TN Low, TN Medium, TN High, PT Low, PT Medium, PT High.
#' The design admits all-comers and does not target specific sample sizes in the
#' individual cohorts.
#' Hyperprior parameters have defaults to match those used in PePS2, but all may
#' be overridden.
#' The returned object includes randomly-sampled outcomes, as well as parameters
#' to run the model. These are all combined in the same list object for passing
#' to RStan, as is the convention.
#' See the accompanying vignette for a full description.
#'
#' @param eff A vector of efficacy outcomes for the patients, where 1=efficacy
#' and 0=no efficacy.
#' @param tox A vector of toxicity outcomes for the patients, where 1=toxicity
#' and 0=no toxicity.
#' @param cohorts A vector of integers from 1 to 6, denoting the cohorts to
#' which the patients belong.
#' @param alpha_mean The prior mean of alpha. Alpha is the efficacy model
#' intercept.
#' @param alpha_sd The prior standard deviation of alpha. Alpha is the efficacy
#' model  intercept.
#' @param beta_mean The prior mean of beta. Beta is the efficacy model term
#' for being previously treated.
#' @param beta_sd The prior standard deviation of beta. Beta is the efficacy
#' model term for being previously treated.
#' @param gamma_mean The prior mean of gamma. Gamma is the efficacy model term
#' for being PD-L1 score = Low.
#' @param gamma_sd The prior standard deviation of gamma. Gamma is the efficacy
#' model term for being PD-L1 score = Low.
#' @param zeta_mean The prior mean of zeta. Zeta is the efficacy model term
#' for being PD-L1 score = Medium.
#' @param zeta_sd The prior standard deviation of zeta. Zeta is the efficacy
#' model term for being PD-L1 score = Medium.
#' @param lambda_mean The prior mean of lambda. Lambda is the toxicity model
#' intercept.
#' @param lambda_sd The prior standard deviation of lambda. Lambda is the
#' toxicity model intercept.
#' @param psi_mean The prior mean of psi. Psi is the joint model association
#' parameter.
#' @param psi_sd The prior standard deviation of psi. Psi is the joint model
#' association parameter.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores}, and
#' \code{control}.
#'
#' @return Object of class \code{\link[rstan:stanfit]{rstan::stanfit}} returned
#' by \code{\link[rstan:sampling]{rstan::sampling}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- stan_peps2(
#'   eff = c(0, 1, 0, 1, 0, 0),
#'   tox = c(0, 0, 1, 1, 0, 0),
#'   cohorts = c(3, 1, 1, 4, 5, 6)
#' )
#' }
stan_peps2 <- function(eff, tox, cohorts,
                       alpha_mean = -2.2, alpha_sd = 2,
                       beta_mean = -0.5, beta_sd = 2,
                       gamma_mean = -0.5, gamma_sd = 2,
                       zeta_mean = -0.5, zeta_sd = 2,
                       lambda_mean = -2.2, lambda_sd = 2,
                       psi_mean = 0, psi_sd = 1, ...) {

  if(length(eff) != length(tox))
    stop('eff, tox and cohort should be the same length')
  if(length(eff) != length(cohorts))
    stop('eff, tox and cohorts should be the same length')

  cohorts = factor(cohorts, levels = 1:6)

  x1 <- as.integer(cohorts %in% 4:6)
  x2 <- as.integer(cohorts == 1 | cohorts == 4)
  x3 <- as.integer(cohorts == 2 | cohorts == 5)

  cohort_sizes <- sapply(1:6, function(i) sum(i == cohorts))
  cohort_eff = unname(tapply(eff, cohorts, sum))
  cohort_eff[is.na(cohort_eff)] = 0
  cohort_tox = unname(tapply(tox, cohorts, sum))
  cohort_tox[is.na(cohort_tox)] = 0

  dat <- list(
    # Data
    J = length(eff),
    eff = eff,
    tox = tox,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    # Cohort counts
    cohort_n = cohort_sizes,
    cohort_eff = cohort_eff,
    cohort_tox = cohort_tox,
    # Hyperparameters
    alpha_mean = alpha_mean,
    alpha_sd = alpha_sd,
    beta_mean = beta_mean,
    beta_sd = beta_sd,
    gamma_mean = gamma_mean,
    gamma_sd = gamma_sd,
    zeta_mean = zeta_mean,
    zeta_sd = zeta_sd,
    lambda_mean = lambda_mean,
    lambda_sd = lambda_sd,
    psi_mean = psi_mean,
    psi_sd = psi_sd
  )

  rstan::sampling(stanmodels$BebopInPeps2, data = dat, ...)
}
