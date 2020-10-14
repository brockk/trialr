#' @title Get data to run the PePS2 trial example
#'
#' @description Get data to run the BEBOP model in the PePS2 trial. The trial
#' investigates pembrolizumab in non-small-cell lung cancer. Patients may be
#' previously treated (PT) or treatment naive (TN). Pembro response rates in
#' lung cancer have been shown to increase with PD-L1 tumour proportion score.
#' PD-L1 score is measured at baseline. Each patient belongs to one of the Low,
#' Medium or High categories. These two baseline variables stratify the patient
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
#' @param num_patients Total number of patients to use, positive integer.
#' @param cohort_probs Probabilities that a patient belongs to each of the 6
#' cohorts, in the order given above; a vector of numbers between 0 and 1 that
#' add up to 1. \code{cohort_probs} or \code{cohort_rho} must be specified.
#' @param prob_eff Probabilities of efficacy in each of the 6 cohorts, in the
#' order given above; a vector of numbers between 0 and 1
#' @param prob_tox Probabilities of toxicity in each of the 6 cohorts, in the
#' order given above; a vector of numbers between 0 and 1
#' @param eff_tox_or Measure of strength of association between efficacy and
#' toxicity, in each of the 6 cohorts, in the order given above; a vector of
#' numbers. Use 1 for no association; numbers increasingly greater than 1 for
#' stronger positive associations, and numbers less than 1 for stronger negative
#' associations
#' @param cohort_rho Concentration parameters for cohort membership, in the
#' order given above, using a Dirichlet distribution. This leads to randomly-
#' sampled cohort sizes distributed Dir(cohort_rho). \code{cohort_probs} or
#' \code{cohort_rho} must be specified.
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
#'
#' @return a \code{list} of parameters
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' dat <- peps2_get_data(num_patients = 60,
#'                       prob_eff = c(0.167, 0.192, 0.5, 0.091, 0.156, 0.439),
#'                       prob_tox = rep(0.1, 6),
#'                       eff_tox_or = rep(1, 6))
#' fit <- stan_peps2(
#'   eff = dat$eff,
#'   tox = dat$tox,
#'   cohorts = dat$cohorts
#' )
#' }
peps2_get_data <- function(num_patients, cohort_probs = NULL,
                           prob_eff, prob_tox, eff_tox_or,
                           cohort_rho = c(15.7, 21.8, 12.4, 20.7, 18.0, 11.4),
                           alpha_mean = -2.2, alpha_sd = 2,
                           beta_mean = -0.5, beta_sd = 2,
                           gamma_mean = -0.5, gamma_sd = 2,
                           zeta_mean = -0.5, zeta_sd = 2,
                           lambda_mean = -2.2, lambda_sd = 2,
                           psi_mean = 0, psi_sd = 1) {

  if (is.null(cohort_probs))
    cohort_probs <- gtools::rdirichlet(1, cohort_rho)
  cohort_sizes <- c(stats::rmultinom(1, size = num_patients,
                                     prob = cohort_probs))
  cohorts <- rep(1:length(cohort_sizes), times = cohort_sizes)
  cohorts = factor(cohorts, levels = 1:length(cohort_rho))

  cohort.params = cbind(cohort_sizes, prob_eff, prob_tox, eff_tox_or)
  cohort.params = split(cohort.params, 1:length(cohort_sizes))
  outcomes <- lapply(cohort.params, function(x) ranBin2(x[1], x[2:3], psi=x[4]))
  outcomes <- do.call(rbind, outcomes)

  eff <- outcomes[, 1]
  tox <- outcomes[, 2]
  x1 <- as.integer(cohorts %in% 4:6)
  x2 <- as.integer(cohorts == 1 | cohorts == 4)
  x3 <- as.integer(cohorts == 2 | cohorts == 5)
  cohort_eff = unname(tapply(eff, cohorts, sum))
  cohort_eff[is.na(cohort_eff)] = 0
  cohort_tox = unname(tapply(tox, cohorts, sum))
  cohort_tox[is.na(cohort_tox)] = 0

  dat <- list(
    # Data
    J = num_patients,
    eff = eff,
    tox = tox,
    cohorts = cohorts,
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
  return(dat)
}
