
# source('LearningRStan/package/SimFuncs.R')
# source('LearningRStan/package/phase2/bebop/peps2/categorical-pdl1/CorrelatedBinaryOutcomes.R')


peps2_get_data <- function(num_patients, cohort_probs = NULL,
                           prob_eff, prob_tox, eff_tox_or,
                           cohort_rho = c(15.7, 21.8, 12.4, 20.7, 18.0, 11.4),
                           alpha_mean = -2.2, alpha_sd = 2,
                           beta_mean = -0.5, beta_sd = 2,
                           gamma_mean = -0.5, gamma_sd = 2,
                           zeta_mean = -0.5, zeta_sd = 2,
                           lambda_mean = -2.2, lambda_sd = 2,
                           psi_mean = 0, psi_sd = 1) {

  # Get data to pass to RStan for a simulated PePS2 trial iteration.
  # Cohort memberships, and efficacy & toxicity outcomes are randomly sampled.
  # Hyperparameters are specifiable but have defaults matching publication.

  # cohort_rho or cohort_probs must be specified

  if (is.null(cohort_probs))
    cohort_probs <- gtools::rdirichlet(1, cohort_rho)
  cohort_sizes <- c(rmultinom(1, size = num_patients, prob = cohort_probs))
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


peps2_trial_decision <- function(dat, fit, min_eff = 0.1, max_tox = 0.3,
                                 eff_cert = 0.7, tox_cert = 0.9) {
  # Provide posterior mean probability of efficacy and toxicity
  # and acceptance decision for each cohort
  acc_eff <- extract(fit, par = 'prob_eff')[[1]] > min_eff
  acc_tox <- extract(fit, par = 'prob_tox')[[1]] < max_tox
  accept <- (apply(acc_eff, 2, mean) > eff_cert) & mean(acc_tox) > tox_cert
  prob_eff <- unname(summary(fit, pars='prob_eff')$summary[, 'mean'])
  prob_tox <- unname(summary(fit, pars='prob_tox')$summary[, 'mean'])
  l <- list(ProbEff=prob_eff, ProbTox=prob_tox, Decision=accept)
  # Append posterior parameter means
  l <- append(l, lapply(extract(fit, pars=c('alpha', 'beta', 'gamma', 'zeta', 'lambda')), mean))
  # Add psi posterior mean
  if('psi' %in% names(fit))
    l <- append(l, lapply(extract(fit, pars=c('psi')), mean))
  return(l)
}


peps2_run_sims <- function(num_sims, sample_data_func, summarise_func, ...) {
  dat <- sample_data_func()
  model <- stanmodels$BebopInPeps2
  return(.run.sims(model, num_sims, sample_data_func, summarise_func, ...))
}
