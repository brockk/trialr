


# EffTox
efftox_solve_p <- function(eff0, tox1, eff_star, tox_star) {
  # Calculate p for the efficacy/toxicity contours that will intersect points
  # (eff0, 0), (eff.star, tox.star), and (1, tox1)

  .objective = function(p, eff0, tox1, eff_star, tox_star) {
    a <- ((1 - eff_star) / (1 - eff0))
    b <- tox_star / tox1
    return(a^p + b^p - 1)
  }

  rt <- uniroot(.objective, interval = c(0, 100),
                eff0 = eff0, tox1 = tox1, eff_star = eff_star,
                tox_star = tox_star)
  return(rt$root)
}

efftox_utility <- function(p, eff0, tox1, prob_eff, prob_tox) {
  a <- ((1 - prob_eff) / (1 - eff0))
  b <- prob_tox / tox1
  r = (a^p + b^p) ^ (1/p)
  return(1 - r)
}

efftox_process <- function(dat, fit, p_e, p_t) {
  # Posterior mean estimates
  prob_eff <- colMeans(rstan::extract(fit, 'prob_eff')[[1]])
  prob_acc_eff <- colMeans(rstan::extract(fit, 'prob_acc_eff')[[1]])
  prob_tox <- colMeans(rstan::extract(fit, 'prob_tox')[[1]])
  prob_acc_tox <- colMeans(rstan::extract(fit, 'prob_acc_tox')[[1]])
  post_utility <- colMeans(rstan::extract(fit, 'utility')[[1]])
  # Derived quantities
  utility = efftox_utility(dat$p, dat$eff0, dat$tox1, prob_eff, prob_tox)
  # Dose admissibility
  dose_indices <- 1:dat$num_doses
  lowest <- min(dat$doses)
  highest <- max(dat$doses)
  in_range <- sapply(dose_indices,
                    function(x) (x >= lowest - 1) & (x <= highest + 1))
  acceptable <- (prob_acc_eff > p_e) & (prob_acc_tox > p_t) & in_range
  if(sum(acceptable) > 0) {
    recommended_dose <- which.max(ifelse(acceptable, utility, NA))  # 2
    # rec <- dose_indices == selected_dose
  } else {
    recommended_dose <- NA
    # rec <- rep(FALSE, length(dose_indices))
  }

  l <- list(dose_indices = dose_indices, recommended_dose = recommended_dose,
            prob_eff = prob_eff, prob_tox = prob_tox,
            prob_acc_eff = prob_acc_eff, prob_acc_tox = prob_acc_tox,
            utility = utility, post_utility = post_utility,
            acceptable = acceptable)
  class(l) <- "efftox_analysis"
  return(l)
}

efftox_analysis_to_df <- function(x) {
  df <- data.frame(DoseLevel = factor(x$dose_indices),
                   ProbEff = x$prob_eff, ProbTox = x$prob_tox,
                   ProbAccEff = x$prob_acc_eff, ProbAccTox = x$prob_acc_tox,
                   Utility = x$utility, Acceptable = x$acceptable)
  return(df)
}

efftox_simulate <- function(dat, num_sims, max_patients, first_dose,
                            eff0, tox1, p, p_e, p_t,
                            true_eff, true_tox,
                            cohort_size = NULL, cohort_sizes = NULL,
                            ...) {

  # Either cohort_size or cohort_sizes must be specified so we now
  # how often to re-evaluate dose.
  # cohort_sizes parameter allows varying cohort sizes
  if(is.null(cohort_size) & is.null(cohort_sizes))
    stop('specify either cohort_size (integer) or cohort_sizes (array of integers)')
  if(is.null(cohort_sizes))
    cohort_sizes = rep(cohort_size, max_patients / cohort_size)

  fit = NULL
  doses_picked = integer(length = num_sims)
  # What is sum(cohort_sizes) != max_patients?
  # efficacies = matrix(nrow = num_sims, ncol = max_patients)
  # toxicities = matrix(nrow = num_sims, ncol = max_patients)
  # doses_given = matrix(nrow = num_sims, ncol = max_patients)
  efficacies = list()
  toxicities = list()
  doses_given = list()

  for(i in 1:num_sims) {
    print(paste('Starting', i))
    this_dat = dat
    dose = first_dose
    cohort_sizes = rep(3, max_patients / 3)
    for(cohort_size in cohort_sizes) {
      if(this_dat$num_patients >= max_patients) break()
      prob_eff = true_eff[dose]
      prob_tox = true_tox[dose]
      # Simulate new efficacy events
      new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)
      new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)
      # And append to trial data
      this_dat$eff = c(this_dat$eff, new_eff)
      this_dat$tox = c(this_dat$tox, new_tox)
      # Also reflect doses delivered
      this_dat$doses <- c(this_dat$doses, rep(dose, cohort_size))
      this_dat$num_patients = this_dat$num_patients + cohort_size
      fit <- efftox_fit(this_dat, fit = fit, ...)
      l <- efftox_process(this_dat, fit, p_e = p_e, p_t = p_e)
      # Select a dose?
      if(sum(l$admissible) > 0) {
        # Select dose
        dose = which.max(ifelse(l$admissible, l$utility, NA))
      } else {
        dose <- NA
        break()
      }
    }
    doses_picked[i] = dose
    efficacies[[i]] = this_dat$eff
    toxicities[[i]] = this_dat$tox
    doses_given[[i]] = this_dat$doses
  }

  return(list(doses_picked = doses_picked,
              efficacies = efficacies,
              toxicities = toxicities,
              doses_given = doses_given))
}

