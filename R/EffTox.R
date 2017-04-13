


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

efftox_parameters_demo <- function() {
  # Demonstration from 'Effective sample size for computing prior
  # hyperparameters in Bayesian phase I-II dose-finding', Thall et al., 2014
  eff0 = 0.5
  tox1 = 0.65
  eff_star = 0.7
  tox_star = 0.25
  p = efftox_solve_p(eff0, tox1, eff_star, tox_star)
  x <- list(
    num_doses = 5,
    real_doses = c(1, 2, 4, 6.6, 10),
    efficacy_hurdle = 0.5,
    toxicity_hurdle = 0.3,
    p = p,
    eff0 = eff0,
    tox1 = tox1,
    eff_star = 0.7,
    tox_star = 0.25,

    alpha_mean = -7.9593, alpha_sd = 3.5487,
    beta_mean = 1.5482, beta_sd = 3.5018,
    gamma_mean = 0.7367, gamma_sd = 2.5423,
    zeta_mean = 3.4181, zeta_sd = 2.4406,
    eta_mean = 0, eta_sd = 0.2,
    psi_mean = 0, psi_sd = 1,

    eff = c(),
    tox = c(),
    doses = c(),
    num_patients = 0

  )
  return(x)
}

efftox_utility <- function(p, eff0, tox1, prob_eff, prob_tox) {
  a <- ((1 - prob_eff) / (1 - eff0))
  b <- prob_tox / tox1
  r = (a^p + b^p) ^ (1/p)
  return(1 - r)
}

efftox_process <- function(dat, fit, p_e, p_t) {
  # fit, an object of class stanfit, yielded by a call to rstan::sampling(stanmodels$EffTox, data = dat)
  #       where dat is a list of your data. See TODO[that explains dat]

  # Posterior mean estimates
  prob_eff <- colMeans(rstan::extract(fit, 'prob_eff')[[1]])
  prob_acc_eff <- colMeans(rstan::extract(fit, 'prob_acc_eff')[[1]])
  prob_tox <- colMeans(rstan::extract(fit, 'prob_tox')[[1]])
  prob_acc_tox <- colMeans(rstan::extract(fit, 'prob_acc_tox')[[1]])
  post_utility <- colMeans(rstan::extract(fit, 'utility')[[1]])
  # Derived quantities
  utility = efftox_utility(dat$p, dat$eff0, dat$tox1,
                           prob_eff, prob_tox)
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

efftox_simulate <- function(dat, num_sims, first_dose, p_e, p_t,
                            true_eff, true_tox, cohort_sizes, ...) {

  recommended_dose <- integer(length = num_sims)
  efficacies <- list()
  toxicities <- list()
  doses_given <- list()

  for(i in 1:num_sims) {
    print(paste('Starting iteration', i))
    this_dat <- dat
    dose <- first_dose
    for(cohort_size in cohort_sizes) {
      prob_eff <- true_eff[dose]
      prob_tox <- true_tox[dose]
      # Simulate new efficacy events
      new_eff <- rbinom(n = cohort_size, size = 1, prob = prob_eff)
      new_tox <- rbinom(n = cohort_size, size = 1, prob = prob_tox)
      # And append to trial data
      this_dat$eff <- c(this_dat$eff, new_eff)
      this_dat$tox <- c(this_dat$tox, new_tox)
      # Also reflect doses delivered
      this_dat$doses <- c(this_dat$doses, rep(dose, cohort_size))
      this_dat$num_patients <- this_dat$num_patients + cohort_size
      samp <- rstan::sampling(stanmodels$EffTox, data = this_dat, ...)
      l <- efftox_process(this_dat, samp, p_e = p_e, p_t = p_e)
      # Select a dose?
      if(sum(l$acceptable) > 0) {
        # Select dose
        dose = which.max(ifelse(l$acceptable, l$utility, NA))
      } else {
        dose <- NA
        break()
      }
    }
    recommended_dose[i] = dose
    efficacies[[i]] = this_dat$eff
    toxicities[[i]] = this_dat$tox
    doses_given[[i]] = this_dat$doses
  }

  return(list(recommended_dose = recommended_dose,
              efficacies = efficacies,
              toxicities = toxicities,
              doses_given = doses_given))
}

efftox_get_tox <- function(eff, util, p, eff0, tox1) {

  a = ((1 - eff) / (1 - eff0))
  return(tox1 * ((1 - util)^p - a^p)^(1 / p))
}

efftox_contour_plot <- function(dat,
                                use_ggplot = FALSE,
                                prob_eff = NULL, prob_tox = NULL,
                                n = 1000, util_lower = -3,
                                util_upper = 3, util_delta = 0.2) {
  eff_vals = seq(0, 1, length.out = n)
  util_vals = seq(util_lower, util_upper, by = util_delta)


  if(use_ggplot) {
    stop('Not implemented')
  } else {
    plot(NULL, ylim = c(0, 1), xlim = c(0, 1), ylab = 'Prob(Toxicity)',
         xlab = 'Prob(Efficacy)')

    for(u in util_vals) {
      tox_vals = efftox_get_tox(eff_vals, u, dat$p, dat$eff0, dat$tox1)
      points(eff_vals, tox_vals, type = 'l', col = 'grey', lwd = 0.2)
    }

    # # Add neutral utility contour
    tox_vals = efftox_get_tox(eff_vals, 0, dat$p, dat$eff0, dat$tox1)
    points(eff_vals, tox_vals, type = 'l', col = 'black', lwd = 2)

    # Add hinge points
    points(dat$eff0, 0, col = 'blue', pch = 2)
    points(1, dat$tox1, col = 'blue', pch = 2)
    points(dat$eff_star, dat$tox_star, col = 'blue', pch = 2)

    # # Add provided eff & tox points
    if(!is.null(prob_eff) & !is.null(prob_tox)) {
      points(prob_eff, prob_tox, col = 'red',
             pch = as.character(1:length(prob_eff)))
    }
  }
}

