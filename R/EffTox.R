
#' @title Calculate the p-index for EffTox utility contours
#'
#' @description
#' Calculate the p-index for EffTox utility contours so that the neutral utility
#' contour intersects the following points in the
#' Prob(Efficacy) - Prob(Toxicity) plane:
#' (\code{eff0}, 0), (1, \code{tox1}) and (\code{eff_star}, \code{tox_star})
#'
#' @param eff0 Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1
#' @param tox1 Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1
#' @param eff_star Efficacy probability of an equi-utility third point
#' @param tox_star Toxicity probability of an equi-utility third point
#'
#' @return The p-index
#'
#' @export
#'
#' @examples
#' efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' @references Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample
#' size for computing prior hyperparameters in Bayesian phase I-II dose-finding
efftox_solve_p <- function(eff0, tox1, eff_star, tox_star) {
  # Calculate p for the efficacy/toxicity contours that will intersect points
  # (eff0, 0), (eff.star, tox.star), and (1, tox1)

  if(any(eff0 <= 0, tox1 <= 0, eff_star <= 0, tox_star <= 0,
         eff0 >= 1, tox1 >= 1, eff_star >= 1, tox_star >= 1)) {
    stop('eff0, tox1, eff_star and tox_star must all be between 0 and 1.')
  }

  .objective = function(p, eff0, tox1, eff_star, tox_star) {
    a <- ((1 - eff_star) / (1 - eff0))
    b <- tox_star / tox1
    return(a^p + b^p - 1)
  }

  rt <- stats::uniroot(.objective, interval = c(0, 100),
                       eff0 = eff0, tox1 = tox1, eff_star = eff_star,
                       tox_star = tox_star)
  return(rt$root)
}

#' Container class for parameters to fit the EffTox model in trialr.
#'
#' @name efftox_params-class
#' @aliases efftox_params
#' @docType class
#'
#' @param real_doses a vector of numbers.The doses under investigation.
#' They should be ordered from lowest to highest and be in consistent units.
#' E.g., to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
#' c(10, 20, 50).
#' @param efficacy_hurdle Minimum acceptable efficacy probability.
#' A number between 0 and 1.
#' @param toxicity_hurdle Maximum acceptable toxicity probability.
#' A number between 0 and 1.
#' @param p_e Certainty required to infer a dose is acceptable with regards to
#' being probably efficacious; a number between 0 and 1.
#' @param p_t Certainty required to infer a dose is acceptable with regards to
#' being probably tolerable; a number between 0 and 1.
#' @param eff0 Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1 (see Details).
#' @param tox1 Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1 (see Details).
#' @param eff_star Efficacy probability of an equi-utility third point (see
#' Details).
#' @param tox_star Toxicity probability of an equi-utility third point (see
#' Details).
#' @param alpha_mean The prior normal mean of the intercept term in the toxicity
#' logit model. A number.
#' @param alpha_sd The prior normal standard deviation of the intercept term in
#' the toxicity logit model. A number.
#' @param beta_mean The prior normal mean of the slope term in the toxicity
#' logit model. A number.
#' @param beta_sd The prior normal standard deviation of the slope term in the
#' toxicity logit model. A number.
#' @param gamma_mean The prior normal mean of the intercept term in the efficacy
#' logit model. A number.
#' @param gamma_sd The prior normal standard deviation of the intercept term in
#' the efficacy logit model. A number.
#' @param zeta_mean The prior normal mean of the slope term in the efficacy
#' logit model. A number.
#' @param zeta_sd The prior normal standard deviation of the slope term in the
#' efficacy logit model. A number.
#' @param eta_mean The prior normal mean of the squared term coefficient in the
#' efficacy logit model. A number.
#' @param eta_sd The prior normal standard deviation of the squared term
#' coefficient in the efficacy logit model. A number.
#' @param psi_mean The prior normal mean of the association term in the combined
#' efficacy-toxicity model. A number.
#' @param psi_sd The prior normal standard deviation of the association term in
#' the combined efficacy-toxicity model. A number.
#'
#' @export
#'
#' @seealso
#' \code{\link{stan_efftox}}
#' \code{\link{stan_efftox_demo}}
efftox_params <- function(real_doses, efficacy_hurdle, toxicity_hurdle,
                          p_e, p_t, eff0, tox1, eff_star, tox_star,
                          alpha_mean, alpha_sd, beta_mean, beta_sd,
                          gamma_mean, gamma_sd, zeta_mean, zeta_sd,
                          eta_mean, eta_sd, psi_mean, psi_sd) {
  # efftox_params class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )

  p <- efftox_solve_p(eff0, tox1, eff_star, tox_star)
  x <- loo::nlist(real_doses, num_doses = length(real_doses),
                  efficacy_hurdle, toxicity_hurdle,
                  p_e, p_t, p, eff0, tox1, eff_star, tox_star,
                  alpha_mean, alpha_sd, beta_mean, beta_sd, gamma_mean, gamma_sd,
                  zeta_mean, zeta_sd, eta_mean, eta_sd, psi_mean, psi_sd, version
  )

  # Initialise with no patients observed
  x$doses = c()
  x$eff = c()
  x$tox = c()
  x$num_patients = 0

  # Set type. This is, at heart, just a list.
  class(x) <- c("efftox_params", "list")
  return(x)
}

#' @title Get parameters to run the EffTox demo
#'
#' @description Get parameters to run the EffTox demo. These match those used
#' to demonstrate EffTox in Thall et al. 2014.
#'
#' @return a \code{list} of parameters, described in \code{efftox_params}
#'
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' names(dat)
#' dat$real_doses == c(1, 2, 4, 6.6, 10)
#'
#' @seealso
#' \code{\link{efftox_params}}
#'
#' @references Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample
#' size for computing prior hyperparameters in Bayesian phase I-II dose-finding
efftox_parameters_demo <- function() {
  # Demonstration from 'Effective sample size for computing prior
  # hyperparameters in Bayesian phase I-II dose-finding', Thall et al., 2014
  x <- efftox_params(real_doses = c(1, 2, 4, 6.6, 10),
                     efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                     p_e = 0.1, p_t = 0.1, eff0 = 0.5, tox1 = 0.65,
                     eff_star = 0.7, tox_star = 0.25,
                     alpha_mean = -7.9593, alpha_sd = 3.5487,
                     beta_mean = 1.5482, beta_sd = 3.5018,
                     gamma_mean = 0.7367, gamma_sd = 2.5423,
                     zeta_mean = 3.4181, zeta_sd = 2.4406,
                     eta_mean = 0, eta_sd = 0.2,
                     psi_mean = 0, psi_sd = 1)

  return(x)
}


#' @title Get the utility of efficacy & toxicity probability pairs
#'
#' @description Get the utility of efficacy & toxicity probability pairs
#'
#' @param p p-index of EffTox utility contours. Use \code{efftox_solve_p}
#' @param eff0 Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1
#' @param tox1 Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1
#' @param prob_eff Probability of efficacy; number between 0 and 1
#' @param prob_tox Probability of toxicity; number between 0 and 1
#'
#' @return Utility value(s)
#'
#' @export
#'
#' @examples
#' p <- efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' u <- efftox_utility(p, 0.5, 0.65, prob_eff = 0.7, prob_tox = 0.25)
#' round(u, 4) == 0
#'
#' u <- efftox_utility(p, 0.5, 0.65, prob_eff = c(0.6, 0.7, 0.8),
#'                     prob_tox = c(0.1, 0.2, 0.3))
#' round(u, 2) == c(0.04, 0.08, 0.12)
#'
#' @seealso \code{\link{efftox_solve_p}}
efftox_utility <- function(p, eff0, tox1, prob_eff, prob_tox) {
  a <- ((1 - prob_eff) / (1 - eff0))
  b <- prob_tox / tox1
  r = (a^p + b^p) ^ (1/p)
  return(1 - r)
}


#' @title Process RStan samples from an EffTox model
#'
#' @description Internal function to process rstan samples from an EffTox model
#' to make inferences about dose-acceptability, dose-utility and which dose
#' should be recommended next.
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param fit An instance of \code{rstan::stanmodel}, derived by fitting the
#' trialr EffTox model.
#' @return An instance of \code{\link{efftox_fit}}.
#'
efftox_process <- function(dat, fit) {
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
  acceptable <- (prob_acc_eff > dat$p_e) & (prob_acc_tox > dat$p_t) & in_range
  if(sum(acceptable) > 0) {
    recommended_dose <- which.max(ifelse(acceptable, utility, NA))  # 2
  } else {
    recommended_dose <- NA
  }

  x <- efftox_fit(dose_indices, recommended_dose, prob_eff, prob_tox,
                  prob_acc_eff, prob_acc_tox, utility, post_utility,
                  acceptable, dat, fit)
  return(x)
}


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
  df <- data.frame(DoseLevel = factor(x$dose_indices),
                   ProbEff = x$prob_eff, ProbTox = x$prob_tox,
                   ProbAccEff = x$prob_acc_eff, ProbAccTox = x$prob_acc_tox,
                   Utility = x$utility, Acceptable = x$acceptable)
  return(df)
}

#' @title Run EffTox simulations
#'
#' @description Run EffTox simulations for assumed true efficacy and toxicity
#' curves.
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param num_sims integer, number of simulated iterations
#' @param first_dose integer, the dose-level to give to patient 1, e.g. 1 for
#' the lowest dose.
#' @param true_eff the true probabilities of efficacy at the doses under
#' investigation; a vector of numbers between 0 and 1.
#' @param true_tox the true probabilities of toxicity at the doses under
#' investigation; a vector of numbers between 0 and 1.
#' @param cohort_sizes a vector of integer cohort sizes. A dose decision is made
#' when each cohort is completed and the next cohort is treated at the
#' recommended dose. To conduct a trial using at most 20 patients, where dose is
#' re-evaluated after every second patient, use \code{rep(2, 10)}. To conduct a
#' trial of 8 patients where dose is re-evaluated after each single patient, use
#' \code{rep(1, 8)}. Cohort size need not be uniform. E.g.
#' \code{c(rep(1, 5), rep(3, 10))} represents a trial where the dose is
#' re-evaluated after each patient for the first 5 patients, and then after
#' every third patient for a further 30 patients.
#' @param ... Extra parameters provided via the ellipsis are passed to
#' \code{stan::sampling}
#'
#' @return A list with named elements \code{recommended_dose},
#' \code{efficacies}, \code{toxicities}, and \code{doses_given}.
#'
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' set.seed(123)
#' # Let's say we want to use only 2 chains. Extra args are passed to stan
#' \dontrun{
#' sims <- efftox_simulate(dat, num_sims = 10, first_dose = 1,
#'                         true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
#'                         true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
#'                         cohort_sizes = rep(3, 13),
#'                         chains = 2)
#' table(sims$recommended_dose) / length(sims$recommended_dose)
#' table(unlist(sims$doses_given)) / length(unlist(sims$doses_given))
#' table(unlist(sims$doses_given)) / length(sims$recommended_dose)
#' }
#' # In real life, we would run thousands of iterations, not 10.
#' # This is an example.
efftox_simulate <- function(dat, num_sims, first_dose, true_eff, true_tox,
                            cohort_sizes, ...) {

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
      new_eff <- stats::rbinom(n = cohort_size, size = 1, prob = prob_eff)
      new_tox <- stats::rbinom(n = cohort_size, size = 1, prob = prob_tox)
      # And append to trial data
      this_dat$eff <- c(this_dat$eff, new_eff)
      this_dat$tox <- c(this_dat$tox, new_tox)
      # Also reflect doses delivered
      this_dat$doses <- c(this_dat$doses, rep(dose, cohort_size))
      this_dat$num_patients <- this_dat$num_patients + cohort_size
      samp <- rstan::sampling(stanmodels$EffTox, data = this_dat, ...)
      l <- efftox_process(this_dat, samp)
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

#' @title Plot EffTox utility contours
#'
#' @description Plot EffTox utility contours. The probability of efficacy is
#' on the x-axis and toxicity on the y-axis. The zero-utility curve is plotted
#' bolder. The three "hinge points" are plotted as blue triangles. Optional
#' Prob(Efficacy) vs Prob(Toxicity) points can be added; these are shown as
#' red numerals, enumerated in the order provided.
#'
#' @param fit An instance of \code{\link{efftox_fit}}.
#' @param use_ggplot logical, TRUE to use ggplot2. Defaults to FALSE to use
#'  standard R graphics.
#' @param prob_eff vector of numbers between 0 and 1, containing the
#' efficacy probabilities of extra points to add to the plot as points,
#' e.g. the posterior mean efficacy probabilities of the doses under
#' investigation. Paired with prob_tox, thus they should be the same length.
#' Defaults to the values fitted by the model. Use NULL to supress.
#' @param prob_tox vector of numbers between 0 and 1, containing the
#' toxicity probabilities of extra points to add to the plot as points,
#' e.g. the posterior mean toxicity probabilities of the doses under
#' investigation. Paired with prob_eff, thus they should be the same length.
#' Defaults to the values fitted by the model. Use NULL to supress.
#' @param num_points integer for number of points to calculate on each curve.
#' The default is 1000 and this should be plenty.
#' @param util_vals A contour is plotted for each of these utility values.
#' The default is contours spaced by 0.2 between from -3 and 3,
#' i.e. \code{seq(-3, 3, by = 0.2)}.
#'
#' @return if \code{use_ggplot = TRUE}, an instance of \code{ggplot}; else no
#' object is returned. Omit assignment in either case to just view the plot.
#'
#' @export
#'
#' @examples
#' fit <- stan_efftox_demo(outcome_str = '1N 2E 3B')
#' efftox_contour_plot(fit)
#' title('EffTox utility contours')
#' # The same with ggplot2
#' efftox_contour_plot(fit, use_ggplot = TRUE) +
#'                     ggplot2::ggtitle('EffTox utility contours')
#'
#' @seealso
#' \code{\link{stan_efftox}}
efftox_contour_plot <- function(fit,
                                use_ggplot = FALSE,
                                prob_eff = fit$prob_eff,
                                prob_tox = fit$prob_tox,
                                num_points = 1000,
                                util_vals = seq(-3, 3, by = 0.2)) {
  eff_vals = seq(0, 1, length.out = num_points)
  dat <- fit$dat

  if(!is.null(prob_eff) & !is.null(prob_tox))
    if(length(prob_eff) != length(prob_tox))
      stop('prob_eff and prob_tox should be the same length')

  if(use_ggplot) {
    tox_vals = sapply(util_vals, function(u) efftox_get_tox(eff_vals, u, dat$p,
                                                            dat$eff0, dat$tox1))
    df = data.frame(eff_vals = rep(eff_vals, times = length(util_vals)),
                    tox_vals = as.numeric(tox_vals),
                    util_vals = rep(util_vals, each = length(eff_vals)))

    plt <- ggplot2::ggplot(df, ggplot2::aes(x = eff_vals, y = tox_vals,
                                            group = as.factor(util_vals))) +
      ggplot2::geom_line(size = 0.5, alpha = 0.25) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::xlab('Prob(Efficacy)') +
      ggplot2::ylab('Prob(Toxicity)')

    # Add neutral utility contour
    tox_vals = efftox_get_tox(eff_vals, 0, dat$p, dat$eff0, dat$tox1)
    df2 = data.frame(eff_vals, tox_vals, util_vals = 0)
    plt <- plt + ggplot2::geom_line(data = df2, size = 1)

    # Add hinge points
    df3 <- data.frame(prob_eff = c(dat$eff0, 1, dat$eff_star),
                      prob_tox = c(0, dat$tox1, dat$tox_star))
    plt <- plt + ggplot2::geom_point(data = df3, ggplot2::aes(x = prob_eff,
                                                              y = prob_tox,
                                                              group = 1),
                                     col = 'blue', shape = 24, size = 3)

    # Add provided eff & tox points
    if(!is.null(prob_eff) & !is.null(prob_tox)) {
      df4 <- data.frame(prob_eff, prob_tox, dl = 1:length(prob_eff))
      dl <- NULL
      plt <- plt + ggplot2::geom_text(data = df4, ggplot2::aes(x = prob_eff,
                                                               y = prob_tox,
                                                               group = 1,
                                                               label = dl),
                                      col = 'red', size = 4)
    }
    return(plt)
  } else {
    graphics::plot(NULL, ylim = c(0, 1), xlim = c(0, 1), ylab = 'Prob(Toxicity)',
                   xlab = 'Prob(Efficacy)')

    for(u in util_vals) {
      tox_vals = efftox_get_tox(eff_vals, u, dat$p, dat$eff0, dat$tox1)
      graphics::points(eff_vals, tox_vals, type = 'l', col = 'grey', lwd = 0.2)
    }

    # Add neutral utility contour
    tox_vals = efftox_get_tox(eff_vals, 0, dat$p, dat$eff0, dat$tox1)
    graphics::points(eff_vals, tox_vals, type = 'l', col = 'black', lwd = 2)

    # Add hinge points
    graphics::points(dat$eff0, 0, col = 'blue', pch = 2)
    graphics::points(1, dat$tox1, col = 'blue', pch = 2)
    graphics::points(dat$eff_star, dat$tox_star, col = 'blue', pch = 2)

    # # Add provided eff & tox points
    if(!is.null(prob_eff) & !is.null(prob_tox)) {
      graphics::points(prob_eff, prob_tox, col = 'red',
                       pch = as.character(1:length(prob_eff)))
    }
  }
}


#' @title Plot densities of EffTox dose utilities
#'
#' @description Plot densities of EffTox dose utilities. Optionally plot only a
#' subset of the doses by specifying the \code{doses} parameter. This function
#' requires ggplot2 be installed.
#'
#' @param fit An instance of \code{efftox_fit}.
#' @param doses optional, vector of integer dose-levels to plot. E.g. to plot
#' only dose-levels 1, 2 & 3 (and suppress the plotting of any other doses), use
#' \code{doses = 1:3}
#'
#' @return an instance of \code{ggplot}. Omit assignment to just view the plot.
#'
#' @export
#'
#' @note This function requires that ggplot2 be installed.
#'
#' @examples
#' fit <- stan_efftox_demo('1N 2E 3B')
#' efftox_utility_density_plot(fit) + ggplot2::ggtitle('My doses')  # Too busy?
#' # Specify subset of doses to make plot less cluttered
#' efftox_utility_density_plot(fit, doses = 1:3) + ggplot2::ggtitle('My doses')
#'
efftox_utility_density_plot <- function(fit, doses = NULL) {
  if(!('ggplot2' %in% utils::installed.packages()))
    stop('This function requires ggplot2 be installed.')

  u <- rstan::extract(fit$fit, par = 'utility')[[1]]
  df <- data.frame(Utility = as.numeric(u),
                   D = rep(1:5, each = nrow(u))
  )
  df$Dose = factor(df$D)
  if(!is.null(doses))
    df = df[df$D %in% doses, ]
  Dose <- Utility <- NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Utility, group = Dose,
                                        colour = Dose)) +
    ggplot2::geom_density()
  return(p)
}


#' @title Get dose-superiority matrix in EffTox
#'
#' @description Get a dose-superiority matrix from an EffTox dose analysis.
#' EffTox seeks to choose the dose with the highest utility, thus superiority
#' is inferred by posterior utility. The item in row i, col j is the posterior
#' probability that the utility of dose j exceeds that of dose i.
#'
#' @param fit An instance of \code{efftox_fit}.
#'
#' @return n by n matrix, where n is number of doses under investigation.
#' The item in row i, col j is the posterior probability that the utility of
#' dose j exceeds that of dose i.
#'
#' @export
#'
#' @examples
#' fit <- stan_efftox_demo('1N 2E 3B')
#' sup_mat <- efftox_superiority(fit)
efftox_superiority <- function(fit) {
  u <- rstan::extract(fit$fit, par = 'utility')[[1]]
  superiority_mat <- sapply(1:ncol(u), function(i) sapply(1:ncol(u), function(j)
    mean(u[ , i] > u[ , j])))
  diag(superiority_mat) <- NA
  dimnames(superiority_mat) = list(1:ncol(u), 1:ncol(u))
  return(superiority_mat)
}

#' @title Parse a string of EffTox outcomes to binary vector notation.
#'
#' @description Parse a string of EffTox outcomes to the binary vector notation
#' required by Stan for model invocation. The outcome string describes the doses
#' given and outcomes observed. The format of the string is described in Brock
#' et al. (2017). The letters E, T, N and B are used to represents patients that
#' experienced (E)fficacy only, (T)oxicity only, (B)oth efficacy and toxicity,
#' and (N)either. These letters are concatenated after numerical dose-levels to
#' convey the outcomes of cohorts of patients. For instance, \code{2ETB}
#' represents a cohort of three patients that were treated at dose-level 2, and
#' experienced efficacy, toxicity and both events, respectively. The results of
#' cohorts are separated by spaces. Thus, \code{2ETB 1NN} extends our previous
#' example, where the next cohort of two were treated at dose-level 1 and both
#' patients experienced neither efficacy nor toxicity. See examples.
#'
#' We present the notation in the EffTox setting but it is applicable in
#' general seamless phase I/II dose-finding scenarios.
#'
#' @param outcome_string character string, conveying doses given and outcomes
#' observed.
#' @param as.list TRUE (be default) to return a \code{list};
#' FALSE to return a \code{data.frame}
#'
#' @return If \code{as.list == TRUE}, a list with elements \code{eff}, \code{tox},
#' \code{doses} and \code{num_patients}. These elements are congruent with those
#' of the same name in \code{efftox_params}.
#' If \code{as.list == FALSE}, a data.frame with columns \code{eff}, \code{tox},
#' and \code{doses}.
#'
#' @export
#'
#' @examples
#' x = efftox_parse_outcomes('1NNE 2EEN 3TBB')
#' x$num_patients == 9
#' x$eff == c(0, 0, 1, 1, 1, 0, 0, 1, 1)
#' sum(x$tox) == 3
#'
#' @references
#' Brock, K., Billingham, L., Copland, M., Siddique, S., Sirovica, M., & Yap, C.
#' (2017). Implementing the EffTox dose-finding design in the Matchpoint trial.
#' BMC Medical Research Methodology, 17(1), 112.
#' https://doi.org/10.1186/s12874-017-0381-x
#'
efftox_parse_outcomes <- function(outcome_string, as.list = TRUE) {

  # Matching is done by regex.

  # This pattern ensures that outcome_string is valid. It is the gate-keeper.
  # It allows leading and trailing white space and demands >0 cohort strings.
  # e.g. "2ENT 3TT 2E "
  valid_str_match <- '^\\s*(\\d+[ETNB]+\\s*)+$'
  # This pattern identifies the individual cohort strings, e.g. 2ENT
  cohort_str_match <- '\\d+[ETNB]+'
  # This pattern extracts the dose-level from a cohort string, e.g. 2
  dl_str_match <- '\\d+'
  # And this pattern extracts the outcomes from a cohort string, e.g ENT
  outcomes_match_str <- '[ETNB]+'

  if(stringr::str_detect(outcome_string, valid_str_match)) {
    doses <- eff <- tox <- c()
    cohort_strs <- stringr::str_extract_all(
      outcome_string, cohort_str_match)[[1]]
    for(cohort_str in cohort_strs) {
      c_dl <- as.integer(stringr::str_extract(cohort_str, dl_str_match))
      if(c_dl <= 0) stop('Dose-levels must be strictly positive integers.')
      c_outcomes <- stringr::str_extract(cohort_str, outcomes_match_str)

      these_doses <- rep(c_dl, nchar(c_outcomes))
      doses <- c(doses, these_doses)

      these_outcomes = stringr::str_split(c_outcomes, '')[[1]]
      these_eff = as.integer((these_outcomes == 'E') | (these_outcomes == 'B'))
      eff = c(eff, these_eff)
      these_tox = as.integer((these_outcomes == 'T') | (these_outcomes == 'B'))
      tox = c(tox, these_tox)
    }
  } else {
    stop(paste0('"', outcome_string, '" is not a valid outcome string.
                A valid example is "1N 2EE 3TB 2BE"'))
  }

  if(as.list) {
    return(list(
      doses = doses, eff = eff, tox = tox, num_patients = length(doses)
    ))
  } else {
    return(data.frame(doses = doses, eff = eff, tox = tox))
  }
}


#' @title Calculate dose-transition pathways for an EffTox study
#'
#' @description Calculate dose-transition pathways for an EffTox study.
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param cohort_sizes vector of future cohort sizes, i.e. positive integers.
#' E.g. Tot calculate paths for the the next cohort of two followed by the next
#' cohort of three, use \code{c(2, 3)}.
#' @param next_dose the dose-level to be given to the immediately next cohort.
#' @param ... extra params passed to \code{rstan::sampling}.
#'
#' @return dose pathways in a \code{data.frame}.
#'
#' @export
#'
#' @examples
#' # Calculate the paths for the first cohort of 3 in Thall et al 2014 example
#' dat <- efftox_parameters_demo()
#' \dontrun{
#' dtps1 <- efftox_dtps(dat = dat, cohort_sizes = c(3), next_dose = 1)
#' }
#' # To calculate future paths in a partially-observed trial
#' dat <- efftox_parameters_demo()
#' dat$doses = array(c(1,1,1))
#' dat$eff = array(c(0,0,0))
#' dat$tox = array(c(1,1,1))
#' dat$num_patients = 3
#' \dontrun{
#' dtps2 <- efftox_dtps(dat = dat, cohort_sizes = c(3), next_dose = 1)
#' }
#'
#' @seealso
#' \code{\link{efftox_params}}
#'
#' \code{\link{efftox_parameters_demo}}
#'
#' @references Brock et al. (submitted 2017), Implementing the EffTox
#' Dose-Finding Design in the Matchpoint Trial.
efftox_dtps <- function(dat, cohort_sizes, next_dose, ...) {
  previous_doses = dat$doses
  previous_eff = dat$eff
  previous_tox = dat$tox
  previous_num_patients = dat$num_patients
  outcomes <- c('E', 'T', 'N', 'B')

  # Calculate feasible outcome combinations by cohort
  cohort_paths <- lapply(cohort_sizes,
                         function(x) gtools::combinations(n = 4, r = x,
                                                          v = outcomes,
                                                          repeats.allowed=TRUE))
  # Flatten cohort outcomes
  cohort_paths <- lapply(cohort_paths, function(x) apply(x, 1, paste0,
                                                         collapse = ''))

  # Calculate pathways
  cohort_paths <- expand.grid(cohort_paths, stringsAsFactors = FALSE)
  # Place to record dose recommendations
  dose_recs = matrix(nrow = nrow(cohort_paths), ncol = ncol(cohort_paths))
  # Cache DTP calculations to avoid needless repetition
  cache <- new.env()
  for(i in 1:nrow(cohort_paths)) {
    cohort_path <- cohort_paths[i,]
    cohort_dose <- next_dose
    dtp <- ""

    for(j in 1:length(cohort_path)) {
      dtp <- ifelse(nchar(dtp) > 0,
                    paste0(dtp, ' ', cohort_dose, cohort_path[j]),
                    paste0(cohort_dose, cohort_path[j])
      )
      if(dtp %in% names(cache)) {
        # Fetch from cache
        print(paste0('Fetching ', dtp, ' from cache'))
        cohort_dose <- cache[[dtp]]
      } else {
        # Calculate
        these_outcomes <- efftox_parse_outcomes(dtp)
        dat$doses <- array(c(previous_doses, these_outcomes$doses))
        dat$eff <- array(c(previous_eff, these_outcomes$eff))
        dat$tox <- array(c(previous_tox, these_outcomes$tox))
        dat$num_patients <- previous_num_patients +
          these_outcomes$num_patients
        print(paste0('Running ', dtp))
        fit <- rstan::sampling(stanmodels$EffTox, data = dat, ...)
        decision <- efftox_process(dat, fit)
        cohort_dose <- decision$recommended_dose
        # Cache
        cache[[dtp]] <- cohort_dose
      }
      dose_recs[i, j] <- cohort_dose
    }
  }

  df <- data.frame(D0 = rep(next_dose, nrow(cohort_paths)))
  for(k in 1:ncol(cohort_paths)) {
    df[, paste0('C', k - 1)] = cohort_paths[, k]
    df[, paste0('D', k)] = dose_recs[, k]
  }
  return(df)
}

#' Class of model fit by \pkg{trialr} using the EffTox dose-finding design.
#'
#' @name efftox_fit-class
#' @aliases efftox_fit
#' @docType class
#'
#' @details
#' See \code{methods(class = "efftox_fit")} for an overview of available
#' methods.
#'
#' @param dose_indices A vector of integers representing the dose-levels under
#' consideration.
#' @param recommended_dose An integer representing the dose-level recommended
#' for the next patient or cohort; or \code{NA} if stopping is recommended.
#' @param prob_eff The posterior mean probabilities of efficacy at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param prob_tox The posterior mean probabilities of toxicity at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @param prob_acc_eff The posterior mean probabilities that efficacy at the
#' doses is acceptable, i.e. that it exceeds the minimum acceptable efficacy
#' threshold; a vector of numbers between 0 and 1.
#' @param prob_acc_tox The posterior mean probabilities that toxicity at the
#' doses is acceptable, i.e. that it is less than the maximum toxicity
#' threshold; a vector of numbers between 0 and 1.
#' @param utility The utilities of doses 1:n, calculated by plugging the
#' posterior mean probabilities of efficacy and toxicity into the utility
#' formula, as advocated by Thall & Cook. Contrast to \code{post_utility};
#' a vector of numbers.
#' @param post_utility The posterior mean utilities of doses 1:n, calculated
#' from the posterior distributions of the utilities. This is in contrast to
#' \code{utility}, which uses plug-in posterior means of efficacy and toxicity,
#' as advocated by Thall & Cook; a vector of numbers.
#' @param acceptable A vector of logical values to indicate whether doses 1:n
#' are acceptable, according to the rules for acceptable efficacy & toxicity,
#' and rules on not skipping untested doses.
#' @param fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
#' @param dat Object \code{\link{efftox_params}} containing data passed to
#' \code{\link[rstan:sampling]{sampling}}.
#'
#' @export
#'
#' @seealso
#' \code{\link{stan_efftox}}
#' \code{\link{stan_efftox_demo}}
efftox_fit <- function(dose_indices, recommended_dose, prob_eff, prob_tox,
                       prob_acc_eff, prob_acc_tox, utility, post_utility,
                       acceptable, dat, fit) {
  # efftox_fit class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )
  x <- loo::nlist(dose_indices, recommended_dose, prob_eff, prob_tox,
                  prob_acc_eff, prob_acc_tox, utility, post_utility, acceptable,
                  dat, fit, version)
  class(x) <- "efftox_fit"
  x
}

#' Fit an EffTox model
#'
#' Fit an EffTox model using Stan for full Bayesian inference.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{efftox_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given}, \code{eff} and
#' \code{tox} parameters. See Details.
#' @param real_doses A vector of numbers.The doses under investigation. They
#' should be ordered from lowest to highest and be in consistent units.
#' E.g., #' to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
#' c(10, 20, 50).
#' @param efficacy_hurdle Minimum acceptable efficacy probability.
#' A number between 0 and 1.
#' @param toxicity_hurdle Maximum acceptable toxicity probability.
#' A number between 0 and 1.
#' @param p_e Certainty required to infer a dose is acceptable with regards to
#' being probably efficacious; a number between 0 and 1.
#' @param p_t Certainty required to infer a dose is acceptable with regards to
#' being probably tolerable; a number between 0 and 1.
#' @param eff0 Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1 (see Details).
#' @param tox1 Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1 (see Details).
#' @param eff_star Efficacy probability of an equi-utility third point (see
#' Details).
#' @param tox_star Toxicity probability of an equi-utility third point (see
#' Details).
#' @param alpha_mean The prior normal mean of the intercept term in the toxicity
#' logit model. A number.
#' @param alpha_sd The prior normal standard deviation of the intercept term in
#' the toxicity logit model. A number.
#' @param beta_mean The prior normal mean of the slope term in the toxicity
#' logit model. A number.
#' @param beta_sd The prior normal standard deviation of the slope term in the
#' toxicity logit model. A number.
#' @param gamma_mean The prior normal mean of the intercept term in the efficacy
#' logit model. A number.
#' @param gamma_sd The prior normal standard deviation of the intercept term in
#' the efficacy logit model. A number.
#' @param zeta_mean The prior normal mean of the slope term in the efficacy logit
#' model. A number.
#' @param zeta_sd The prior normal standard deviation of the slope term in the
#' efficacy logit model. A number.
#' @param eta_mean The prior normal mean of the squared term coefficient in the
#' efficacy logit model. A number.
#' @param eta_sd The prior normal standard deviation of the squared term
#' coefficient in the efficacy logit model. A number.
#' @param psi_mean The prior normal mean of the association term in the combined
#' efficacy-toxicity model. A number.
#' @param psi_sd The prior normal standard deviation of the association term in
#' the combined efficacy-toxicity model. A number.
#' @param doses_given A optional vector of dose-levels given to patients
#' 1:num_patients, where 1=lowest dose, 2=second dose, etc. Only required when
#' \code{outcome_str} is not provided.
#' @param eff An optional vector of efficacy outcomes for patients
#' 1:num_patients, where 1=efficacy and 0=no efficacy. Only required when
#' \code{outcome_str} is not provided.
#' @param tox An optional vector of toxicity outcomes for patients
#' 1:num_patients, where 1=toxicity and 0=no toxicity. Only required when
#' \code{outcome_str} is not provided.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. \code{\link[rstan:sampling]{sampling}}.
#'
#' @details
#' The quickest and easiest way to fit an EffTox model to some observed outcomes
#' is to describe the outcomes using \pkg{trialr}'s syntax for efficacy-toxicity
#' dose-finding outcomes. See \code{\link{efftox_parse_outcomes}} for full
#' details and examples.
#'
#' Utility or attractivess scores are calculated in EffTox using L^p norms.
#' Imagine the first quadrant of a scatter plot with prob_eff along the x-axis
#' and prob_tox along the y-axis.
#' The point (1, 0) (i.e. guaranteed efficacy & no toxicity) is the holy grail.
#' The neutral contour intersects the points (eff0, 0), (1, tox1) and
#' (eff_star, tox_star). A unique curve intersects these three points and
#' identifies a value for p, the exponent in the L^p norm. On this neutral-
#' utility contour, scores are equal to zero. A family of curves with different
#' utility scores is defined that are "parallel" to this neutral curve.
#' Points with probabilities of efficacy and toxicity that are nearer to (1, 0)
#' will yield greater scores, and vice-versa.
#'
#' @return An object of class \code{\link{efftox_fit}}
#'
#' @author Kristian Brock \email{kristian.brock@@gmail.com}
#'
#' @references
#'   Thall, P., & Cook, J. (2004). Dose-Finding Based on Efficacy-Toxicity
#'     Trade-Offs. Biometrics, 60(3), 684-693.
#'
#'   Thall, P., Herrick, R., Nguyen, H., Venier, J., & Norris, J. (2014).
#'     Effective sample size for computing prior hyperparameters in Bayesian
#'     phase I-II dose-finding. Clinical Trials, 11(6), 657-666.
#'     https://doi.org/10.1177/1740774514547397
#'
#'   Brock, K., Billingham, L., Copland, M., Siddique, S., Sirovica, M., &
#'     Yap, C. (2017). Implementing the EffTox dose-finding design in the
#'     Matchpoint trial. BMC Medical Research Methodology, 17(1), 112.
#'     https://doi.org/10.1186/s12874-017-0381-x
#'
#' @seealso
#'   \code{\link{efftox_fit}}
#'   \code{\link{stan_efftox_demo}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This model is presented in Thall et al. (2014)
#' mod1 <- stan_efftox('1N 2E 3B',
#'                      real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'                      efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
#'                      p_e = 0.1, p_t = 0.1,
#'                      eff0 = 0.5, tox1 = 0.65,
#'                      eff_star = 0.7, tox_star = 0.25,
#'                      alpha_mean = -7.9593, alpha_sd = 3.5487,
#'                      beta_mean = 1.5482, beta_sd = 3.5018,
#'                      gamma_mean = 0.7367, gamma_sd = 2.5423,
#'                      zeta_mean = 3.4181, zeta_sd = 2.4406,
#'                      eta_mean = 0, eta_sd = 0.2,
#'                      psi_mean = 0, psi_sd = 1, seed = 123)
#'
#' # Shorthand for the above is:
#' mod2 <- stan_efftox_demo('1N 2E 3B', seed = 123)
#'
#' # the seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
#' }
stan_efftox <- function(outcome_str = NULL,
                        real_doses, efficacy_hurdle, toxicity_hurdle, p_e, p_t,
                        eff0, tox1, eff_star, tox_star,
                        alpha_mean, alpha_sd, beta_mean, beta_sd,
                        gamma_mean, gamma_sd, zeta_mean, zeta_sd,
                        eta_mean, eta_sd, psi_mean, psi_sd,
                        doses_given = NULL,
                        eff = NULL,
                        tox = NULL,
                        ...) {

  # Create parameters object to pass to Stan
  dat <- efftox_params(real_doses, efficacy_hurdle, toxicity_hurdle,
                       p_e, p_t, eff0, tox1, eff_star, tox_star,
                       alpha_mean, alpha_sd, beta_mean, beta_sd,
                       gamma_mean, gamma_sd, zeta_mean, zeta_sd,
                       eta_mean, eta_sd, psi_mean, psi_sd)

  # Add outcomes
  if(is.null(outcome_str)) {
    if(length(doses_given) != length(eff))
      stop('doses_given and eff vectors should have same length')
    if(length(tox) != length(eff))
      stop('tox and eff vectors should have same length')
    dat$doses <- doses_given
    dat$eff <- eff
    dat$tox <- tox
    dat$num_patients <- length(doses_given)
  } else {
    outcomes_df <- efftox_parse_outcomes(outcome_str, as.list = TRUE)
    dat$num_patients <- outcomes_df$num_patients
    dat$doses <- outcomes_df$doses
    dat$eff <- outcomes_df$eff
    dat$tox <- outcomes_df$tox
  }

  # Fit data to model using Stan
  samp <- rstan::sampling(stanmodels$EffTox, data = dat, ...)
  # Create useful output from posterior samples
  decision <- efftox_process(dat, samp)

  return(decision)
}

#' Fit the EffTox model presented in Thall et al. (2014)
#'
#' Fit the EffTox model presented in Thall et al. (2014) using Stan for full
#' Bayesian inference.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{efftox_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given}, \code{eff} and
#' \code{tox} parameters. See Details.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. \code{\link[rstan:sampling]{sampling}}.
#'
#' @return An object of class \code{\link{efftox_fit}}
#'
#' @author Kristian Brock \email{kristian.brock@@gmail.com}
#'
#' @references
#'   Thall, P., & Cook, J. (2004). Dose-Finding Based on Efficacy-Toxicity
#'     Trade-Offs. Biometrics, 60(3), 684-693.
#'
#'   Thall, P., Herrick, R., Nguyen, H., Venier, J., & Norris, J. (2014).
#'     Effective sample size for computing prior hyperparameters in Bayesian
#'     phase I-II dose-finding. Clinical Trials, 11(6), 657-666.
#'     https://doi.org/10.1177/1740774514547397
#'
#'   Brock, K., Billingham, L., Copland, M., Siddique, S., Sirovica, M., &
#'     Yap, C. (2017). Implementing the EffTox dose-finding design in the
#'     Matchpoint trial. BMC Medical Research Methodology, 17(1), 112.
#'     https://doi.org/10.1186/s12874-017-0381-x
#'
#' @seealso
#'   \code{\link{efftox_fit}}
#'   \code{\link{stan_efftox}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This model is presented in Thall et al. (2014)
#' mod2 <- stan_efftox_demo('1N 2E 3B', seed = 123)
#'
#' # The seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
#' }
stan_efftox_demo <- function(outcome_str, ...) {
  stan_efftox(outcome_str,
              real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
              efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
              p_e = 0.1, p_t = 0.1,
              eff0 = 0.5, tox1 = 0.65,
              eff_star = 0.7, tox_star = 0.25,
              alpha_mean = -7.9593, alpha_sd = 3.5487,
              beta_mean = 1.5482, beta_sd = 3.5018,
              gamma_mean = 0.7367, gamma_sd = 2.5423,
              zeta_mean = 3.4181, zeta_sd = 2.4406,
              eta_mean = 0, eta_sd = 0.2,
              psi_mean = 0, psi_sd = 1, ...)
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
  treated <- data.frame(
    Patient = 1:length(x$dat$doses),
    Dose = x$dat$doses,
    Toxicity = x$dat$tox,
    Efficacy = x$dat$eff
  )
  print(treated)
  cat('\n')

  # Dose-level data
  df <- efftox_analysis_to_df(x)
  print(df)
  cat('\n')

  # Extras
  if(sum(x$acceptable) == 0) {
    cat('The model advocates stopping.')

  } else {
    cat(paste0('The model recommends selecting dose-level ',
               x$recommended_dose, '.'))
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
