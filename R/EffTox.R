


# EffTox

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
#' @export
#'
#' @examples
#' efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' @references Thall et al. 2014, Effective sample size for computing prior
#' hyperparameters in Bayesian phase I-II dose-finding
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

#' @title Get parameters to run the EffTox demo
#'
#' @description Get parameters to run the EffTox demo. These match those used
#' to demonstrate EffTox in Thall et al. 2014.
#'
#' @return a \code{list} of parameters
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' names(dat)
#' dat$real_doses == c(1, 2, 4, 6.6, 10)
#'
#' @seealso \link{\code{efftox_params}}
#' @references Thall et al. 2014, Effective sample size for computing prior
#' hyperparameters in Bayesian phase I-II dose-finding
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
#' @export
#'
#' @examples
#' p <- efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' u <- efftox_utility(p, 0.5, 0.65, prob_eff = 0.7, prob_tox = 0.25)
#' round(u, 4) == 0
#'
#' u <- efftox_utility(p, 0.5, 0.65, prob_eff = c(0.6, 0.7, 0.8), prob_tox = c(0.1, 0.2, 0.3))
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
#' @description Process RStan samples from an EffTox model to make inferences about
#' dose-acceptability, dose-utility and which dose should be recommended next.
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param fit An instance of \code{rstan::stanmodel}, derived by sampling an
#' EffTox model. Use \code{stan::sampling(stanmodels$EffTox, data = dat)}.
#' @param p_e Certainty required to infer a dose is acceptable with regards to
#' being probably efficacious; a number between 0 and 1.
#' @param p_t Certainty required to infer a dose is acceptable with regards to
#' being probably tolerable; a number between 0 and 1.
#' @return An instance of \code{\link{efftox_analysis}}.
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' dat$num_patients <- 3
#' dat$eff <- c(0, 1, 1)
#' dat$tox <- c(0, 0, 1)
#' dat$doses <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = dat)
#' decision <- efftox_process(dat, fit, p_e = 0.1, p_t = 0.1)
#' decision$recommended_dose == 3
#' @seealso
#' \code{\link{efftox_params}}
#'
#' \code{\link{efftox_parameters_demo}}
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


#' @title EffTox analysis to data.frame
#'
#' @description Convenient function to turn an \code{\link{efftox_analysis}}
#' into a \code{data.frame}.
#'
#' @param x An \code{\link{efftox_analysis}}
#'
#' @return a \code{data.frame}
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' dat$num_patients <- 3
#' dat$eff <- c(0, 1, 1)
#' dat$tox <- c(0, 0, 1)
#' dat$doses <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = dat)
#' decision <- efftox_process(dat, fit, p_e = 0.1, p_t = 0.1)
#' df = efftox_analysis_to_df(decision)
#' round(df$Utility, 2) == c(-0.64, 0.04, 0.24, -0.05, -0.19)
#'
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
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param num_sims integer, number of simulated iterations
#' @param first_dose integer, the dose-level to give to patient 1, e.g. 1 for
#' the lowest dose.
#' @param p_e Certainty required to infer a dose is acceptable with regards to
#' being probably efficacious; a number between 0 and 1.
#' @param p_t Certainty required to infer a dose is acceptable with regards to
#' being probably tolerable; a number between 0 and 1.
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
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' set.seed(123)
#' # Let's say we want to use only 2 chains. Extra args are passed to stan
#' sims = efftox_simulate(dat, num_sims = 2, first_dose = 1, p_e, p_t,
#'                        true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
#'                        true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
#'                        cohort_sizes = rep(3, 13),
#'                        chains = 2)
#' table(sims$recommended_dose) / length(sims$recommended_dose)
#' table(unlist(sims$doses_given)) / length(unlist(sims$doses_given))
#' table(unlist(sims$doses_given)) / length(sims$recommended_dose)
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

#' @title Get the probability of toxicity for probability-of-efficacy and utility pairs
#'
#' @description Get the probability of toxicity for probability-of-efficacy and utility pairs
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
#' @export
#'
#' @examples
#' p <- efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' prob_tox <- efftox_get_tox(0.7, 0, p, eff0 = 0.5, tox1 = 0.65)
#' round(prob_tox, 2) == 0.25
#'
#' prob_tox <- efftox_get_tox(0.7, seq(-0.5, 0.25, by = 0.25), p, eff0 = 0.5, tox1 = 0.65)
#' round(prob_tox, 2) == c(0.57, 0.41, 0.25, 0.09)
#'
#' prob_tox <- efftox_get_tox(c(0.5, 0.7, 0.8), 0.25, p, eff0 = 0.5, tox1 = 0.65)
#' round(prob_tox, 2) == c(NaN, 0.09, 0.22)
#'
#' prob_tox <- efftox_get_tox(c(0.5, 0.7, 0.8), c(-1, 0, 1), p, eff0 = 0.5, tox1 = 0.65)
#' round(prob_tox, 2) == c(0.63, 0.25, NaN)
#'
#' @note Various ways of vectorising the function are demonstrated in the examples
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
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param use_ggplot logical, TRUE to use ggplot2. Defaults to FALSE to use
#'  standard R graphics.
#' @param prob_eff an optional vector of numbers between 0 and 1, containing the
#' efficacy probabilities of extra points to add to the plot as points,
#' e.g. the posterior mean efficacy probabilities of the doses under
#' investigation. Paired with prob_tox, thus they should be the same length.
#' @param prob_tox an optional vector of numbers between 0 and 1, containing the
#' toxicity probabilities of extra points to add to the plot as points,
#' e.g. the posterior mean toxicity probabilities of the doses under
#' investigation. Paired with prob_eff, thus they should be the same length.
#' @param num_points integer for number of points to calculate on each curve.
#' The default is 1000 and this should be plenty.
#' @param util_vals A contour is plotted for each of these utility values.
#' The default is contours spaced by 0.2 between from -3 and 3,
#' i.e. \code{seq(-3, 3, by = 0.2)}.
#'
#' @return if \code{use_ggplot = TRUE}, an instance of \code{ggplot}; else no
#' object is returned. Omit assignment in either case to just view the plot.
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' efftox_contour_plot(dat)
#' # Add posterior beliefs
#' dat$num_patients <- 3
#' dat$eff <- c(0, 1, 1)
#' dat$tox <- c(0, 0, 1)
#' dat$doses <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = dat)
#' decision <- efftox_process(dat, fit, p_e = 0.1, p_t = 0.1)
#' efftox_contour_plot(dat, prob_eff = decision$prob_eff, prob_tox = decision$prob_tox)
#' title('EffTox utility contours')
#' # The same with ggplot2
#' efftox_contour_plot(dat, prob_eff = decision$prob_eff, prob_tox = decision$prob_tox,
#' use_ggplot = TRUE) + ggtitle('EffTox utility contours')
#'
#' @seealso
#' \code{\link{efftox_params}}
#'
#' \code{\link{efftox_parameters_demo}}
efftox_contour_plot <- function(dat,
                                use_ggplot = FALSE,
                                prob_eff = NULL, prob_tox = NULL,
                                num_points = 1000,
                                util_vals = seq(-3, 3, by = 0.2)) {
  eff_vals = seq(0, 1, length.out = num_points)

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
      ggplot2::xlab('Prob(Efficacy)') + ggplot2::ylab('Prob(Toxicity)')

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
      # prob_eff = decision$prob_eff
      # prob_tox = decision$prob_tox
      df4 <- data.frame(prob_eff, prob_tox, dl = 1:length(prob_eff))
      plt <- plt + ggplot2::geom_text(data = df4, ggplot2::aes(x = prob_eff,
                                                               y = prob_tox,
                                                               group = 1,
                                                               label = dl),
                                      col = 'red', size = 4)
    }
    return(plt)
  } else {
    plot(NULL, ylim = c(0, 1), xlim = c(0, 1), ylab = 'Prob(Toxicity)',
         xlab = 'Prob(Efficacy)')

    for(u in util_vals) {
      tox_vals = efftox_get_tox(eff_vals, u, dat$p, dat$eff0, dat$tox1)
      points(eff_vals, tox_vals, type = 'l', col = 'grey', lwd = 0.2)
    }

    # Add neutral utility contour
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


#' @title Plot densities of EffTox dose utilities
#'
#' @description Plot densities of EffTox dose utilities. Optionally plot only a
#' subset of the doses by specifying the \code{doses} parameter. This function
#' requires ggplot2 be installed.
#'
#' @param fit An instance of \code{rstan::stanmodel}, derived by sampling an
#' EffTox model. Use \code{stan::sampling(stanmodels$EffTox, data = dat)}.
#' @param doses optional, vector of integer dose-levels to plot. E.g. to plot
#' only dose-levels 1, 2 & 3 (and suppress the plotting of any other doses), use
#' \code{doses = 1:3}
#'
#' @return an instance of \code{ggplot}. Omit assignment to just view the plot.
#' @export
#'
#' @note This function requires that ggplot2 be installed.
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' dat$num_patients <- 3
#' dat$eff <- c(0, 1, 1)
#' dat$tox <- c(0, 0, 1)
#' dat$doses <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = dat)
#' efftox_utility_density_plot(fit) + ggtitle('My doses')  # Bit too busy?
#' efftox_utility_density_plot(fit, doses = 1:3) + ggtitle('My doses') # Clearer
#'
#' @seealso
efftox_utility_density_plot <- function(fit, doses = NULL) {
  if(!('ggplot2' %in% installed.packages()))
    stop('THis function requires ggplot2 be installed.')

  u <- rstan::extract(fit, par = 'utility')[[1]]
  df <- data.frame(Utility = as.numeric(u),
                   D = rep(1:5, each = nrow(u))
  )
  df$Dose = factor(df$D)
  if(!is.null(doses))
    df = df[df$D %in% doses, ]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Utility, group = Dose,
                                        colour = Dose)) + ggplot2::geom_density()
  return(p)
}


#' @title Get dose-superiority matrix in EffTox
#'
#' @description Get a dose-superiority matrix from an EffTox dose analysis.
#' EffTox seeks to choose the dose with the highest utility, thus superiority
#' is inferred by posterior utility. The item in row i, col j is the posterior
#' probability that the utility of dose j exceeds that of dose i.
#'
#' @param fit An instance of \code{rstan::stanmodel}, derived by sampling an
#' EffTox model. Use \code{stan::sampling(stanmodels$EffTox, data = dat)}.
#'
#' @return n by n matrix, where n is number of doses under investigation.
#' The item in row i, col j is the posterior probability that the utility of
#' dose j exceeds that of dose i.
#' @export
#'
#' @examples
#' dat <- efftox_parameters_demo()
#' dat$num_patients <- 3
#' dat$eff <- c(0, 1, 1)
#' dat$tox <- c(0, 0, 1)
#' dat$doses <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = dat)
#' sup_mat <- efftox_superiority(fit)
efftox_superiority <- function(fit) {
  u <- rstan::extract(fit, par = 'utility')[[1]]
  superiority_mat <- sapply(1:ncol(u), function(i) sapply(1:ncol(u), function(j)
    mean(u[ , i] > u[ , j])))
  diag(superiority_mat) <- NA
  dimnames(superiority_mat) = list(1:ncol(u), 1:ncol(u))
  return(superiority_mat)
}
