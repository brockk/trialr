
library(trialr)
library(dplyr)
library(loo)

# Code ----
#' Class of models fit by \pkg{trialr} using the EffTox design.
#'
#' @name efftox_fit-class
#' @aliases efftox_fit
#' @docType class
#'
#' @details
#' See \code{methods(class = "efftox_fit")} for an overview of available
#' methods.
#'
#' @slot
#'
#' @slot dose_indices A vector of integers representing the dose-levels under
#' consideration.
#' @slot recommended_dose An integer representing the dose-level recommended
#' for the next patient or cohort; or \code{NA} stopping is recommended.
#' @slot prob_eff The posterior mean probabilities of efficacy at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @slot prob_tox The posterior mean probabilities of toxicity at doses 1:n;
#' a vector of numbers between 0 and 1.
#' @slot prob_acc_eff The posterior mean probabilities that efficacy at the
#' doses is acceptable, i.e. that it exceeds the minimum acceptable efficacy
#' threshold; a vector of numbers between 0 and 1.
#' @slot prob_acc_tox The posterior mean probabilities that toxicity at the
#' doses is acceptable, i.e. that it is less than the maximum toxicity
#' threshold; a vector of numbers between 0 and 1.
#' @slot utility The utilities of doses 1:n, calculated by plugging the
#' posterior mean probabilities of efficacy and toxicity into the utility
#' formula, as advocated by Thall & Cook. Contrast to \code{post_utility};
#' a vector of numbers.
#' @slot post_utility The posterior mean utilities of doses 1:n, calculated
#' from the posterior distributions of the utilities. This is in contrast to
#' \code{utility}, which uses plug-in posterior means of efficacy and toxicity,
#' as advocated by Thall & Cook; a vector of numbers.
#' @slot acceptable A vector of logical values to indicate whether doses 1:n
#' are acceptable, according to the rules for acceptable efficacy & toxicity,
#' and rules on not skipping untested doses.
#' @slot fit An object of class \code{\link[rstan:stanfit]{stanfit}},
#' containing the posterior samples.
#'
#' @seealso
#' \code{\link{stan_efftox}}
#' \code{\link{stan_efftox_demo}}
#' \code{\link{efftox_process}}
efftox_fit <- function(dose_indices, recommended_dose, prob_eff, prob_tox,
                       prob_acc_eff, prob_acc_tox, utility, post_utility,
                       acceptable, fit) {
  # efftox_fit class
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )
  x <- nlist(dose_indices, recommended_dose, prob_eff, prob_tox, prob_acc_eff,
             prob_acc_tox, utility, post_utility, acceptable, fit, version)
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
#' @param gamma_mean The prior mean of the intercept term in the efficacy logit model. A number.
#' @param gamma_sd The prior standard deviation of the intercept term in the efficacy logit model. A number.
#' @param zeta_mean The prior mean of the slope term in the efficacy logit model. A number.
#' @param zeta_sd The prior standard deviation of the slope term in the efficacy logit model. A number.
#' @param eta_mean The prior mean of the squared term coefficient in the efficacy logit model. A number.
#' @param eta_sd The prior standard deviation of the squared term coefficient in the efficacy logit model. A number.
#' @param psi_mean The prior mean of the association term in the combined efficacy-toxicity model. A number.
#' @param psi_sd The prior standard deviation of the association term in the combined efficacy-toxicity model. A number.
#' @param doses_given A optional vector of dose-levels given to patients
#' 1:num_patients, where 1=lowest dose, 2=second dose, etc. Only required when
#' \code{outcome_str} is not provided.
#' @param eff An optional vector of efficacy outcomes for patients
#' 1:num_patients, where 1=efficacy and 0=no efficacy. Only required when
#' \code{outcome_str} is not provided.
#' @param tox An optional vector of toxicity outcomes for patients
#' 1:num_patients, where 1=toxicity and 0=no toxicity. Only required when
#' \code{outcome_str} is not provided.
#' @param ...Extra parameters are passed to \code{rstan::sampling}. Commonly
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
#'     Trade-Offs. Biometrics, 60(3), 684–693.
#'
#'   Thall, P., Herrick, R., Nguyen, H., Venier, J., & Norris, J. (2014).
#'     Effective sample size for computing prior hyperparameters in Bayesian
#'     phase I-II dose-finding. Clinical Trials, 11(6), 657–666.
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
#'   \code{\link{efftox_process}}
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

  p <- efftox_solve_p(eff0, tox1, eff_star, tox_star)
  # Create data object to pass to Stan. Add parameters
  dat <- nlist(real_doses, num_doses = length(real_doses),
    efficacy_hurdle, toxicity_hurdle,
    p_e, p_t, p, eff0, tox1, eff_star, tox_star,
    alpha_mean, alpha_sd, beta_mean, beta_sd, gamma_mean, gamma_sd,
    zeta_mean, zeta_sd, eta_mean, eta_sd, psi_mean, psi_sd
  )

  # Add outcomes
  if(is.null(outcome_str)) {
    if(length(doses_given) != length(efficacy))
      stop('doses_given and efficacy vectors should have same length')
    if(length(toxicity) != length(efficacy))
      stop('toxicity and efficacy vectors should have same length')
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

#' Fit the EffTox model presented in Thal et al. (2014)
#'
#' Fit the EffTox model presented in Thal et al. (2014) using Stan for full
#' Bayesian inference.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{efftox_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given}, \code{eff} and
#' \code{tox} parameters. See Details.
#' @param ...Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores},
#' \code{control}. \code{\link[rstan:sampling]{sampling}}.
#'
#' @return An object of class \code{\link{efftox_fit}}
#'
#' @author Kristian Brock \email{kristian.brock@@gmail.com}
#'
#' @references
#'   Thall, P., & Cook, J. (2004). Dose-Finding Based on Efficacy-Toxicity
#'     Trade-Offs. Biometrics, 60(3), 684–693.
#'
#'   Thall, P., Herrick, R., Nguyen, H., Venier, J., & Norris, J. (2014).
#'     Effective sample size for computing prior hyperparameters in Bayesian
#'     phase I-II dose-finding. Clinical Trials, 11(6), 657–666.
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
#'   \code{\link{efftox_process}}
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
              p = 0.9773632,
              eff0 = 0.5, tox1 = 0.65,
              eff_star = 0.7, tox_star = 0.25,
              alpha_mean = -7.9593, alpha_sd = 3.5487,
              beta_mean = 1.5482, beta_sd = 3.5018,
              gamma_mean = 0.7367, gamma_sd = 2.5423,
              zeta_mean = 3.4181, zeta_sd = 2.4406,
              eta_mean = 0, eta_sd = 0.2,
              psi_mean = 0, psi_sd = 1, ...)
}

print.efftox_fit <- function(x) {
  df <- efftox_analysis_to_df(x)
  print(df)
  if(sum(x$acceptable) == 0) {
    cat('The model advocates stopping.')

  } else {
    cat(paste0('The model recommends selecting dose-level ',
               x$recommended_dose, '.'))
  }
}

as.data.frame.efftox_fit <- function(x, ...) {
  as.data.frame(x$fit, ...)
}

plot.efftox_fit <- function(x,  pars = 'utility', ...) {
  plot(x$fit, pars = pars, ...)
}



# Usage ----
mod1 <- stan_efftox_demo('1N 2E 3B', seed = 123)
names(mod1)
print(mod1)
as.data.frame(mod1) %>% head
as.data.frame(mod1, pars = c('utility')) %>% head
plot(mod1)
plot(mod1, pars = 'prob_eff')

mod2 <- stan_efftox('1N 2E 3B',
            real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
            efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
            p_e = 0.1, p_t = 0.1,
            p = 0.9773632, eff0 = 0.5, tox1 = 0.65,
            eff_star = 0.7, tox_star = 0.25,
            alpha_mean = -7.9593, alpha_sd = 3.5487,
            beta_mean = 1.5482, beta_sd = 3.5018,
            gamma_mean = 0.7367, gamma_sd = 2.5423,
            zeta_mean = 3.4181, zeta_sd = 2.4406,
            eta_mean = 0, eta_sd = 0.2,
            psi_mean = 0, psi_sd = 1, seed = 123)
# chains = 6, iter = 5000
mod2

dat <- efftox_parameters_demo()
dat$num_patients <- 3
dat$eff <- c(0, 1, 1)
dat$tox <- c(0, 0, 1)
dat$doses <- c(1, 2, 3)
samp3 <- rstan::sampling(stanmodels$EffTox, data = dat, seed = 123)
mod3 <- efftox_process(dat, samp3)
mod3

mod1
mod2
mod3

mod4 <- stan_efftox_demo('1NNN 2ENN 3BTE', chains = 6, iter = 5000)
mod4

stan_efftox_demo('1TT')
stan_efftox_demo('1TT 2TT')
stan_efftox_demo('1TT 2TT 3TT')
stan_efftox_demo('1TT 2TT 3TT 4TT') # Stops, at last
