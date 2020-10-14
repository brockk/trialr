
#' @importFrom gtools inv.logit
.pi_T_sample <- function(scaled_doses, n = 10000,
                         alpha_mean, alpha_sd,
                         beta_mean, beta_sd) {

  alpha <- rnorm(n, alpha_mean, alpha_sd)
  beta <- rnorm(n, beta_mean, beta_sd)
  A <- as.matrix(data.frame(alpha, beta))
  B <- as.matrix(data.frame(1, scaled_doses))
  inv.logit(A %*% t(B))
}

#' @importFrom stats sd
.pi_T_moments <- function(scaled_doses, n = 10000,
                          alpha_mean, alpha_sd,
                          beta_mean, beta_sd) {
  samp <- .pi_T_sample(scaled_doses, n, alpha_mean, alpha_sd, beta_mean,
                       beta_sd)
  list(
    mean = colMeans(samp),
    sd = apply(samp, 2, sd)
  )
}

#' @importFrom gtools inv.logit
.pi_E_sample <- function(scaled_doses, n = 10000,
                         gamma_mean, gamma_sd,
                         zeta_mean, zeta_sd,
                         eta_mean, eta_sd) {

  gamma <- rnorm(n, gamma_mean, gamma_sd)
  zeta <- rnorm(n, zeta_mean, zeta_sd)
  eta <- rnorm(n, eta_mean, eta_sd)
  A <- as.matrix(data.frame(gamma, zeta, eta))
  B <- as.matrix(data.frame(1, scaled_doses, scaled_doses^2))
  inv.logit(A %*% t(B))
}

#' @importFrom stats sd
.pi_E_moments <- function(scaled_doses, n = 10000,
                          gamma_mean, gamma_sd,
                          zeta_mean, zeta_sd,
                          eta_mean, eta_sd) {
  samp <- .pi_E_sample(scaled_doses, n, gamma_mean, gamma_sd, zeta_mean,
                       zeta_sd, eta_mean, eta_sd)
  list(mean = colMeans(samp), sd = apply(samp, 2, sd))
}

.normal_ess <- function(mean, sd) {
  a <- (1 - mean) * (mean / sd)^2 - mean
  b <- a * (1 - mean) / mean
  list(a = a, b = b, ess = a + b)
}

.ess <- function(pi_star) {
  x <- sapply(
    1:length(pi_star$mean),
    function(i) .normal_ess(pi_star$mean[i], pi_star$sd[i])$ess
  )
  mean(x)
}

#' Get normal prior hyperparameters for the EffTox model.
#'
#' Get normal prior hyperparameters for the EffTox model using the algorithm
#' presented in Thall et al. (2014) that targets a family of priors with a
#' pre-specified effective sample size (ESS).
#'
#' @param doses A vector of numbers, the doses under investigation. They
#' should be ordered from lowest to highest and be in consistent units.
#' E.g. to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
#' c(10, 20, 50). Specify \code{doses} or \code{scaled_doses}.
#' @param scaled_doses Optional, vector of numbers, representing the scaled
#' doses under investigation. Thall et al. advocate
#' \code{scaled_doses = log(doses) - mean(log(doses))}, and that is what we use
#' here. Specify \code{doses} or \code{scaled_doses}.
#' @param pi_T Vector of prior expectations of probabilities of toxicity at the
#' doses. Should be congruent to \code{doses} or \code{scaled_doses}.
#' @param ess_T Numerical, sought total effective sample size for priors on
#' parameters in the toxicity sub-model. Thall et al. (2014) advocate values in
#' (0.3, 1.0) but stress that stress-testing with values outside this range may
#' be necessary.
#' @param pi_E Vector of prior expectations of probabilities of efficacy at the
#' doses. Should be congruent to \code{doses} or \code{scaled_doses}.
#' @param ess_E Numerical, sought total effective sample size for priors on
#' parameters in the efficacy sub-model. Thall et al. (2014) advocate values in
#' (0.3, 1.0) but stress that stress-testing with values outside this range may
#' be necessary.
#' @param num_samples Number of samples to draw from priors. The default 10^4
#' seems to be a nice compromise between accuracy and speed. Orders of magnitude
#' larger take a long time to run.
#' @param seed Optional seed. This process involves randomness so seeds are used
#' for repeatable results.
#'
#' @return An instance of class \code{\link{efftox_priors}}.
#'
#' @importFrom stats optim
#'
#' @export
#'
#' @references
#'   Thall, P., Herrick, R., Nguyen, H., Venier, J., & Norris, J. (2014).
#'     Effective sample size for computing prior hyperparameters in Bayesian
#'     phase I-II dose-finding. Clinical Trials, 11(6), 657-666.
#'     https://doi.org/10.1177/1740774514547397
#'
#' @examples
#' \dontrun{
#' # Reproduce the priors calculated in Thall et al. (2014)
#' p <- get_efftox_priors(
#'   doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'   pi_T = c(0.02, 0.04, 0.06, 0.08, 0.10), ess_T = 0.9,
#'   pi_E = c(0.2, 0.4, 0.6, 0.8, 0.9), ess_E = 0.9
#' )
#' p
#' # These are close to the published example. They do not match exactly because
#' # the process of deriving them is iterative.
#' }
get_efftox_priors <- function(doses = NULL, scaled_doses = NULL,
                              pi_T, ess_T, pi_E, ess_E, num_samples = 10^4,
                              seed = 123) {
  if(is.null(scaled_doses)) {
    if(is.null(doses)) {
      stop('Either doses or scaled_doses should be a numerical vector.')
    } else {
      y <- log(doses) - mean(log(doses))
    }
  } else {
    y <- scaled_doses
  }

  .f <- function(x) {
    set.seed(seed)
    pi_T_star <- .pi_T_moments(scaled_doses = y, n = num_samples,
                               alpha_mean = x[1], alpha_sd = x[2],
                               beta_mean = x[3], beta_sd = x[4])
    ess_T_star <- .ess(pi_T_star)
    pi_E_star <- .pi_E_moments(scaled_doses = y, n = num_samples,
                               gamma_mean = x[5], gamma_sd = x[6],
                               zeta_mean = x[7], zeta_sd = x[8],
                               eta_mean = 0, eta_sd = 0.2)
    ess_E_star <- .ess(pi_E_star)

    obj <- sum((pi_T_star$mean - pi_T)^2) + sum((pi_E_star$mean - pi_E)^2)
    obj <- obj + 0.1 * ((ess_T_star - ess_T)^2 + (ess_E_star - ess_E)^2)
    obj <- obj + 0.2 * ((x[2] - x[4])^2 + (x[6] - x[8])^2)
    obj
  }

  par <- c(0, 1, 0, 1, 0, 1, 0, 1)
  x <- optim(par, fn = .f, method = "Nelder-Mead")
  efftox_priors(alpha_mean = x$par[1], alpha_sd = x$par[2],
                beta_mean = x$par[3], beta_sd = x$par[4],
                gamma_mean = x$par[5], gamma_sd = x$par[6],
                zeta_mean = x$par[7], zeta_sd = x$par[8],
                eta_mean = 0, eta_sd = 0.2,
                psi_mean = 0, psi_sd = 1)
}
