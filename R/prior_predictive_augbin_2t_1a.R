
#' Sample data from the Augmented Binary model prior predictive distribution.
#'
#' Sample data from the prior predictive distributions of the two-period, single
#' arm Augmented Binary model, subject to chosen prior parameters.
#'
#' @param num_samps Number of samples.
#' @param alpha_mean Prior mean of alpha parameter.
#' @param alpha_sd Prior sd of alpha parameter.
#' @param beta_mean Prior mean of beta parameter.
#' @param beta_sd Prior sd of beta parameter.
#' @param gamma_mean Prior mean of gamma parameter.
#' @param gamma_sd Prior sd of gamma parameter.
#' @param sigma_mean Prior mean of sigma parameter.
#' @param sigma_sd Prior sd of sigma parameter.
#' @param omega_lkj_eta Prior eta parameter for LKJ prior on covariance matrix
#' of log tumour sizes.
#' @param alpha_d1_mean Prior mean of alpha_D1 parameter.
#' @param alpha_d1_sd Prior sd of alpha_D1 parameter.
#' @param gamma_d1_mean Prior mean of gamma_D1 parameter.
#' @param gamma_d1_sd Prior sd of gamma_D1 parameter.
#' @param alpha_d2_mean Prior mean of alpha_D2 parameter.
#' @param alpha_d2_sd Prior sd of alpha_D2 parameter.
#' @param gamma_d2_mean Prior mean of gamma_D2 parameter.
#' @param gamma_d2_sd Prior sd of gamma_D2 parameter.
#'
#' @return Object of class \code{\link[tibble]{tibble}}
#'
#' @export
#'
#' @seealso
#'   \code{\link{stan_augbin}}
#'
#' @importFrom tibble tibble
#' @importFrom gtools inv.logit
#' @importFrom stats runif rnorm
#'
#' @examples
#' prior_predictive_augbin_2t_1a(num_samps = 1000,
#'                               alpha_mean = 0, alpha_sd = 1,
#'                               beta_mean = 0, beta_sd = 1,
#'                               gamma_mean = 0, gamma_sd = 1,
#'                               sigma_mean = 0, sigma_sd = 1,
#'                               omega_lkj_eta = 1,
#'                               alpha_d1_mean = 0, alpha_d1_sd = 1,
#'                               gamma_d1_mean = 0, gamma_d1_sd = 1,
#'                               alpha_d2_mean = 0, alpha_d2_sd = 1,
#'                               gamma_d2_mean = 0, gamma_d2_sd = 1)
prior_predictive_augbin_2t_1a <- function(
  num_samps,
  alpha_mean, alpha_sd, beta_mean, beta_sd,
  gamma_mean, gamma_sd, sigma_mean, sigma_sd,
  omega_lkj_eta,
  alpha_d1_mean, alpha_d1_sd,
  gamma_d1_mean, gamma_d1_sd,
  alpha_d2_mean, alpha_d2_sd,
  gamma_d2_mean, gamma_d2_sd
) {
  N <- num_samps
  alpha_samp <- rnorm(n = N, mean = alpha_mean, sd = alpha_sd)
  beta_samp <- rnorm(n = N, mean = beta_mean, sd = beta_sd)
  gamma_samp <- rnorm(n = N, mean = gamma_mean, sd = gamma_sd)
  # LKJ
  Omega_samp <- rlkjcorr(n = N, K = 2, eta = omega_lkj_eta)
  sigma1_samp <- abs(rnorm(n = N, mean = sigma_mean, sd = sigma_sd))
  sigma2_samp <- abs(rnorm(n = N, mean = sigma_mean, sd = sigma_sd))
  alpha_D1_samp <- rnorm(n = N, mean = alpha_d1_mean, sd = alpha_d1_sd)
  gamma_D1_samp <- rnorm(n = N, mean = gamma_d1_mean, sd = gamma_d1_sd)
  alpha_D2_samp <- rnorm(n = N, mean = alpha_d2_mean, sd = alpha_d2_sd)
  gamma_D2_samp <- rnorm(n = N, mean = gamma_d2_mean, sd = gamma_d2_sd)

  z0_samp <- runif(n = N, min = 5, max = 10)
  mu1_samp <- alpha_samp + gamma_samp * z0_samp
  mu2_samp <- beta_samp + gamma_samp * z0_samp
  y_samp <- t(sapply(
    1:N,
    function(i) MASS::mvrnorm(n = 1, mu = c(mu1_samp[i], mu2_samp[i]),
                              Sigma = matrix(c(sigma1_samp[i], 0,
                                               0, sigma2_samp[i]), ncol = 2) %*%
                                Omega_samp[i, , ]  %*%
                                matrix(c(sigma1_samp[i], 0,
                                         0, sigma2_samp[i]), ncol = 2))
  ))

  z1_samp <- z0_samp * exp(y_samp[, 1])
  z2_samp <- z0_samp * exp(y_samp[, 2])
  prob_d1_samp <- inv.logit(alpha_D1_samp + gamma_D1_samp * z0_samp)
  prob_d2_samp <- inv.logit(alpha_D2_samp + gamma_D2_samp * z1_samp)
  tibble(
    id = 1:length(z0_samp),
    z0 = z0_samp,
    z1 = z1_samp,
    z2 = z2_samp,
    y0 = rep(0, length(z0_samp)),
    y1 = log(z1_samp / z0_samp),
    y2 = log(z2_samp / z0_samp),
    prob_d1 = prob_d1_samp,
    prob_d2 = prob_d2_samp
  )
}
