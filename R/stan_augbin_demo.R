
#' Simple helper function to demonstrate fitting of an Augmented Binary model.
#'
#' This function exist mostly to demonstrate things you can do to instances of
#' \code{\link{augbin_fit}} without having to paste into each example the not
#' inconsiderable blob of code to sample outcomes and fit the model.
#'
#' @return instance of \code{\link{augbin_fit}}
#'
#' @seealso
#'   \code{\link{stan_augbin}}
#'   \code{\link{augbin_fit}}
#'   \code{\link{prior_predictive_augbin_2t_1a}}
#'   \code{\link[rstan:sampling]{sampling}}
#'
#' @importFrom stats runif rbinom
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- stan_augbin_demo()
#' # I told you it was simple.
#' }
stan_augbin_demo <- function() {
  priors <- list(alpha_mean = 0, alpha_sd = 1,
                 beta_mean = 0, beta_sd = 1,
                 gamma_mean = 0, gamma_sd = 1,
                 sigma_mean = 0, sigma_sd = 1,
                 omega_lkj_eta = 1,
                 alpha_d1_mean = 0, alpha_d1_sd = 1,
                 gamma_d1_mean = 0, gamma_d1_sd = 1,
                 alpha_d2_mean = 0, alpha_d2_sd = 1,
                 gamma_d2_mean = 0, gamma_d2_sd = 1)
  # Scenario 1 of Table 1 in Wason & Seaman (2013)
  N <- 50
  sigma <- 1
  delta1 <- -0.356
  mu <- c(0.5 * delta1, delta1)
  Sigma = matrix(c(0.5 * sigma^2, 0.5 * sigma^2, 0.5 * sigma^2, sigma^2),
                 ncol = 2)
  alphaD <- -1.5
  gammaD <- 0
  set.seed(123456)
  y <- MASS::mvrnorm(n = N, mu, Sigma)
  z0 <- runif(N, min = 5, max = 10)
  z1 <- exp(y[, 1]) * z0
  z2 <- exp(y[, 2]) * z0
  d1 <- rbinom(N, size = 1, prob = gtools::inv.logit(alphaD + gammaD * z0))
  d2 <- rbinom(N, size = 1, prob = gtools::inv.logit(alphaD + gammaD * z1))
  tumour_size <- data.frame(z0, z1, z2) # Sizes in cm
  non_shrinkage_failure <- data.frame(d1, d2)
  # Fit
  ## Not run:
  fit <- stan_augbin(tumour_size, non_shrinkage_failure,
                     prior_params = priors, model = '2t-1a', seed = 123)
  fit
}
