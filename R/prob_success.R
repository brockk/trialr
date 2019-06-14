
#' Calculate the probability of success.
#'
#' @param x an R object of class \code{"augbin_fit"}
#' @param ... arguments passed to other methods
#' @export
#' @rdname prob_success
#' @examples
#' \dontrun{
#' # TODO
#' }
prob_success <- function(x, ...) {
  UseMethod("prob_success", x)
}


#' Calculate the probability of success.
#'
#' @param y1_lower numeric, minimum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 1 to baseline.
#' Defaults to negative infinity.
#' @param y1_upper numeric, maximum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 1 to baseline.
#' Defaults to positive infinity.
#' @param y2_lower numeric, minimum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 2 to baseline.
#' @param y2_upper numeric, maximum threshold to constitute success,
#' scrutinising the log of the tumour size ratio comparing time 2 to baseline.
#' Defaults to log(0.7), i.e. tumour shrinkage of at least 30% is success.
#' @param probs pair of probabilities to use to calculate the credible interval
#' for the probability of success.
#'
#' @return numerical vector of probabilities.
#'
#' @rdname prob_success
#'
#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom gtools inv.logit
#' @export
prob_success.augbin_2t_1a_fit <- function(fit,
                                          y1_lower = -Inf, y1_upper = Inf,
                                          y2_lower = -Inf, y2_upper = log(0.7),
                                          probs = c(0.025, 0.975),
                                          ...) {

  if(length(probs) != 2)
    stop('probs should be a pair of probabilities between 0 and 1.')
  if(any(probs > 1) | any(probs < 0))
    stop('probs should be a pair of probabilities between 0 and 1.')

  alpha <- rstan::extract(fit$fit)$alpha
  beta <- rstan::extract(fit$fit)$beta
  gamma <- rstan::extract(fit$fit)$gamma
  sigma1 <- rstan::extract(fit$fit)$sigma[,1]
  sigma2 <- rstan::extract(fit$fit)$sigma[,2]
  alphaD1 <- rstan::extract(fit$fit)$alphaD1
  gammaD1 <- rstan::extract(fit$fit)$gammaD1
  alphaD2 <- rstan::extract(fit$fit)$alphaD2
  gammaD2 <- rstan::extract(fit$fit)$gammaD2

  .get_prob_success <- function(z0, z1) {
    mu1 <- alpha + gamma * z0
    mu2 <- beta + gamma * z0
    prob_success_samp_tno <- (
      (1 - inv.logit(alphaD1 + z0 * gammaD1)) *
        (1 - inv.logit(alphaD2 + z1 * gammaD2)) *
        ( pnorm(q = y1_upper, mean = mu1, sd = sigma1) -
            pnorm(q = y1_lower, mean = mu1, sd = sigma1) ) *
        ( pnorm(q = y2_upper, mean = mu2, sd = sigma2) -
            pnorm(q = y2_lower, mean = mu2, sd = sigma2) )
    )
    prob_success_samp_tno
  }

  tibble(tno = 1:fit$num_patients,
         z0 = fit$tumour_size[, 1],
         z1 = fit$tumour_size[, 2]) %>%
    mutate(
      prob_success_samp = map2(.data$z0, .data$z1, .get_prob_success),
      prob_success = map_dbl(.data$prob_success_samp, mean),
      lower = map_dbl(.data$prob_success_samp, quantile, probs = min(probs)),
      upper = map_dbl(.data$prob_success_samp, quantile, probs = max(probs)),
      ci_width = .data$upper - .data$lower
    ) -> inferences
  inferences
}
