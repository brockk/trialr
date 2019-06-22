
#' Calculate the probability of success.
#'
#' @param x an R object of class \code{"augbin_fit"}
#' @param ... arguments passed to other methods
#' @export
#' @rdname prob_success
prob_success <- function(x, ...) {
  UseMethod("prob_success", x)
}

#' Calculate the probability of success for an augbin_2t_1a_fit object.
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
#' Defaults to log(0.7).
#' @param probs pair of probabilities to use to calculate the credible interval
#' for the probability of success.
#' @param newdata data for which to infer the probability of success.
#' A dataframe-like object with baseline tumour sizes in first column, and first
#' and second post-baseline tumour sizes in columns 2 and 3. Omitted by default.
#' When omitted, newdata is set to be the \code{fit$tumour_size}.
#'
#' @return Object of class \code{\link[tibble]{tibble}}
#'
#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom gtools inv.logit
#' @importFrom purrr map2
#' @importFrom stats quantile pnorm
#' @export
#' @rdname prob_success
#'
#' @examples
#' \dontrun{
#' fit <- stan_augbin_demo()
#' prob_success(fit, y2_upper = log(0.7))
#' }
prob_success.augbin_2t_1a_fit <- function(x,
                                          y1_lower = -Inf, y1_upper = Inf,
                                          y2_lower = -Inf, y2_upper = log(0.7),
                                          probs = c(0.025, 0.975),
                                          newdata = NULL,
                                          ...) {

  if(length(probs) != 2)
    stop('probs should be a pair of probabilities between 0 and 1.')
  if(any(probs > 1) | any(probs < 0))
    stop('probs should be a pair of probabilities between 0 and 1.')
  if(is.null(newdata))
    newdata <- x$tumour_size
  num_patients <- nrow(newdata)

  alpha <- rstan::extract(x$fit)$alpha
  beta <- rstan::extract(x$fit)$beta
  gamma <- rstan::extract(x$fit)$gamma
  sigma1 <- rstan::extract(x$fit)$sigma[,1]
  sigma2 <- rstan::extract(x$fit)$sigma[,2]
  alphaD1 <- rstan::extract(x$fit)$alphaD1
  gammaD1 <- rstan::extract(x$fit)$gammaD1
  alphaD2 <- rstan::extract(x$fit)$alphaD2
  gammaD2 <- rstan::extract(x$fit)$gammaD2

  .get_prob_success <- function(z0, z1) {
    mu1 <- alpha + gamma * z0
    mu2 <- beta + gamma * z0
    prob_success_samp_id <- (
      (1 - inv.logit(alphaD1 + z0 * gammaD1)) *
        (1 - inv.logit(alphaD2 + z1 * gammaD2)) *
        ( pnorm(q = y1_upper, mean = mu1, sd = sigma1) -
            pnorm(q = y1_lower, mean = mu1, sd = sigma1) ) *
        ( pnorm(q = y2_upper, mean = mu2, sd = sigma2) -
            pnorm(q = y2_lower, mean = mu2, sd = sigma2) )
    )
    prob_success_samp_id
  }

  tibble(id = 1:num_patients,
         z0 = newdata[, 1],
         z1 = newdata[, 2]) %>%
    mutate(
      prob_success_samp = map2(.data$z0, .data$z1, .get_prob_success),
      prob_success = map_dbl(.data$prob_success_samp, mean),
      lower = map_dbl(.data$prob_success_samp, quantile, probs = min(probs)),
      upper = map_dbl(.data$prob_success_samp, quantile, probs = max(probs)),
      ci_width = .data$upper - .data$lower
    ) -> inferences

  inferences
}
