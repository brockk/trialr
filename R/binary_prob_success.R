

#' Calculate the binary probability of success.
#'
#' @param x an R object of class \code{"augbin_fit"}
#' @param ... arguments passed to other methods
#' @export
#' @rdname binary_prob_success
binary_prob_success <- function(x, ...) {
  UseMethod("binary_prob_success", x)
}


#' Calculate the binary probability of success from an augbin_2t_1a_fit object.
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
#' @param conf.level confidence level for interval.
#'
#' @return data.frame
#'
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate summarise
#' @importFrom binom binom.confint
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- stan_augbin_demo()
#' binary_prob_success(fit, y2_upper = log(0.7))
#' }
binary_prob_success.augbin_2t_1a_fit <- function(fit,
                                                 y1_lower = -Inf,
                                                 y1_upper = Inf,
                                                 y2_lower = -Inf,
                                                 y2_upper = log(0.7),
                                                 conf.level = 0.95,
                                                 ...) {

  as_tibble(fit) %>%
    mutate(success = (.data$d1 == 0) & (.data$d2 == 0) &
             (.data$y1 < y1_upper) & (.data$y1 > y1_lower) &
             (.data$y2 < y2_upper) & (.data$y2 > y2_lower)) %>%
    summarise(num_success = sum(.data$success)) %>% .[[1]] -> num_success
  binom.confint(x = num_success, n = fit$num_patients,
                conf.level = conf.level, ...) %>%
    mutate(ci_width = .data$upper - .data$lower)
}
