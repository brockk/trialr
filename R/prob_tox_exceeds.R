
#' Calculate the probability that the rate of toxicity exceeds some threshold
#'
#' @param x an R object of class \code{"dose_finding_fit"}
#' @param ... arguments passed to other methods
#' @return numerical vector of probabilities
#' @export
#' @rdname prob_tox_exceeds
#' @examples
#' \dontrun{
#' # CRM example
#' target <- 0.2
#' fit <- stan_crm('1N 2N 3T', skeleton = c(0.1, 0.2, 0.35, 0.6),
#'                  target = target, model = 'empiric', beta_sd = sqrt(1.34),
#'                  seed = 123)
#' prob_tox_exceeds(fit, target)
#' }
prob_tox_exceeds <- function(x, ...) {
  UseMethod("prob_tox_exceeds", x)
}

#' Calculate the probability that the rate of toxicity exceeds some threshold
#'
#' @param threshold numeric, threshold value.
#' @return numerical vector of probabilities
#' @rdname prob_tox_exceeds
#' @importFrom magrittr "%>%"
#' @importFrom tidybayes gather_draws
#' @importFrom dplyr mutate summarise ungroup select
#' @export
prob_tox_exceeds.dose_finding_fit <- function(x, threshold, ...) {
  prob_tox <- dose <- .value <- TooToxic <- ProbToxExceeds <- . <- NULL
  x %>%
    gather_draws(prob_tox[dose]) %>%
    mutate(TooToxic = .value > threshold) %>%
    summarise(ProbToxExceeds = mean(TooToxic)) %>%
    ungroup() %>%
    select(ProbToxExceeds) %>% .[[1]]
}
