
#' Fit an EffTox model to the incrementally observed outcomes on a trial pathway.
#'
#' Fit a EffTox model to the outcomes cumulatively observed at the end of each
#' cohort in a trial pathway. E.g. if the trial pathway is 1EN 2NN 3BT, we have
#' three cohorts of two patients. This function will fit the model to the
#' following four states: before any patients have been evaluated; after 1EN;
#' after 1EN 2NN; and finally after 1EN 2NN 3BT. This allows us to analyse how
#' the trial model is evolving in its estimation as trial data is accumulated.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{efftox_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given} and \code{tox}
#' parameters. See Details.
#' @param verbose logical, TRUE to get log messages.
#' @param ... All other parameters are passed to \code{\link{stan_efftox}}.
#'
#'
#' @return A \code{\link{list}} of \code{\link{dose_finding_path_node}} objects.
#'
#' @author Kristian Brock
#'
#' @seealso
#'   \code{\link{efftox_parse_outcomes}},
#'   \code{\link{stan_efftox}},
#'   \code{\link{dose_finding_path_node}}
#'
#' @importFrom stringr str_trim
#' @export
#'
#' @examples
#' \dontrun{
#' # EffTox example
#' paths <- efftox_path_analysis(
#'   outcome_str = '1NNN 2NEN 3NEB',
#'   real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'   efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
#'   p_e = 0.1, p_t = 0.1,
#'   eff0 = 0.5, tox1 = 0.65,
#'   eff_star = 0.7, tox_star = 0.25,
#'   alpha_mean = -7.9593, alpha_sd = 3.5487,
#'   beta_mean = 1.5482, beta_sd = 3.5018,
#'   gamma_mean = 0.7367, gamma_sd = 2.5423,
#'   zeta_mean = 3.4181, zeta_sd = 2.4406,
#'   eta_mean = 0, eta_sd = 0.2,
#'   psi_mean = 0, psi_sd = 1, seed = 123, refresh = 0)
#'
#' length(paths)  # 4
#' names(paths)[1]  # ""
#' names(paths)[2]  # "1NNN"
#' names(paths)[3]  # "1NNN 2NEN"
#' names(paths)[4]  # "1NNN 2NEN 3NEB"
#' # Each node is an analysis fit to the cumulative outcomes
#' # Converting to a tibble presents some nice tidyverse-related opportunities
#' library(tibble)
#' df <- as_tibble(paths)
#' df
#' }
efftox_path_analysis <- function(outcome_str,
                                 verbose = FALSE,
                                 ...) {

  # Break outcomes pathway into cohorts.
  cohorts = parse_eff_tox_dose_finding_outcomes(outcome_str)

  cache <- list()
  # Root node is the current scenario. Create and cache.
  root_node_id <- 1
  fit <- stan_efftox(outcome_str = '', ...)
  next_dose <- fit$recommended_dose
  root <- dose_finding_path_node(node_id = root_node_id,
                                 parent_node_id = NA,
                                 depth = 0,
                                 outcomes = '',
                                 next_dose = next_dose,
                                 fit = fit,
                                 parent_fit = NULL)
  cache[['']] <- root
  parent <- root

  # Fit model to each successive cohort, accumulating outcomes as we go.
  cumulative_outcome_str = ''
  node_id <- root_node_id + 1
  for(cohort in cohorts) {
    cumulative_outcome_str <- str_trim(paste0(cumulative_outcome_str, ' ',
                                              cohort$dose, cohort$outcomes))
    if(verbose) print(paste0('Running ', cumulative_outcome_str))
    fit <- stan_efftox(outcome_str = cumulative_outcome_str, ...)
    # Cache
    node <- dose_finding_path_node(node_id = node_id,
                                   parent_node_id = parent$.node,
                                   depth = node_id - 1,
                                   outcomes = cumulative_outcome_str,
                                   next_dose = fit$recommended_dose,
                                   fit = fit,
                                   parent_fit = parent$fit)
    cache[[cumulative_outcome_str]] <- node
    parent <- node
    node_id <- node_id + 1
  }

  # Add type and return
  class(cache) <- c("dose_finding_paths", "list")
  cache
}
