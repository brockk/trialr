
#' Fit a CRM model to the incrementally observed outcomes on a trial pathway.
#'
#' Fit a continuous reassessment method (CRM) model to the outcomes cumulatively
#' observed at the end of each cohort in a trial pathway. E.g. if the trial
#' pathway is 1NN 2NN 3NT, we have three cohorts of two patients. This function
#' will fit the model to the following four states: before any patients have
#' been evaluated; after 1NN; after 1NN 2NN; and finally after 1NN 2NN 3NT.
#' This allows us to analyse how the trial model is evolving in its estimation
#' as trial data is accumulated.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{df_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given} and \code{tox}
#' parameters. See Details.
#' @param skeleton a vector of the prior guesses of toxicity at doses.
#' This should be a monotonically-increasing vector of numbers between 0 and 1.
#' @param target the target toxicity probability, a number between 0 and 1.
#' This value would normally be one of the values in \code{skeleton}, but that
#' is not a requirement.
#' @param model Character string to denote desired model. One of \code{empiric},
#' \code{logistic}, \code{logistic_gamma}, or \code{logistic2}.
#' The choice of model determines which extra parameters are required by
#' \code{...}. See Details.
#' @param verbose logical, TRUE to get log messages.
#' @param ... Extra parameters passed to \code{\link{stan_crm}}.
#'
#' @details
#' Different model choices require that different parameters are provided.
#' See below.
#'
#' @section Parameter requirements of \code{empiric} model:
#' \itemize{
#'   \item \code{beta_sd}
#' }
#'
#' @section Parameter requirements of \code{logistic} model:
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_mean}
#'   \item \code{beta_sd}
#' }
#'
#' @section Parameter requirements of \code{logistic_gamma} model:
#' \itemize{
#'   \item \code{a0}
#'   \item \code{beta_shape}
#'   \item \code{beta_inverse_scale}
#' }
#'
#' @section Parameter requirements of \code{logistic2} model:
#' \itemize{
#'   \item \code{alpha_mean}
#'   \item \code{alpha_sd}
#'   \item \code{beta_mean}
#'   \item \code{beta_sd}
#' }
#'
#' @return A \code{\link{list}} of \code{\link{dose_finding_path_node}} objects.
#'
#' @author Kristian Brock
#'
#' @seealso
#'   \code{\link{df_parse_outcomes}},
#'   \code{\link{stan_crm}},
#'   \code{\link{dose_finding_path_node}}
#'
#' @importFrom stringr str_trim
#' @export
#'
#' @examples
#' \dontrun{
#' # CRM example
#' target <- 0.25
#' skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
#' paths <- crm_path_analysis(
#'   outcome_str = '1NNN 2NTN 2NNN',
#'   skeleton = skeleton, target = target, model = 'empiric',
#'   beta_sd = 1, seed = 123, refresh = 0)
#' length(paths)  # 4
#' names(paths)[1]  # ""
#' names(paths)[2]  # "1NNN"
#' names(paths)[3]  # "1NNN 2NTN"
#' names(paths)[4]  # "1NNN 2NTN 2NNN"
#' # Each node is an analysis fit to the cumulative outcomes
#' # Converting to a tibble presents some nice tidyverse-related opportunities
#' library(tibble)
#' df <- as_tibble(paths)
#' df
#' }
crm_path_analysis <- function(outcome_str,
                              skeleton,
                              target,
                              model,
                              verbose = FALSE,
                              ...) {

  # Break outcomes pathway into cohorts.
  cohorts = parse_dose_finding_outcomes(outcome_str)

  num_doses <- length(skeleton)
  cache <- list()
  # Root node is the current scenario. Create and cache.
  root_node_id <- 1
  fit <- stan_crm(outcome_str = '', skeleton = skeleton, target = target,
                  model = model, ...)
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
    fit <- stan_crm(outcome_str = cumulative_outcome_str, skeleton = skeleton,
                    target = target, model = model, ...)
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
