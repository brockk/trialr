
#' @title Calculate dose-transition pathways for an EffTox study
#'
#' @description Calculate dose-transition pathways for an EffTox study.
#' The function \code{\link{efftox_dtps_to_dataframe}} performs a similar
#' function, but is much less-flexible.
#'
#' @param cohort_sizes vector of future cohort sizes, i.e. positive integers.
#' E.g. To calculate paths for the the next cohort of two followed by another
#' cohort of three, use \code{cohort_sizes = c(2, 3)}.
#' @param previous_outcomes Outcomes observed hitherto in the syntax required
#' by \code{\link{efftox_parse_outcomes}}.
#' @param next_dose the dose-level to be given to the immediately next cohort.
#' @param user_dose_func optional delegate for deciding dose. A function that
#' takes a \code{\link{efftox_fit}} as the sole argument and returns the integer
#' (1-based) dose-level to be given next, or NA to show that no dose should be
#' chosen and the trial stopped. This function gives the user the opportunity to
#' build in custom behaviour to tailor the dose selection decision in response
#' to the insights garnered by the fit model, or recommend that a trial path
#' be halted immediately. If omitted, the dose ordinarily chosen by the model is
#' used. An example is given below.
#' @param verbose logical, TRUE to get progress messages.
#' @param i_am_patient logical, TRUE to show your tolerance for waiting for over 100
#' models to fit. Set to FALSE by default.
#' @param ... extra params passed to \code{rstan::sampling}.
#'
#' @return dose pathways in a \code{data.frame}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate paths for the first cohort of 3 in Thall et al 2014 example
#' paths1 <- efftox_dtps(cohort_sizes = c(3), next_dose = 1,
#'                       real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'                       efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
#'                       p_e = 0.1, p_t = 0.1,
#'                       eff0 = 0.5, tox1 = 0.65,
#'                       eff_star = 0.7, tox_star = 0.25,
#'                       alpha_mean = -7.9593, alpha_sd = 3.5487,
#'                       beta_mean = 1.5482, beta_sd = 3.5018,
#'                       gamma_mean = 0.7367, gamma_sd = 2.5423,
#'                       zeta_mean = 3.4181, zeta_sd = 2.4406,
#'                       eta_mean = 0, eta_sd = 0.2,
#'                       psi_mean = 0, psi_sd = 1, seed = 123)
#'
#'
#'
#' # Calculate paths for the next two cohorts of 2, in an in-progress trial
#' # Warning: this create 100 paths. It will run for a minute or two.
#' paths2 <- efftox_dtps(cohort_sizes = c(2, 2),
#'                       previous_outcomes = '1NN 2EE',
#'                       next_dose = 1,
#'                       real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'                       efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
#'                       p_e = 0.1, p_t = 0.1,
#'                       eff0 = 0.5, tox1 = 0.65,
#'                       eff_star = 0.7, tox_star = 0.25,
#'                       alpha_mean = -7.9593, alpha_sd = 3.5487,
#'                       beta_mean = 1.5482, beta_sd = 3.5018,
#'                       gamma_mean = 0.7367, gamma_sd = 2.5423,
#'                       zeta_mean = 3.4181, zeta_sd = 2.4406,
#'                       eta_mean = 0, eta_sd = 0.2,
#'                       psi_mean = 0, psi_sd = 1, seed = 123,
#'                       i_am_patient = TRUE)
#'
#' # Paths can be converted to a tibble
#' library(tibble)
#' library(dplyr)
#' df <- as_tibble(paths2)
#' df %>% print(n = 200)
#'
#' # And shaped in a wide format
#' spread_paths(df %>% select(-fit, -parent_fit, -dose_index)) %>%
#'   print(n = 100)
#' # Incredibly, there are 100 ways these two cohorts of two can end up.
#'
#'
#'
#' # An example with a custom dose selection function.
#' # Define a function to select the maximal utility dose, no matter what.
#' # Note: this diverges from the original authors' intentions; we provide this
#' # for illustration only!
#' max_utility_dose <- function(efftox_fit) {
#'   return(which.max(efftox_fit$utility))
#' }
#' # Fit the paths, providing the user_dose_func parameter
#' # Warning: this create 100 paths. It will run for a minute or two.
#' paths3 <- efftox_dtps(cohort_sizes = c(2, 2),
#'                       previous_outcomes = '1NN 2EE',
#'                       next_dose = 1,
#'                       real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
#'                       efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
#'                       p_e = 0.1, p_t = 0.1,
#'                       eff0 = 0.5, tox1 = 0.65,
#'                       eff_star = 0.7, tox_star = 0.25,
#'                       alpha_mean = -7.9593, alpha_sd = 3.5487,
#'                       beta_mean = 1.5482, beta_sd = 3.5018,
#'                       gamma_mean = 0.7367, gamma_sd = 2.5423,
#'                       zeta_mean = 3.4181, zeta_sd = 2.4406,
#'                       eta_mean = 0, eta_sd = 0.2,
#'                       psi_mean = 0, psi_sd = 1,
#'                       user_dose_func = max_utility_dose,
#'                       seed = 123, i_am_patient = TRUE)
#'
#' # We can see where the dose-selections differ at the second future cohort
#' # by joining these paths to those calculated in the previous example:
#' left_join(
#'   as_tibble(paths2)%>%
#'     select(.node, .parent, .depth, outcomes, model_dose = next_dose),
#'   as_tibble(paths3) %>%
#'     select(.node, user_dose = next_dose),
#'   by = '.node'
#' ) %>% spread_paths() %>%
#'   filter(model_dose2 != user_dose2)
#' # They differ in many places. The user defined functions sometimes selects
#' # higher doses; sometimes lower.
#' }
#'
#' @seealso
#'   \code{\link{efftox_parse_outcomes}},
#'   \code{\link{stan_efftox}},
#'   \code{\link{efftox_path_analysis}},
#'   \code{\link{dose_finding_path_node}}
#'
#' @references
#' Yap C, Billingham LJ, Cheung YK, Craddock C, O’Quigley J.
#' Dose transition pathways: The missing link between complex dose-finding
#' designs and simple decision-making. Clinical Cancer Research.
#' 2017;23(24):7440-7447. doi:10.1158/1078-0432.CCR-17-0582
#'
#' Brock K, Billingham L, Copland M, Siddique S, Sirovica M, Yap C.
#' Implementing the EffTox dose-finding design in the Matchpoint trial.
#' BMC Medical Research Methodology. 2017;17(1):112.
#' doi:10.1186/s12874-017-0381-x
#'
efftox_dtps <- function(cohort_sizes,
                        previous_outcomes = '',
                        next_dose = NULL,
                        user_dose_func = NULL,
                        verbose = FALSE,
                        i_am_patient = FALSE,
                        ...) {

  if(!all(cohort_sizes == ceiling(cohort_sizes)))
    stop('cohort_sizes must be stricly positive integers.')
  if(!all(cohort_sizes > 0))
    stop('cohort_sizes must be stricly positive integers.')

  if(nchar(previous_outcomes) > 0)
    dat <- efftox_parse_outcomes(previous_outcomes)
  else
    dat <- list(doses = c(), eff = c(), tox = c(), num_patients = 0)
  previous_doses <- dat$doses
  previous_eff <- dat$eff
  previous_tox <- dat$tox
  previous_num_patients <- dat$num_patients

  # Calculate feasible outcome combinations by cohort
  outcomes <- c('E', 'T', 'N', 'B')
  cohort_paths <- lapply(cohort_sizes,
                         function(x) gtools::combinations(n = 4, r = x,
                                                          v = outcomes,
                                                          repeats.allowed=TRUE))
  # Flatten cohort outcomes
  cohort_paths <- lapply(cohort_paths, function(x) apply(x, 1, paste0,
                                                         collapse = ''))

  # Calculate pathways
  cohort_paths <- expand.grid(cohort_paths, stringsAsFactors = FALSE)

  num_paths <- nrow(cohort_paths)
  if(num_paths >= 50 & num_paths < 100) {
    message(paste0('You have requested ', num_paths,
                   ' model evaluations. Be patient.'))
  }
  if(num_paths >= 100 & !i_am_patient) {
    stop(paste0('You have requested ', num_paths,
                ' model evaluations but also flagged your impatience.',
                ' Run again with i_am_patient = TRUE'))
  }

  # Cache DTP calculations to avoid needless repetition
  cache <- list()
  # Root node is the current scenario
  root_node_id <- 1
  fit <- stan_efftox(outcome_str = previous_outcomes, ...)
  if(is.null(next_dose)) {
    if(is.null(user_dose_func))
      next_dose <- fit$recommended_dose
    else
      next_dose <- user_dose_func(fit)
  }
  root <- dose_finding_path_node(node_id = root_node_id,
                                 parent_node_id = NA,
                                 depth = 0,
                                 outcomes = '',
                                 next_dose = next_dose,
                                 fit = fit,
                                 parent_fit = NULL)
  cache[['']] <- root
  node_id <- root_node_id + 1

  for(i in 1:nrow(cohort_paths)) {
    cohort_path <- cohort_paths[i, ]
    cohort_dose <- next_dose
    dtp <- ""
    parent <- root

    for(j in 1:length(cohort_path)) {
      # If the dose is NA, this trial has stopped.
      if(!is.na(cohort_dose)) {
        dtp <- ifelse(nchar(dtp) > 0,
                      paste0(dtp, ' ', cohort_dose, cohort_path[j]),
                      paste0(cohort_dose, cohort_path[j])
        )
        if(dtp %in% names(cache)) {
          # Fetch from cache
          if(verbose) print(paste0('Fetching ', dtp, ' from cache'))
          parent <- cache[[dtp]]
          cohort_dose <- parent$next_dose
        } else {
          # Fit model for path, and cache.
          these_outcomes <- efftox_parse_outcomes(dtp)
          dat$doses <- array(c(previous_doses, these_outcomes$doses))
          dat$eff <- array(c(previous_eff, these_outcomes$eff))
          dat$tox <- array(c(previous_tox, these_outcomes$tox))
          dat$num_patients <- previous_num_patients +
            these_outcomes$num_patients
          if(verbose) print(paste0('Running ', dtp))
          fit <- stan_efftox(doses_given = dat$doses, eff = dat$eff,
                             tox = dat$tox, ...)
          if(is.null(user_dose_func))
            cohort_dose <- fit$recommended_dose
          else
            cohort_dose <- user_dose_func(fit)

          # Cache
          node <- dose_finding_path_node(node_id = node_id,
                                         parent_node_id = parent$.node,
                                         depth = j,
                                         outcomes = as.character(cohort_path[j]),
                                         next_dose = cohort_dose,
                                         fit = fit,
                                         parent_fit = parent$fit)
          cache[[dtp]] <- node
          parent <- node
          node_id <- node_id + 1
        }
      }
    }
  }

  class(cache) <- c("dose_finding_paths", "list")
  cache
}
