
#' Calculate dose-transition pathways for a CRM study
#'
#' Calculate dose-transition pathways (DTPs, Yap et al, 2017) for a dose-finding
#' trial using the continual reassessment method (CRM) design. DTPs are a
#' glimpse into the future for an in-progress trial. They tell us what the model
#' would advise for all feasible future outcomes. They can be used in the design
#' stages to detect possible undesirable behaviour. They can be used during the
#' trial to aid planning and understanding.
#'
#' @param skeleton a vector of the prior guesses of toxicity at doses.
#' This should be a monotonically-increasing vector of numbers between 0 and 1.
#' @param target the target toxicity probability, a number between 0 and 1.
#' This value would normally be one of the values in \code{skeleton}, but that
#' is not a requirement.
#' @param model Character string to denote desired model. One of \code{empiric},
#' \code{logistic}, \code{logistic_gamma}, or \code{logistic2}.
#' The choice of model determines which extra parameters are required by
#' \code{...}. See Details.
#' @param cohort_sizes vector of future cohort sizes, i.e. positive integers.
#' E.g. To calculate paths for the the next cohort of two followed by another
#' cohort of three, use \code{cohort_sizes = c(2, 3)}.
#' @param previous_outcomes Outcomes observed hitherto in the syntax required
#' by \code{\link{df_parse_outcomes}}.
#' @param next_dose optional, integer (1-based) dose-level to be given to the
#' next cohort. If omitted, the dose suggested by the model is used.
#' @param user_dose_func optional delegate for deciding dose. A function that
#' takes a \code{\link{crm_fit}} as the sole argument and returns the integer
#' (1-based) dose-level to be given next, or NA to show that no dose should be
#' chosen and the trial stopped. This function gives the user the opportunity to
#' build in custom behaviour to tailor the dose selection decision in response
#' to the insights garnered by the fit model, or recommend that a trial path
#' be halted immediately. If omitted, the dose ordinarily chosen by the model is
#' used. An example is given below.
#' @param verbose logical, TRUE to get log messages.
#' @param i_am_patient logical. The number of paths to analyse grows faster than
#' linearly in the number of future cohorts to resolve. Fitting many models by
#' MCMC can take a long time. This function will not proceed unless you signify
#' your patience when the number of paths to reolve exceeds 100.
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
#'   \code{\link{crm_path_analysis}},
#'   \code{\link{dose_finding_path_node}}
#'
#' @export
#'
#' @references
#' Yap C, Billingham LJ, Cheung YK, Craddock C, O’Quigley J.
#' Dose transition pathways: The missing link between complex dose-finding
#' designs and simple decision-making. Clinical Cancer Research.
#' 2017;23(24):7440-7447. doi:10.1158/1078-0432.CCR-17-0582
#'
#' @examples
#' \dontrun{
#' target <- 0.25
#' skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
#'
#' # Run DTPs for the first two cohorts of two for new a trial:
#' paths <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric',
#'                   cohort_sizes = c(2, 2), next_dose = 3, beta_sd = 1)
#' length(paths)  # 13
#'
#' library(tibble)
#' df <- as_tibble(paths)
#' df
#'
#'
#' # Run DTPs for the next cohort of three in a trial that has already treated
#' # six patients, seeing some toxicity at dose-level 3:
#' paths2 <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric',
#'                    cohort_sizes = c(3), previous_outcomes = '2NNN 3TTN',
#'                    beta_sd = 1)
#' length(paths2)  # 5
#' as_tibble(paths2)
#' # We see that de-escalation to dose-level 2 should occur now, and that any
#' # further toxicity will result in advice for further de-escalation to
#' # dose-level 1.
#'
#'
#' # An example with a custom dose selection function
#' paths3 <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric',
#'                    cohort_sizes = c(3, 3), previous_outcomes = '2NN 3TN',
#'                    next_dose = 2, beta_sd = 1,
#'                    user_dose_func = function(x) {
#'                      careful_escalation(x, tox_threshold = target + 0.1,
#'                                         certainty_threshold = 0.7)
#'                    }, seed = 123, refresh = 0)
#' spread_paths(as_tibble(paths3) %>% select(-fit, -parent_fit, -dose_index))
#' # Stopping is recommended when the dose selection function returns NA.
#' }
crm_dtps <- function(skeleton,
                     target,
                     model,
                     cohort_sizes,
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

  max_depth <- length(cohort_sizes)
  num_paths = 1 + sum(sapply(1:max_depth,
                             function(i) prod((cohort_sizes + 1)[1:i])))
  if(num_paths >= 50 & num_paths < 100) {
    message(paste0('You have requested ', num_paths,
                   ' model evaluations. Be patient.'))
  }
  if(num_paths >= 100 & !i_am_patient) {
    stop(paste0('You have requested ', num_paths,
                ' model evaluations but also flagged your impatience.',
                ' Run again with i_am_patient = TRUE'))
  }

  if(nchar(previous_outcomes) > 0)
    dat <- df_parse_outcomes(previous_outcomes)
  else
    dat <- list(doses = c(), tox = c(), num_patients = 0)
  num_doses <- length(skeleton)
  previous_doses <- dat$doses
  previous_tox <- dat$tox
  previous_num_patients <- dat$num_patients
  outcomes <- c('T', 'N')

  # Calculate feasible outcome combinations by cohort
  cohort_paths <- lapply(cohort_sizes,
                         function(x) gtools::combinations(n = 2, r = x,
                                                          v = outcomes,
                                                          repeats.allowed=TRUE))
  # Flatten cohort outcomes
  cohort_paths <- lapply(cohort_paths, function(x) apply(x, 1, paste0,
                                                         collapse = ''))

  # Calculate pathways
  cohort_paths <- expand.grid(cohort_paths, stringsAsFactors = FALSE)

  # Cache DTP calculations to avoid needless repetition
  cache <- list()
  # Root node is the current scenario
  root_node_id <- 1
  fit <- stan_crm(outcome_str = previous_outcomes, skeleton = skeleton,
                  target = target, model = model, ...)
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
          these_outcomes <- df_parse_outcomes(dtp)
          dat$doses <- array(c(previous_doses, these_outcomes$doses))
          dat$tox <- array(c(previous_tox, these_outcomes$tox))
          dat$num_patients <- previous_num_patients +
            these_outcomes$num_patients
          if(verbose) print(paste0('Running ', dtp))
          fit <- stan_crm(skeleton = skeleton, target = target, model = model,
                          doses_given = dat$doses, tox = dat$tox, ...)
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
