#' @title Parse a string of phase I/II dose-finding trial outcomes.
#'
#' @description Parse a string of phase I/II dose-finding trial outcomes.
#' Phase I/II trials conduct dose-finding by efficacy and toxicity outcomes.
#'
#' @param outcome_string character representing doses given, outcomes
#' observed, and timing of analyses. See Description.
#'
#' @description Parse a string of phase I/II dose-finding outcomes to a list.
#' The outcome string describes the doses given, efficacy and toxicity outcomes
#' observed and the timing of analyses that recommend a dose. The format of the
#' string is described in Brock _et al_. (2017). The letters E, T, N & B are
#' used to represents patients that experienced (E)fficacy, (T)oxicity,
#' (N)either and (B)oth. These letters are concatenated after numerical
#' dose-levels to convey the outcomes of cohorts of patients.
#' For instance, \code{2NET} represents a cohort of three patients that were
#' treated at dose-level 2, one of whom experienced toxicity only, one that
#' experienced efficacy only, and one that had neither.
#' The results of cohorts are separated by spaces and it is assumed that a
#' dose-finding decision takes place at the end of a cohort. Thus,
#' \code{2NET 1NN} builds on our previous example, where the next cohort of two
#' were treated at dose-level 1 and neither of these patients experienced
#' either event. See examples.
#'
#' @return a list with a slot for each cohort. Each cohort slot is itself a
#' list, containing elements:
#' * \code{dose}, the integer dose delivered to the cohort;
#' * \code{outcomes}, a character string representing the \code{E}, \code{T}
#' \code{N} or \code{B} outcomes for the patients in this cohort.
#'
#' @export
#'
#' @examples
#' x = parse_eff_tox_dose_finding_outcomes('1NEN 2ENT 3TB')
#' length(x)
#' x[[1]]$dose
#' x[[1]]$outcomes
#' x[[2]]$dose
#' x[[2]]$outcomes
#' x[[3]]$dose
#' x[[3]]$outcomes
#'
#' @references
#' Brock, K., Billingham, L., Copland, M., Siddique, S., Sirovica, M., & Yap, C.
#' (2017). Implementing the EffTox dose-finding design in the Matchpoint trial.
#' BMC Medical Research Methodology, 17(1), 112.
#' https://doi.org/10.1186/s12874-017-0381-x
parse_eff_tox_dose_finding_outcomes <- function(outcome_string) {

  if(outcome_string == '') return(list())

  # Matching is done by regex.
  # This pattern ensures that outcome_string is valid. It is the gate-keeper.
  # It allows leading and trailing white space and demands >0 cohort strings.
  # e.g. "2NNT 3TT 2N "
  valid_str_match <- '^\\s*(\\d+[ETNB]+\\s*)+$'
  # This pattern identifies the individual cohort strings, e.g. "2NET"
  cohort_str_match <- '\\d+[ETNB]+'
  # This pattern extracts the dose-level from a cohort string, e.g. "2"
  dl_str_match <- '\\d+'
  # And this pattern extracts the outcomes from a cohort string, e.g "NET"
  outcomes_match_str <- '[ETNB]+'

  cohorts <- list()
  cohort_id <- 1

  if(stringr::str_detect(outcome_string, valid_str_match)) {
    cohort_strs <- stringr::str_extract_all(
      outcome_string, cohort_str_match)[[1]]
    for(cohort_str in cohort_strs) {
      c_dl <- as.integer(stringr::str_extract(cohort_str, dl_str_match))
      if(c_dl <= 0) stop('Dose-levels must be strictly positive integers.')
      c_outcomes <- stringr::str_extract(cohort_str, outcomes_match_str)
      cohorts[[cohort_id]] <- list(dose = c_dl, outcomes = c_outcomes)
      cohort_id <- cohort_id + 1
    }
  } else {
    stop(paste0('"', outcome_string, '" is not a valid outcome string.
                A valid example is "1N 2NE 3BB 2NT"'))
  }

  cohorts
}
