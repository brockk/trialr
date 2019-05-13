
#' @title Parse a string of EffTox outcomes to binary vector notation.
#'
#' @description Parse a string of EffTox outcomes to the binary vector notation
#' required by Stan for model invocation. The outcome string describes the doses
#' given and outcomes observed. The format of the string is described in Brock
#' et al. (2017). The letters E, T, N and B are used to represents patients that
#' experienced (E)fficacy only, (T)oxicity only, (B)oth efficacy and toxicity,
#' and (N)either. These letters are concatenated after numerical dose-levels to
#' convey the outcomes of cohorts of patients. For instance, \code{2ETB}
#' represents a cohort of three patients that were treated at dose-level 2, and
#' experienced efficacy, toxicity and both events, respectively. The results of
#' cohorts are separated by spaces. Thus, \code{2ETB 1NN} extends our previous
#' example, where the next cohort of two were treated at dose-level 1 and both
#' patients experienced neither efficacy nor toxicity. See examples.
#'
#' We present the notation in the EffTox setting but it is applicable in
#' general seamless phase I/II dose-finding scenarios.
#'
#' @param outcome_string character string, conveying doses given and outcomes
#' observed.
#' @param as.list TRUE (be default) to return a \code{list};
#' FALSE to return a \code{data.frame}
#'
#' @return If \code{as.list == TRUE}, a list with elements \code{eff}, \code{tox},
#' \code{doses} and \code{num_patients}. These elements are congruent with those
#' of the same name in \code{efftox_params}.
#' If \code{as.list == FALSE}, a data.frame with columns \code{eff}, \code{tox},
#' and \code{doses}.
#'
#' @export
#'
#' @examples
#' x = efftox_parse_outcomes('1NNE 2EEN 3TBB')
#' x$num_patients == 9
#' x$eff == c(0, 0, 1, 1, 1, 0, 0, 1, 1)
#' sum(x$tox) == 3
#'
#' @references
#' Brock, K., Billingham, L., Copland, M., Siddique, S., Sirovica, M., & Yap, C.
#' (2017). Implementing the EffTox dose-finding design in the Matchpoint trial.
#' BMC Medical Research Methodology, 17(1), 112.
#' https://doi.org/10.1186/s12874-017-0381-x
#'
efftox_parse_outcomes <- function(outcome_string, as.list = TRUE) {

  if(outcome_string == '') {
    if(as.list) {
      return(list(
        doses = array(integer(length = 0)), eff = array(integer(length = 0)),
        tox = array(integer(length = 0)), num_patients = 0
      ))
    } else {
      return(data.frame(doses = array(integer(length = 0)),
                        eff = array(integer(length = 0)),
                        tox = array(integer(length = 0)))
      )
    }
  }

  # Matching is done by regex.

  # This pattern ensures that outcome_string is valid. It is the gate-keeper.
  # It allows leading and trailing white space and demands >0 cohort strings.
  # e.g. "2ENT 3TT 2E "
  valid_str_match <- '^\\s*(\\d+[ETNB]+\\s*)+$'
  # This pattern identifies the individual cohort strings, e.g. 2ENT
  cohort_str_match <- '\\d+[ETNB]+'
  # This pattern extracts the dose-level from a cohort string, e.g. 2
  dl_str_match <- '\\d+'
  # And this pattern extracts the outcomes from a cohort string, e.g ENT
  outcomes_match_str <- '[ETNB]+'

  if(stringr::str_detect(outcome_string, valid_str_match)) {
    doses <- eff <- tox <- c()
    cohort_strs <- stringr::str_extract_all(
      outcome_string, cohort_str_match)[[1]]
    for(cohort_str in cohort_strs) {
      c_dl <- as.integer(stringr::str_extract(cohort_str, dl_str_match))
      if(c_dl <= 0) stop('Dose-levels must be strictly positive integers.')
      c_outcomes <- stringr::str_extract(cohort_str, outcomes_match_str)

      these_doses <- rep(c_dl, nchar(c_outcomes))
      doses <- c(doses, these_doses)

      these_outcomes = stringr::str_split(c_outcomes, '')[[1]]
      these_eff = as.integer((these_outcomes == 'E') | (these_outcomes == 'B'))
      eff = c(eff, these_eff)
      these_tox = as.integer((these_outcomes == 'T') | (these_outcomes == 'B'))
      tox = c(tox, these_tox)
    }
  } else {
    stop(paste0('"', outcome_string, '" is not a valid outcome string.
                A valid example is "1N 2EE 3TB 2BE"'))
  }

  if(as.list) {
    return(list(
      doses = doses, eff = eff, tox = tox, num_patients = length(doses)
    ))
  } else {
    return(data.frame(doses = doses, eff = eff, tox = tox))
  }
}

