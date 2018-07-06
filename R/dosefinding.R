
#' @title Parse a string of dose-finding trial outcomes to binary vector notation.
#'
#' @description Parse a string of dose-finding trial outcomes to the binary
#' vector notation required by Stan for model invocation. The outcome string
#' describes the doses given and outcomes observed. The format of the string is
#' the pure phase I analogue to that described in Brock et al. (2017).
#' The letters T and N are used to represents patients that experienced
#' (T)oxicity and (N)o toxicity. These letters are concatenated after numerical
#' dose-levels to convey the outcomes of cohorts of patients.
#' For instance, \code{2NNT} represents a cohort of three patients that were
#' treated at dose-level 2, one of whom experienced toxicity, and two that did
#' not. The results of cohorts are separated by spaces. Thus, \code{2NNT 1NN}
#' extends our previous example, where the next cohort of two were treated at
#' dose-level 1 and neither experienced toxicity. See examples.
#'
#' @param outcome_string character string, conveying doses given and outcomes
#' observed.
#' @param as.list TRUE (be default) to return a \code{list};
#' FALSE to return a \code{data.frame}
#'
#' @return If \code{as.list == TRUE}, a list with elements \code{tox},
#' \code{doses} and \code{num_patients}. These elements are congruent with those
#' of the same name in \code{crm_params}, for example.
#' If \code{as.list == FALSE}, a data.frame with columns \code{tox} and
#' \code{doses}.
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
df_parse_outcomes <- function(outcome_string, as.list = TRUE) {

  # Matching is done by regex.

  # This pattern ensures that outcome_string is valid. It is the gate-keeper.
  # It allows leading and trailing white space and demands >0 cohort strings.
  # e.g. "2NNT 3TT 2N "
  valid_str_match <- '^\\s*(\\d+[NT]+\\s*)+$'
  # This pattern identifies the individual cohort strings, e.g. 2NNT
  cohort_str_match <- '\\d+[NT]+'
  # This pattern extracts the dose-level from a cohort string, e.g. 2
  dl_str_match <- '\\d+'
  # And this pattern extracts the outcomes from a cohort string, e.g NNT
  outcomes_match_str <- '[NT]+'

  if(stringr::str_detect(outcome_string, valid_str_match)) {
    doses <- tox <- c()
    cohort_strs <- stringr::str_extract_all(
      outcome_string, cohort_str_match)[[1]]
    for(cohort_str in cohort_strs) {
      c_dl <- as.integer(stringr::str_extract(cohort_str, dl_str_match))
      if(c_dl <= 0) stop('Dose-levels must be strictly positive integers.')
      c_outcomes <- stringr::str_extract(cohort_str, outcomes_match_str)

      these_doses <- rep(c_dl, nchar(c_outcomes))
      doses <- c(doses, these_doses)

      these_outcomes = stringr::str_split(c_outcomes, '')[[1]]
      these_tox = as.integer((these_outcomes == 'T'))
      tox = c(tox, these_tox)
    }
  } else {
    stop(paste0('"', outcome_string, '" is not a valid outcome string.
                A valid example is "1N 2NN 3TT 2NT"'))
  }

  if(as.list) {
    return(list(
      doses = doses, tox = tox, num_patients = length(doses)
    ))
  } else {
    return(data.frame(doses = doses, tox = tox))
  }
}
