
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

  cohorts <- parse_eff_tox_dose_finding_outcomes(outcome_string)
  doses = integer(length = 0)
  eff = integer(length = 0)
  tox = integer(length = 0)
  for(cohort in cohorts) {
    c_dl <- cohort$dose
    c_outcomes <- cohort$outcomes

    these_outcomes <- stringr::str_split(c_outcomes, '')[[1]]
    these_eff = as.integer((these_outcomes == 'E') | (these_outcomes == 'B'))
    these_tox = as.integer((these_outcomes == 'T') | (these_outcomes == 'B'))
    these_doses <- rep(c_dl, length(these_tox))

    doses <- c(doses, these_doses)
    eff = c(eff, these_eff)
    tox = c(tox, these_tox)
  }

  if(as.list) {
    return(list(
      doses = doses, eff = eff, tox = tox, num_patients = length(doses)
    ))
  } else {
    return(data.frame(doses = doses, eff = eff, tox = tox))
  }
}

