
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
#' @param as.list TRUE (the default) to return a \code{list};
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
#' x = df_parse_outcomes('1NNN 2NTN 3TTT')
#' x$num_patients  # 9
#' x$doses         # c(1, 1, 1, 2, 2, 2, 3, 3, 3)
#' x$tox           # c(0, 0, 0, 0, 1, 0, 1, 1, 1)
#' sum(x$tox)      # 4
#'
#' @references
#' Brock, K., Billingham, L., Copland, M., Siddique, S., Sirovica, M., & Yap, C.
#' (2017). Implementing the EffTox dose-finding design in the Matchpoint trial.
#' BMC Medical Research Methodology, 17(1), 112.
#' https://doi.org/10.1186/s12874-017-0381-x
#'
df_parse_outcomes <- function(outcome_string, as.list = TRUE) {

  cohorts <- parse_dose_finding_outcomes(outcome_string)
  doses = integer(length = 0)
  tox = integer(length = 0)
  for(cohort in cohorts) {
    c_dl <- cohort$dose
    c_outcomes <- cohort$outcomes

    these_outcomes <- stringr::str_split(c_outcomes, '')[[1]]
    these_tox = as.integer((these_outcomes == 'T'))
    these_doses <- rep(c_dl, length(these_tox))

    doses <- c(doses, these_doses)
    tox = c(tox, these_tox)
  }

  if(as.list) {
    return(list(
      doses = doses, tox = tox, num_patients = length(doses)
    ))
  } else {
    return(data.frame(doses = doses, tox = tox))
  }
}
