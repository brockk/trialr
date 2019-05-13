
#' Dose selection function that practices careful escalation.
#'
#' Dose selection function that avoids dose-skipping in escalation and
#' advocates stopping when there is sufficient evidence that the risk of
#' toxicity at a reference dose exceeds some threshold.
#'
#' @param dose_finding_fit Instance of \code{\link{dose_finding_fit}}.
#' @param tox_threshold numeric, the toxicity threshold.
#' @param certainty_threshold numeric, the required confidence that the risk of
#' toxicity exceeds `tox_threshold` to advocate stopping.
#' @param reference_dose the integer index of the reference dose. 1 by default,
#' i.e. the lowest dose-level.
#' @param start_dose the integer index of the desired starting dose. 1 by
#' default. This is required for the function to give the desired answer when
#' no patients have yet been treated.
#'
#' @return an integer dose-level
#' @export
#'
#' @examples
#' \dontrun{
#' # CRM example
#' fit <- stan_crm('1N 2N 3T', skeleton = c(0.1, 0.2, 0.35, 0.6),
#'                 target = 0.2, model = 'empiric', beta_sd = 1,
#'                 seed = 123)
#' }
careful_escalation <- function(dose_finding_fit, tox_threshold,
                               certainty_threshold,
                               reference_dose = 1,
                               start_dose = 1) {

  prob_too_toxic <- prob_tox_exceeds(dose_finding_fit, tox_threshold)
  if(prob_too_toxic[reference_dose] > certainty_threshold) {
    # Stop
    return(NA)
  } else {
    if(dose_finding_fit$num_patients > 0) {
      # Select recommended dose, without skipping untested doses in escalation
      max_dose <- max(dose_finding_fit$doses)
      return(min(dose_finding_fit$recommended_dose, max_dose + 1))
    } else {
      # The trial has not started so select start dose
      return(start_dose)
    }
  }
}
