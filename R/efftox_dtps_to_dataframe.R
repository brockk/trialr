
#' @title Calculate dose-transition pathways for an EffTox study
#'
#' @description Calculate dose-transition pathways for an EffTox study.
#' Note that TODO TODO TODO
#'
#' @param dat An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param cohort_sizes vector of future cohort sizes, i.e. positive integers.
#' E.g. To calculate paths for the the next cohort of two followed by another
#' cohort of three, use \code{cohort_sizes = c(2, 3)}.
#' @param next_dose the dose-level to be given to the immediately next cohort.
#' @param ... extra params passed to \code{rstan::sampling}.
#'
#' @return dose pathways in a \code{data.frame}.
#'
#' @export
#'
#' @examples
#' # Calculate the paths for the first cohort of 3 in Thall et al 2014 example
#' dat <- efftox_parameters_demo()
#' \dontrun{
#' dtps1 <- efftox_dtps_to_dataframe(dat = dat, cohort_sizes = c(3),
#'                                   next_dose = 1)
#' }
#' # To calculate future paths in a partially-observed trial
#' dat <- efftox_parameters_demo()
#' dat$doses = array(c(1,1,1))
#' dat$eff = array(c(0,0,0))
#' dat$tox = array(c(1,1,1))
#' dat$num_patients = 3
#' \dontrun{
#' dtps2 <- efftox_dtps_to_dataframe(dat = dat, cohort_sizes = c(3),
#'                                   next_dose = 1)
#' }
#'
#' @seealso
#' \code{\link{efftox_dtps}}, \code{\link{efftox_params}},
#' \code{\link{efftox_parameters_demo}}
#'
#' @references
#' Brock K, Billingham L, Copland M, Siddique S, Sirovica M, Yap C.
#' Implementing the EffTox dose-finding design in the Matchpoint trial.
#' BMC Medical Research Methodology. 2017;17(1):112.
#' doi:10.1186/s12874-017-0381-x
#'
#'
efftox_dtps_to_dataframe <- function(dat, cohort_sizes, next_dose, ...) {

  if(!all(cohort_sizes == ceiling(cohort_sizes)))
    stop('cohort_sizes must be stricly positive integers.')
  if(!all(cohort_sizes > 0))
    stop('cohort_sizes must be stricly positive integers.')

  previous_doses = dat$doses
  previous_eff = dat$eff
  previous_tox = dat$tox
  previous_num_patients = dat$num_patients
  outcomes <- c('E', 'T', 'N', 'B')

  # Calculate feasible outcome combinations by cohort
  cohort_paths <- lapply(cohort_sizes,
                         function(x) gtools::combinations(n = 4, r = x,
                                                          v = outcomes,
                                                          repeats.allowed=TRUE))
  # Flatten cohort outcomes
  cohort_paths <- lapply(cohort_paths, function(x) apply(x, 1, paste0,
                                                         collapse = ''))

  # Calculate pathways
  cohort_paths <- expand.grid(cohort_paths, stringsAsFactors = FALSE)
  # Place to record dose recommendations
  dose_recs = matrix(nrow = nrow(cohort_paths), ncol = ncol(cohort_paths))
  # Cache DTP calculations to avoid needless repetition
  cache <- new.env()
  for(i in 1:nrow(cohort_paths)) {
    cohort_path <- cohort_paths[i,]
    cohort_dose <- next_dose
    dtp <- ""

    for(j in 1:length(cohort_path)) {
      dtp <- ifelse(nchar(dtp) > 0,
                    paste0(dtp, ' ', cohort_dose, cohort_path[j]),
                    paste0(cohort_dose, cohort_path[j])
      )
      if(dtp %in% names(cache)) {
        # Fetch from cache
        print(paste0('Fetching ', dtp, ' from cache'))
        cohort_dose <- cache[[dtp]]
      } else {
        # Calculate
        these_outcomes <- efftox_parse_outcomes(dtp)
        dat$doses <- array(c(previous_doses, these_outcomes$doses))
        dat$eff <- array(c(previous_eff, these_outcomes$eff))
        dat$tox <- array(c(previous_tox, these_outcomes$tox))
        dat$num_patients <- previous_num_patients +
          these_outcomes$num_patients
        print(paste0('Running ', dtp))
        fit <- rstan::sampling(stanmodels$EffTox, data = dat, ...)
        decision <- efftox_process(dat, fit)
        cohort_dose <- decision$recommended_dose
        # Cache
        cache[[dtp]] <- cohort_dose
      }
      dose_recs[i, j] <- cohort_dose
    }
  }

  df <- data.frame(D0 = rep(next_dose, nrow(cohort_paths)))
  for(k in 1:ncol(cohort_paths)) {
    df[, paste0('C', k - 1)] = cohort_paths[, k]
    df[, paste0('D', k)] = dose_recs[, k]
  }
  return(df)
}
