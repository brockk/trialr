#' Fit a Neuenschwander, Branson & Gsponer logit dose-finding model
#'
#' Fit Neuenschwander, Branson & Gsponer logit model for dose-finding using Stan
#' for full Bayesian inference.
#'
#' @param outcome_str A string representing the outcomes observed hitherto.
#' See \code{\link{df_parse_outcomes}} for a description of syntax and
#' examples. Alternatively, you may provide \code{doses_given} and \code{tox}
#' parameters. See Details.
#' @param real_doses A vector of numbers, the doses under investigation. They
#' should be ordered from lowest to highest and be in consistent units.
#' E.g. to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
#' c(10, 20, 50).
#' @param d_star d_star, numeric reference dose-level. The linear covariate
#' in this logit model is \code{dose / d_star}.
#' @param target the target toxicity probability, a number between 0 and 1.
#' @param alpha_mean Prior mean of intercept variable for normal prior.
#' See Details.
#' @param alpha_sd Prior standard deviation of intercept variable for normal prior.
#' See Details.
#' @param beta_mean Prior mean of gradient variable for normal prior.
#' See Details.
#' @param beta_sd Prior standard deviation of slope variable for normal prior.
#' See Details.
#' @param doses_given A optional vector of dose-levels given to patients
#' 1:num_patients, where 1=lowest dose, 2=second dose, etc. Only required when
#' \code{outcome_str} is not provided.
#' @param tox An optional vector of toxicity outcomes for patients
#' 1:num_patients, where 1=toxicity and 0=no toxicity. Only required when
#' \code{outcome_str} is not provided.
#' @param weights An optional vector of numeric weights for the observations
#' for patients 1:num_patients, thus facilitating a time-to-event (TITE) design.
#' Can be used with \code{outcome_str}, or with \code{doses_given} and
#' \code{tox}. It is generally tidier to specify \code{doses_given},
#' \code{tox} and \code{weights} when a TITE-analysis is desired.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores}, and
#' \code{control}.
#'
#' @details
#' The quickest and easiest way to fit this model to some observed outcomes
#' is to describe the outcomes using \pkg{trialr}'s syntax for dose-finding
#' outcomes. See \code{\link{df_parse_outcomes}} for full details and examples.
#'
#' The two-parameter model form is:
#'
#' \eqn{F(x_{i}, \alpha, \beta) = 1 / (1 + \exp{-(\alpha + \exp{(\beta)} log(x_i / d_star))}) }
#'
#' and the required parameters are:
#'
#' \itemize{
#'   \item \code{alpha_mean}
#'   \item \code{alpha_sd}
#'   \item \code{beta_mean}
#'   \item \code{beta_sd}
#' }
#'
#' @return An object of class \code{nbg_fit}, which inherits behaviour from
#' \code{\link{crm_fit}}.
#'
#' @author Kristian Brock \email{kristian.brock@gmail.com}
#'
#' @references
#'
#'   Neuenschwander, B., Branson, M., & Gsponer, T. (2008).
#'   Critical aspects of the Bayesian approach to phase I cancer trials.
#'   Statistics in Medicine, 27, 2420–2439. https://doi.org/10.1002/sim
#'
#' @seealso
#'   \code{\link{crm_fit}}
#'   \code{\link[rstan:sampling]{sampling}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Non-TITE example:
#' fit1 <- stan_nbg('1NNN 2NNN 3TTT', real_doses = c(10, 20, 50, 100, 200),
#'                  d_star = 200, target = 0.25,
#'                  alpha_mean = -1, alpha_sd = 2,
#'                  beta_mean = 0, beta_sd = 1,
#'                  seed = 123)
#' fit1$recommended_dose
#'
#' # The seed is passed to the Stan sampler. The usual Stan sampler params like
#' # cores, iter, chains etc are passed on too via the ellipsis operator.
#'
#' # TITE-CRM example
#' fit2 <-stan_nbg(real_doses = c(10, 20, 50, 100, 200), d_star = 200,
#'                 target = 0.25,
#'                 doses_given = c(3, 3, 3, 3),
#'                 tox = c(0, 0, 0, 0),
#'                 weights = c(73, 66, 35, 28) / 126,
#'                 alpha_mean = -1, alpha_sd = 2,
#'                 beta_mean = 0, beta_sd = 1,
#'                 seed = 123)
#' fit2$recommended_dose
#' }
stan_nbg <- function(outcome_str = NULL, real_doses, d_star, target,
                     alpha_mean = NULL, alpha_sd = NULL,
                     beta_mean = NULL, beta_sd = NULL,
                     doses_given = NULL,
                     tox = NULL,
                     weights = NULL,
                     ...) {

  # Model requires monotonically increasing doses.
  # Stop if this is not the case.
  if(length(real_doses) > 1) {
    if(any(real_doses[-1] <= real_doses[-length(real_doses)])) {
      stop('Real doses must be strictly monotonically-increasing.')
    }
  }

  # Create parameters object to pass to Stan
  version <- list(
    trialr = utils::packageVersion("trialr"),
    rstan = utils::packageVersion("rstan")
  )
  dat <- list(num_doses = length(real_doses),
              real_doses = real_doses,
              d_star = d_star,
              target = target,
              alpha_mean = alpha_mean,
              alpha_sd = alpha_sd,
              beta_mean = beta_mean,
              beta_sd = beta_sd,
              version = version)
  # Initialise with no patients observed
  dat$num_patients = 0
  dat$doses = integer(length = 0)
  dat$tox = integer(length = 0)
  dat$weights = numeric(length = 0)

  # Add outcomes
  if(is.null(outcome_str)) {
    if(length(doses_given) != length(tox))
      stop('doses_given and tox vectors should have same length')

    dat$doses <- array(doses_given)
    dat$tox <- array(tox)
    dat$num_patients <- length(doses_given)
  } else {
    outcomes_df <- df_parse_outcomes(outcome_str, as.list = TRUE)
    dat$num_patients <- outcomes_df$num_patients
    dat$doses <- array(outcomes_df$doses)
    dat$tox <- array(outcomes_df$tox)
  }
  # Add weights if specified; infer all to be 1 if not.
  if(is.null(weights))
    dat$weights <- array(rep(1, dat$num_patients))
  else
    dat$weights <- array(weights)

  # Check parameters
  if(is.null(alpha_mean)) stop('alpha_mean parameter must be specified.')
  if(is.null(alpha_sd)) stop('alpha_sd parameter must be specified.')
  if(alpha_sd <= 0) stop('alpha_sd parameter must be strictly positive.')
  if(is.null(beta_mean)) stop('beta_mean parameter must be specified.')
  if(is.null(beta_sd)) stop('beta_sd parameter must be specified.')
  if(beta_sd <= 0) stop('beta_sd parameter must be strictly positive.')

  # Fit data to model using Stan, after performing model-specific checks.
  samp <- rstan::sampling(stanmodels$NeuenschwanderTwoParamLogit,
                          data = dat, ...)

  # Create useful output from posterior samples.
  # In this regard, the model perfectly mimics the CRM:
  decision <- crm_process(dat, samp)
  # However, prepend type to allow specialisation of generics, if needed:
  class(decision) <- c("nbg_fit", class(decision))


  return(decision)
}
