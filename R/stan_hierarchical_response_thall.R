
#' Fit the hierarchical response model described by
#' Thall \emph{et al.} (2003).
#'
#' Fit the hierarchical response model to exchangeable groups described by
#' Thall \emph{et al.} (2003).
#'
#' @param group_responses vector of integers, number of responses in each group
#' @param group_sizes vector of integers, number of patients in each group
#' @param mu_mean mean parameter of normal prior distribution on mu. See details.
#' @param mu_sd standard deviation parameter of normal prior distribution on mu.
#' See details.
#' @param tau_alpha parameter alpha of inverse gamma prior distribution on tau. See
#' details.
#' @param tau_beta beta parameter of inverse gamma prior distribution on tau. See
#' details.
#' @param ... Extra parameters are passed to \code{rstan::sampling}. Commonly
#' used options are \code{iter}, \code{chains}, \code{warmup}, \code{cores}, and
#' \code{control}.
#'
#' @return Object of class \code{\link[rstan:stanfit]{rstan::stanfit}} returned
#' by \code{\link[rstan:sampling]{rstan::sampling}}
#'
#' @details
#' Thall \emph{et al.} (2003) describe hierarchical methods for analysing
#' treatment effects of a common intervention in several sub-types of a disease.
#' The treatment effects are assumed to be different but exchangeable and
#' correlated. Observing efficacy in one cohort, for example, increases one's
#' expectations of efficacy in others.
#' They demonstrate the hierarchical approach in a trial with binary response
#' outcomes and in another with time-to-event outcomes.
#' This function fits their model for binary response outcomes.
#'
#' Let the probability of response in group \eqn{i} be \eqn{\pi[i]} for
#' \eqn{i = 1,...,N}.
#' They assume a logistic model so that
#' \eqn{\theta_{i} = \log{\pi_{i} / (1 - \pi_{i})}}
#' is the log-odds of response in group \eqn{i}.
#' They assume that \eqn{\theta_{i} \sim N(\mu, \sigma^2)}.
#'
#' The authors implemented their model in BUGS.
#' As is the convention in BUGS, the authors define normal distributions by a
#' precision parameter \eqn{\tau} as opposed to the standard deviation parameter
#' \eqn{\sigma} used here. We have re-specified their model to comply with the
#' Stan convention of using standard deviation. The authors use a normal
#' prior on \eqn{\mu}, and a gamma prior on \eqn{\tau}, equivalent to
#' an inverse gamma prior on \eqn{\tau^{-1} = \sigma^2}.
#'
#' The authors provide WinBUGS code in their publication.
#' We implement their model here in Stan.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example from p.778 of Thall et al. (2003)
#' mod0 <- stan_hierarchical_response_thall(
#'   group_responses = c(0, 0, 1, 3, 5, 0, 1, 2, 0, 0),
#'   group_sizes = c(0, 2 ,1, 7, 5, 0, 2, 3, 1, 0),
#'   mu_mean = -1.3863,
#'   mu_sd = sqrt(1 / 0.1),
#'   tau_alpha = 2,
#'   tau_beta = 20)
#' }
#'
#' @references Thall, Wathen, Bekele, Champlin, Baker, and Benjamin. 2003.
#' “Hierarchical Bayesian approaches to phase II trials in diseases with
#' multiple subtypes.” Statistics in Medicine 22 (5): 763–80.
#' https://doi.org/10.1002/sim.1399.
#'
#' @seealso
#'   \code{\link[rstan:stanfit]{rstan::stanfit}},
#'   \code{\link[rstan:sampling]{rstan::sampling}}
#'
stan_hierarchical_response_thall <- function(
  group_responses, group_sizes, mu_mean, mu_sd, tau_alpha, tau_beta,...) {

  if(length(group_sizes) != length(group_responses))
    stop('group_sizes and group_responses must be the same length.')
  if(any(group_responses > group_sizes))
    stop('Cannot have more responses than patients in any group.')

  num_groups <- length(group_sizes)
  dat <- list(
    num_groups = num_groups,
    group_responses = group_responses,
    group_sizes = group_sizes,
    mu_mean = mu_mean,
    mu_sd = mu_sd,
    tau_alpha = tau_alpha,
    tau_beta = tau_beta
  )
  fit <- rstan::sampling(stanmodels$ThallHierarchicalBinary,
                         data = dat, ...)
  fit
}
