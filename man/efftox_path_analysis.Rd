% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efftox_path_analysis.R
\name{efftox_path_analysis}
\alias{efftox_path_analysis}
\title{Fit an EffTox model to the incrementally observed outcomes on a trial pathway.}
\usage{
efftox_path_analysis(outcome_str, verbose = FALSE, ...)
}
\arguments{
\item{outcome_str}{A string representing the outcomes observed hitherto.
See \code{\link{efftox_parse_outcomes}} for a description of syntax and
examples. Alternatively, you may provide \code{doses_given} and \code{tox}
parameters. See Details.}

\item{verbose}{logical, TRUE to get log messages.}

\item{...}{All other parameters are passed to \code{\link{stan_efftox}}.}
}
\value{
A \code{\link{list}} of \code{\link{dose_finding_path_node}} objects.
}
\description{
Fit a EffTox model to the outcomes cumulatively observed at the end of each
cohort in a trial pathway. E.g. if the trial pathway is 1EN 2NN 3BT, we have
three cohorts of two patients. This function will fit the model to the
following four states: before any patients have been evaluated; after 1EN;
after 1EN 2NN; and finally after 1EN 2NN 3BT. This allows us to analyse how
the trial model is evolving in its estimation as trial data is accumulated.
}
\examples{
\dontrun{
# EffTox example
paths <- efftox_path_analysis(
  outcome_str = '1NNN 2NEN 3NEB',
  real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
  efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
  p_e = 0.1, p_t = 0.1,
  eff0 = 0.5, tox1 = 0.65,
  eff_star = 0.7, tox_star = 0.25,
  alpha_mean = -7.9593, alpha_sd = 3.5487,
  beta_mean = 1.5482, beta_sd = 3.5018,
  gamma_mean = 0.7367, gamma_sd = 2.5423,
  zeta_mean = 3.4181, zeta_sd = 2.4406,
  eta_mean = 0, eta_sd = 0.2,
  psi_mean = 0, psi_sd = 1, seed = 123, refresh = 0)

length(paths)  # 4
names(paths)[1]  # ""
names(paths)[2]  # "1NNN"
names(paths)[3]  # "1NNN 2NEN"
names(paths)[4]  # "1NNN 2NEN 3NEB"
# Each node is an analysis fit to the cumulative outcomes
# Converting to a tibble presents some nice tidyverse-related opportunities
library(tibble)
df <- as_tibble(paths)
df
}
}
\seealso{
\code{\link{efftox_parse_outcomes}},
  \code{\link{stan_efftox}},
  \code{\link{dose_finding_path_node}}
}
\author{
Kristian Brock
}
