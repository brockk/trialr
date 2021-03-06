% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trialr_simulate.R
\name{trialr_simulate}
\alias{trialr_simulate}
\title{Run a simulation study.}
\usage{
trialr_simulate(
  N,
  get_data_func,
  fit_model_func,
  summarise_func = NULL,
  num_logs = 10,
  num_saves = NULL,
  save_func = NULL
)
}
\arguments{
\item{N}{integer, number of simulated iterations to run.}

\item{get_data_func}{Function that takes no parameters and returns a sampled
dataset to be analysed. I.e. the call signature is f().}

\item{fit_model_func}{Function that accepts the output of
\code{get_data_func} as the sole parameter and fits the model or performs the
analysis, returning an object of arbitrary type.}

\item{summarise_func}{Optional. If provided, this function should accept the
ouputs of \code{get_data_func} and \code{fit_model_func} as parameters 1 & 2
and perform some post-fit processing or simplification. The result of this
call is the output from iteration i. If omitted, the fit object from
\code{fit_model_func} is simply used as the output from iteration i.}

\item{num_logs}{Number of log messages to receive about progress. NULL to
suppress logging. E.g. if N=100 and num_logs=10, you will get log messages
when i=10, 20, 30, etc.}

\item{num_saves}{Number of interimittent saves to attempt. NULL to
suppress saving E.g. if N=100 and num_saves=10, the save_func delegate will
be called after iteration i=10, 20, 30, etc.}

\item{save_func}{Optional. Function that takes the interim list of simulated
objects as the sole parameter and saves them somehow. This, combined with
\code{num_saves}, allows periodic saving of in-progress results to avoid
complete data loss if the simulation study fails for some reason.}
}
\value{
\code{list} of length \code{N}. The items in the list are as returned
by \code{summarise_func} or \code{fit_model_func}.
}
\description{
This function is a fairly flexible way of running simulation studies in
trialr, and beyond. It essentially uses delegates to perform this pattern:
\preformatted{
for i in 1:N:
  data = get_data_func()
  fit = fit_model_func(data)
  if summarise_func is null:
    sims[i] = fit
  else
    sims[i] = summarise_func(data, fit)
  end
loop
return sims
}
}
\examples{
get_data_func <- function() {
  group_sizes <- rbinom(n = 5, size = 50, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
  group_responses <- rbinom(n = 5, size = group_sizes,
                            prob = c(0.2, 0.5, 0.2, 0.2, 0.2))
  list(
    group_responses = group_responses, group_sizes = group_sizes,
    mu_mean = gtools::logit(0.1), mu_sd = 1, tau_alpha = 2, tau_beta = 20
  )
}
fit_model_func <- function(data) {
  data <- append(data, list(refresh = 0))
  do.call(stan_hierarchical_response_thall, args = data)
}
summarise_func <- function(data, fit) {
  # Probability that estimate response rate exceeds 30\%
  unname(colMeans(as.data.frame(fit, 'prob_response') > 0.3))
}
\dontrun{
sims <- trialr_simulate(N = 20, get_data_func, fit_model_func, summarise_func)
# Posterior probabilities that the response rate in each cohort exceeds 30\%:
do.call(rbind, sims)
# Cohorts are in columns; simulated iterations are in rows.
}
}
