\name{efftox_analysis}
\alias{efftox_analysis}
\title{Processed results of an EffTox dose-update analysis}
\description{
Processed results of an EffTox dose-update analysis, yielded by
\code{\link{efftox_process}}.
The model is invoked on data and RStan performs sampling from the posterior
distribution.
\code{\link{efftox_process}} then conducts the required processing to ultimately
arrive at a dose recommendation.
}
\details{
The object is essentially a list with slots:
\itemize{
   \item \code{dose_indices}, a vector of integers representing the dose-levels
   under consideration.
   \item \code{recommended_dose}, an integer representing the dose recommended
   for the next patient or cohort. TODO: what happens when stop recommended?
   \item \code{prob_eff}, the posterior mean probabilities of efficacy at doses
   1:n; a vector of numbers between 0 and 1.
   \item \code{prob_eff}, the posterior mean probabilities of toxicity at doses
   1:n; a vector of numbers between 0 and 1.
   \item \code{prob_acc_eff}, the posterior mean probabilities that efficacy at
   the doses is acceptable, i.e. that it exceeds the minimum acceptable efficacy
   threshold; a vector of numbers between 0 and 1.
   \item \code{prob_acc_eff}, the posterior mean probabilities that toxicity at
   the doses is acceptable, i.e. that it is less than the maximum toxicity
   threshold; a vector of numbers between 0 and 1.
   \item \code{utility}, the utilities of doses 1:n, calculated by plugging the
   posterior mean probabilities of efficacy and toxicity into the utility
   formula, as advocated by Thall & Cook. Contrast to \code{post_utility};
   a vector of numbers.
   \item \code{post_utility}, the posterior mean utilities of doses 1:n,
   calculated from the posterior distributions of the utilities. This is in
   contrast to \code{utility}, which uses plug-in posterior means of efficacy
   and toxicity, as advocated by Thall & Cook; a vector of numbers.
   \item \code{acceptable}, a vector of logical values to indicate whether
   doses 1:n are acceptable, according to the rules for acceptable efficacy &
   toxicity, and rules on not skipping untested doses.
   \item \code{fit}, An object of class \code{\link[rstan:stanfit]{stanfit}},        containing the posterior samples.
  }
}

\value{
A list with elements identified above.
}

\seealso{
  \code{\link{efftox_process}}
}
