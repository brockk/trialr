\name{model_ThallHierarchicalBinary}
\alias{model_ThallHierarchicalBinary}
\alias{ThallHierarchicalBinary}
\title{Stan model for Thall et al.'s hierarchical Bayesian model for binary data}
\description{
This is the Stan implementation of Thall et al.'s hierarchical model for
analysing the response rate of a common drug in multiple sub-types of a single
disease.

The model is compiled when \code{trialr} is installed and is available at
run-time under \code{stanmodels$ThallHierarchicalBinary}.

The design was first published by Thall et al. in 2003.

See the Hierarchical Bayesian Response vignette for a detailed description of
the probability model and a demonstration of the implementation in \code{trial}.
}
\references{
Thall, Wathen, Bekele, Champlin, Baker, & Benjamin. 2003. "Hierarchical Bayesian
approaches to phase II trials in diseases with multiple subtypes".
Statistics in Medicine, 22(5), 763-780. https://doi.org/10.1002/sim.1399
}
\examples{
# Hierarchical model for responses in a disease with multiple subtypes
dat <- thallhierarchicalbinary_parameters_demo()
samp <- rstan::sampling(stanmodels$ThallHierarchicalBinary, data = dat)
rstan::plot(samp, pars = 'p') # Plot the modelled response rates in the subtypes
}
\seealso{
\code{\link{thallhierarchicalbinary_params}}

\code{\link{thallhierarchicalbinary_parameters_demo}}
}
