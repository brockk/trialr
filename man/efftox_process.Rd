% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efftox_process.R
\name{efftox_process}
\alias{efftox_process}
\title{Process RStan samples from an EffTox model}
\usage{
efftox_process(dat, fit)
}
\arguments{
\item{dat}{An instance of \code{\link{efftox_params}}, a list of EffTox
parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.}

\item{fit}{An instance of \code{rstan::stanmodel}, derived by fitting the
trialr EffTox model.}
}
\value{
An instance of \code{\link{efftox_fit}}.
}
\description{
Internal function to process rstan samples from an EffTox model
to make inferences about dose-acceptability, dose-utility and which dose
should be recommended next.
}
