% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EffTox.R
\name{as.mcmc.list.efftox_fit}
\alias{as.mcmc.list.efftox_fit}
\title{Convert \code{\link{efftox_fit}} to instance of
\code{\link[coda]{mcmc.list}}}
\usage{
\method{as.mcmc.list}{efftox_fit}(x, ...)
}
\arguments{
\item{x}{Object of class \code{\link{efftox_fit}}}

\item{...}{Extra variables that are passed onwards.}
}
\value{
Object of class \code{\link[coda]{mcmc.list}}
}
\description{
This function allows trialr to use tidybayes functions.
}
