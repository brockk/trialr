\name{thallhierarchicalbinary_params}
\alias{thallhierarchicalbinary_params}
\title{Parameters to be passed to the ThallHierachicalBinary model in Stan}
\description{
The parameters required by the \code{trialr} implementation of Thall et al's
hierarchical Bayesian model for binary outcomes are bundled up in a list, for
convenience.
The model analyses the binary response rates to a common drug of multiple
subtypes of a single disease.
See the included vignette for a full description and demonstration.
}
\details{The list object has the following elements:
  \itemize{
   \item \code{m}, the number of disease sub-types under investigation,
   an integer.
   \item \code{x}, the counts of responses in the \code{m} cohorts,
   a vector of non-negative integers.
   \item \code{n}, the counts of patients in the \code{m} cohorts,
   a vector of non-negative integers.
   \item \code{target_resp} the threshold response rate sought in the cohorts,
   a number between 0 and 1.
   \item \code{mu_mean}, mean of the normal prior distribution on the treatment
   effects.
   \item \code{mu_sd}, standard deviation of the normal prior distribution
   on the treatment effects.
   \item \code{tau_alpha}, alpha parameter of the gamma prior distribution on
   the precision, tau.
   \item \code{tau_beta}, beta parameter of the gamma prior distribution on
   the precision, tau.
  }

  An example of the thallhierarchicalbinary_params list is returned by
  \code{\link{thallhierarchicalbinary_parameters_demo}}.
}

\value{
A list with elements identified above.
}

\seealso{
  \code{\link{thallhierarchicalbinary_parameters_demo}}
}
