\name{peps2_params}
\alias{peps2_params}
\title{Parameters to be passed to \code{BebopInPeps2} model in Stan}
\description{
The parameters required by the \code{BebopInPeps2} Stan model are bundled up in
a list, for convenience. This list contains elements that we might call "data"
(e.g. trial outcomes), and "parameters" (e.g. modelling choices).
}
\details{The list object has the following elements:
  \itemize{
   \item \code{J}, Total number of patients used, positive integer.
   \item \code{eff}, a vector of logical values to represent the presence (=1)
   or absence (=0) of efficacy events in the patients.
   \item \code{tox}, a vector of logical values to represent the presence (=1)
   or absence (=0) of toxicity events in the patients.
   \item \code{x1}, a vector of logical values to represent whether the
   patients are pre-treated (=1) or treatment-naive (=0).
   \item \code{x2}, a vector of logical values to represent whether the
   patients are Low PD-L1 (=1) or not (=0).
   \item \code{x3}, a vector of logical values to represent whether the
   patients are Medium PD-L1 (=1) or not (=0). Logically, if \code{x2} and
   \code{x3} are both zero, then the patient is High PD-L1.
   \item \code{cohort_n}, vector of cohort sizes, provided for convenience.
   \item \code{cohort_eff}, vector of counts of efficacy events by cohort,
   provided for convenience.
   \item \code{cohort_tox}, vector of counts of toxicity events by cohort,
   provided for convenience.
   \item \code{alpha_mean}, The prior mean of alpha. Alpha is the efficacy model
   intercept.
   \item \code{alpha_sd}, The prior standard deviation of alpha. Alpha is the
   efficacy model  intercept.
   \item \code{beta_mean}, The prior mean of beta. Beta is the efficacy model
   term for being previously treated.
   \item \code{beta_sd}, The prior standard deviation of beta. Beta is the
   efficacy model term for being previously treated.
   \item \code{gamma_mean}, The prior mean of gamma. Gamma is the efficacy model
   term for being PD-L1 score = Low.
   \item \code{gamma_sd}, The prior standard deviation of gamma. Gamma is the
   efficacy model term for being PD-L1 score = Low.
   \item \code{zeta_mean}, The prior mean of zeta. Zeta is the efficacy model
   term for being PD-L1 score = Medium.
   \item \code{zeta_sd}, The prior standard deviation of zeta. Zeta is the
   efficacy model term for being PD-L1 score = Medium.
   \item \code{lambda_mean}, The prior mean of lambda. Lambda is the toxicity
   model intercept.
   \item \code{lambda_sd}, The prior standard deviation of lambda. Lambda is the
   toxicity model intercept.
   \item \code{psi_mean}, The prior mean of psi. Psi is the joint model
   association parameter.
   \item \code{psi_sd}, The prior standard deviation of psi. Psi is the joint
   model association parameter.
  }

  An example of the peps2_params list is returned by
  \code{\link{peps2_get_data}}.
}

\value{
A list with elements identified above.
}

\seealso{
  \code{\link{peps2_get_data}}
}
