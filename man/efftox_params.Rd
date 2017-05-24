\name{efftox_params}
\alias{efftox_params}
\title{Parameters to be passed to EffTox model in Stan}
\description{
The parameters required by the EffTox Stan model are bundled up in a list, for
convenience.
}
\details{The list object has the following elements:
  \itemize{
   \item \code{num_doses}, an integer.
   \item \code{real_doses}, a vector of numbers.The doses under investigation.
   They should be ordered from lowest to highest and be in consistent units.
   E.g., to conduct a dose-finding trial of doses 10mg, 20mg and 50mg, use
   c(10, 20, 50).
   \item \code{efficacy_hurdle}, acceptable doses should have efficacy
   probability exceeding this value, a number between 0 and 1.
   \item \code{toxicity hurdle}, acceptable doses should have toxicity
   probability less than this value, a number between zero and 1.
   \item \code{p_e}, certainty required to infer a dose is acceptable with
   regards to being probably efficacious; a number between 0 and 1.
   \item \code{p_t}, certainty required to infer a dose is acceptable with
   regards to being probably tolerable; a number between 0 and 1.
   \item \code{p}, the p-index of EffTox utility contours.
   Use \code{efftox_solve_p}. A number.
   \item \code{eff0}, Efficacy probability required when toxicity is impossible;
    a number between 0 and 1.
   \item \code{tox1}, Toxicity probability permitted when efficacy is guaranteed;
    a number between 0 and 1.
   \item \code{eff_star}, Efficacy probability of a third point on
   the neutral utility curve that joins (\code{eff0}, 0) to (1, \code{tox1}).
   A number between 0 and 1.
   \item \code{tox_star}, Toxicity probability of a third point on
   the neutral utility curve that joins (\code{eff0}, 0) to (1, \code{tox1}).
   A number between 0 and 1.
   \item \code{alpha_mean}, The prior mean of the intercept term in the toxicity
   logit model. A number.
   \item \code{alpha_sd}, The prior standard deviation of the intercept term in
   the toxicity logit model. A number.
   \item \code{beta_mean}, The prior mean of the slope term in the toxicity
   logit model. A number.
   \item \code{beta_sd}, The prior standard deviation of the slope term in
   the toxicity logit model. A number.
   \item \code{gamma_mean}, The prior mean of the intercept term in the efficacy
   logit model. A number.
   \item \code{gamma_sd}, The prior standard deviation of the intercept term in
   the efficacy logit model. A number.
   \item \code{zeta_mean}, The prior mean of the slope term in the efficacy
   logit model. A number.
   \item \code{zeta_sd}, The prior standard deviation of the slope term in
   the efficacy logit model. A number.
   \item \code{eta_mean}, The prior mean of the squared term coefficient in the
   efficacy logit model. A number.
   \item \code{eta_sd}, The prior standard deviation of the squared term
   coefficient in the efficacy logit model. A number.
   \item \code{psi_mean}, The prior mean of the association term in the combined
   efficacy-toxicity model. A number.
   \item \code{psi_sd}, The prior standard deviation of the association term in
   the combined efficacy-toxicity model. A number.
   \item \code{eff}, A vector of efficacy outcomes for patients 1:num_patients,
   where 1=efficacy and 0=no efficacy.
   \item \code{tox}, A vector of toxicity outcomes for patients 1:num_patients,
   where 1=toxicity and 0=no toxicity.
   \item \code{doses}, A vector of dose-levels given to patients 1:num_patients,
   where 1=lowest dose, 2=second dose, etc.
   \item \code{num_patients}, number of patients, positive integer
  }

  An example of the efftox_params list is returned by
  \code{\link{efftox_parameters_demo}}.
}

\value{
A list with elements identified above.
}

\seealso{
  \code{\link{efftox_parameters_demo}}

  \code{\link{efftox_solve_p}}
}
