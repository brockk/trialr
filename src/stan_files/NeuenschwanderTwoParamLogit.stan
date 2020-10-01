// The two-parameter logit model for dose-finding presented by Neuenschwander,
// Branson and Gsponer (2009), using logistic link function and normal priors on
// the intercept and slope terms.
//
// i.e. F(x, alpha, beta) = 1 / (1 + exp{-(alpha + exp(beta) log(x / x*))})
// where x is dose, x* is a reference dose, the intercept parameter is
// alpha ~ N(alpha_mean, alpha_sd)
// and the slope parameter is
// beta ~ N(beta_mean, beta_sd).
//
// The slope is exponentiated because it is required to be positive for
// monotonic dose-toxicity curves. However, the intercept can take any value so
// it is not exponentiated. In their paper, the authors use the log of the
// intercept and gradient terms presented above. I have done that to promote
// consistency with models already implemented in trialr.
//
// This model is similar to the 2-param logit CRM. It differs in that it uses
// dose as an explanatory covariate, rather than a skeleton.
//
// References:
// 1) Neuenschwander, B., Branson, M., & Gsponer, T. (2008).
// Critical aspects of the Bayesian approach to phase I cancer trials.
// Statistics in Medicine, 27, 2420–2439. https://doi.org/10.1002/sim

functions {
  real log_joint_pdf(int num_patients, int[] tox, real[] codified_doses,
                     real[] weights, real alpha, real beta) {
    real p;
    p = 0;
    for(j in 1:num_patients) {
      real prob_tox;
      real p_j;
      prob_tox = inv_logit(alpha + exp(beta) * codified_doses[j]);
      p_j = (weights[j] * prob_tox)^tox[j] *
              (1 - weights[j] * prob_tox)^(1 - tox[j]);
      p += log(p_j);
    }
    return p;
  }
}

data {
  // Hyperparameters
  real alpha_mean;
  real<lower=0> alpha_sd;
  real beta_mean;
  real<lower=0> beta_sd;

  // Fixed trial parameters
  int<lower=1> num_doses;

  // Doses under investigation, e.g. 10, 20, 30 for 10mg, 20mg, 30mg:
  real<lower=0> real_doses[num_doses];
  // Reference dose
  real<lower=0> d_star;

  // Observed trial outcomes
  int<lower=0> num_patients;
  // Binary toxicity event for patients j=1,..,num_patients
  int<lower=0, upper=1> tox[num_patients];
  // Dose-levels given for patients j=1,..,num_patients.
  // Dose-levels are 1-based indices of real_doses.
  // E.g. 1 means 1st dose in real_doses was given.
  int<lower=1, upper=num_doses> doses[num_patients];
  // Weights given to observations j=1, .., num_patients.
  // Weights between 0 and 1 suggest use of TITE-CRM.
  // However, weights can take whatever real value you want.
  real weights[num_patients];
}

transformed data {
  // Codified doses in this model are simply doses given over the reference dose
  real codified_doses[num_patients];
  for(j in 1:num_patients) {
    codified_doses[j] = log(real_doses[doses[j]] / d_star);
  }
}

parameters {
  // Parameters in toxicity model:
  real alpha;
  real beta;
}

transformed parameters {
  // Posterior probability of toxicity at doses i=1,...,num_doses
  real<lower=0, upper=1> prob_tox[num_doses];
  for(i in 1:num_doses) {
    prob_tox[i] = inv_logit(alpha + exp(beta) * log(real_doses[i] / d_star));
  }
}

model {
  target += normal_lpdf(alpha | alpha_mean, alpha_sd);
  target += normal_lpdf(beta | beta_mean, beta_sd);
  target += log_joint_pdf(num_patients, tox, codified_doses, weights,
                          alpha, beta);
}

generated quantities {
  vector[num_patients] log_lik;
  for (j in 1:num_patients) {
    real p_j;
    p_j = inv_logit(alpha + exp(beta) * codified_doses[j]);
    log_lik[j] = log((weights[j] * p_j)^tox[j] *
                  (1 - weights[j] * p_j)^(1 - tox[j]));
  }
}
