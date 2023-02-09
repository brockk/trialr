// The two-parameter Bayesian Continual Reassessment Method (CRM) model for
// dose-finding using the logistic link function and a normal prior on the
// parameters.
//
// i.e. F(x, alpha, beta) = 1 / (1 + exp{-(alpha + exp(beta) x)})
// where x is the skeleton (the prior dose-toxicity curve), the intercept
// parameter is
// alpha ~ N(alpha_mean, alpha_sd)
// and the slope parameter is
// beta ~ N(beta_mean, beta_sd).
//
// The slope is exponentiated because it is required to be positive for
// monotonic dose-toxicity curves. See p.18 Cheung (2011).
// However, the intercept can take any value so it is not exponentiated.
//
// References:
// 1) O’Quigley, J., Pepe, M., & Fisher, L. (1990).
// Continual reassessment method: a practical design for phase 1 clinical trials
// in cancer. Biometrics, 46(1), 33–48. http://doi.org/10.2307/2531628
// 2) Cheung, Y. K. (2011).
// Dose Finding by the Continual Reassessment Method.
// New York: Chapman & Hall / CRC Press.

functions {
  real log_joint_pdf(int num_patients, int[] tox, int[] doses, real[] weights,
                     real[] codified_doses, real alpha, real beta) {
    real p;
    p = 0;
    for(j in 1:num_patients) {
      real prob_tox;
      real p_j;
      prob_tox = inv_logit(alpha + exp(beta) * codified_doses[doses[j]]);
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

  // Prior probability of toxicity at each dose, commonly referred to as the
  // skeleton. Should be monotonically increasing.
  real<lower=0, upper = 1> skeleton[num_doses];

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
  // Codified doses, aka dose-labels, are obtained by backwards substitution
  // of skeleton values into dose-toxicity relationship. I.e. given assumed
  // dose-toxicity relationship and parameters, what "dose" gives me prob_tox =
  // skeleton[1], skeleton[2], etc? See p.18 of Cheung (2011).
  real codified_doses[num_doses];
  for(i in 1:num_doses) {
    codified_doses[i] = (logit(skeleton[i]) - alpha_mean) / exp(beta_mean);
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
    prob_tox[i] = inv_logit(alpha + exp(beta) * codified_doses[i]);
  }
}

model {
  target += normal_lpdf(alpha | alpha_mean, alpha_sd);
  target += normal_lpdf(beta | beta_mean, beta_sd);
  target += log_joint_pdf(num_patients, tox, doses, weights, codified_doses,
                          alpha, beta);
}

generated quantities {
  vector[num_patients] log_lik;
  for (j in 1:num_patients) {
    real p_j;
    p_j = inv_logit(alpha + exp(beta) * codified_doses[doses[j]]);
    log_lik[j] = log((weights[j] * p_j)^tox[j] *
                  (1 - weights[j] * p_j)^(1 - tox[j]));
  }
}
