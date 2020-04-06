// Dose-Finding Based on Efficacy-Toxicity Trade-Offs, by Thall, Cook
// Effective sample size for computing prior hyperparameters in Bayesian phase I-II dose-finding,
//  by Thall, Herrick, Nguyen, Venier, Norris

functions {
  real log_joint_pdf(real[] coded_doses, real[] coded_doses_squ,
                     int num_patients, int[] eff, int[] tox, int[] doses,
                     real alpha, real beta, real gamma, real zeta, real eta, real psi) {
    real p;
    p = 0;
    for(j in 1:num_patients) {
      real prob_eff;
      real prob_tox;
      real p_j;
      prob_eff = inv_logit(gamma + zeta * coded_doses[doses[j]] + eta * coded_doses_squ[doses[j]]);
      prob_tox = inv_logit(alpha + beta * coded_doses[doses[j]]);
      p_j = prob_eff^eff[j] * (1. - prob_eff)^(1. - eff[j]) * prob_tox^tox[j] *
              (1. - prob_tox)^(1. - tox[j]) + (-1.)^(eff[j] + tox[j]) * prob_eff *
              prob_tox * (1. - prob_eff) * (1. - prob_tox) *
              (exp(psi) - 1.) / (exp(psi) + 1.);
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
  real gamma_mean;
  real<lower=0> gamma_sd;
  real zeta_mean;
  real<lower=0> zeta_sd;
  real eta_mean;
  real<lower=0> eta_sd;
  real psi_mean;
  real<lower=0> psi_sd;
  // Fixed trial parameters
  int<lower=1> num_doses;
  real<lower=0> real_doses[num_doses]; // Doses under investigation, e.g. 10, 20, 30 for 10mg, 20mg, 30mg
  real p;  // The p of the L^p norm used to model the efficacy-toxicity indifference contours.
           // See Efficacy-Toxicity trade-offs based on L-p norms: Technical Report UTMDABTR-003-06
  real eff0; // Minimum required Pr(Efficacy) when Pr(Toxicity) = 0
  real tox1; // Maximum permissable Pr(Toxicity) when Pr(Efficacy) = 1
  real efficacy_hurdle; // A dose is acceptable if prob(eff) exceeds this hurdle...
  real toxicity_hurdle; //  ... and prob(tox) is less than this hurdle
  // Observed trial outcomes
  int<lower=0> num_patients;
  int<lower=0, upper=1> eff[num_patients]; // Binary efficacy event for patients j=1,..,num_patients
  int<lower=0, upper=1> tox[num_patients]; // Binary toxicity event for patients j=1,..,num_patients
  int<lower=1, upper=num_doses> doses[num_patients];  // Dose-levels given for patients j=1,..,num_patients.
                                   // Dose-levels are 1-based indices of real_doses
                                   // E.g. 1 means 1st dose in real_doses was given
}

transformed data {
  // Thall & Cook transform the actual doses by logging and centralising:
  real coded_doses[num_doses];
  real coded_doses_squ[num_doses]; // The square of coded_doses
  real mean_log_dose; // Variable created for convenience
  mean_log_dose = 0.0;
  for(i in 1:num_doses)
    mean_log_dose = mean_log_dose + log(real_doses[i]);
  mean_log_dose = mean_log_dose / num_doses;
  for(i in 1:num_doses)
  {
    coded_doses[i] = log(real_doses[i]) - mean_log_dose;
    coded_doses_squ[i] = coded_doses[i]^2;
  }
}

parameters {
  // Coefficients in toxicity logit model:
  real alpha;
  real beta;
  // Coefficients in efficacy logit model:
  // Why not delta and epsilon? They have such common, alternative usage that I skipped them.
  real gamma;
  real zeta;
  real eta;
  // Association:
  real psi;
}

transformed parameters {
  real<lower=0, upper=1> prob_eff[num_doses]; // Posterior probability of efficacy at doses i=1,...,num_doses
  real<lower=0, upper=1> prob_tox[num_doses]; // Posterior probability of toxicity at doses i=1,...,num_doses
  real utility[num_doses]; // Posterior utility of doses i=1,...,num_doses
  // Calculate the utility of each dose using the method described in
  // "Efficacy-Toxicity trade-offs based on L-p norms: Technical Report UTMDABTR-003-06", John Cook
  for(i in 1:num_doses)
  {
    real r_to_the_p; // Convenience variable, as in (2) of Cook.
    prob_tox[i] = inv_logit(alpha + beta * coded_doses[i]);
    prob_eff[i] = inv_logit(gamma + zeta * coded_doses[i] + eta * coded_doses_squ[i]);
    r_to_the_p = ((1 - prob_eff[i]) / (1 - eff0))^p + (prob_tox[i] / tox1)^p;
    utility[i] = 1 - r_to_the_p ^ (1.0/p);
  }
}

model {
  target += normal_lpdf(alpha | alpha_mean, alpha_sd);
  target += normal_lpdf(beta | beta_mean, beta_sd);
  target += normal_lpdf(gamma | gamma_mean, gamma_sd);
  target += normal_lpdf(zeta | zeta_mean, zeta_sd);
  target += normal_lpdf(eta | eta_mean, eta_sd);
  target += normal_lpdf(psi | psi_mean, psi_sd);
  target += log_joint_pdf(coded_doses, coded_doses_squ, num_patients, eff, tox, doses,
                          alpha, beta, gamma, zeta, eta, psi);
}

