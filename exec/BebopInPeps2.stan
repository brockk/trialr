// BEBOP is a phase II design that admits predictive baseline information to
// study co-primary efficacy and toxicity outcomes. It was developed by Brock
// for the PePS2 trial but it uses the probability model at the core of EffTox
// by Thall and Cook.

functions {
  real log_joint_pdf(int J, int[] eff, int[] tox, int[] x1, int[] x2, int[] x3,
                     real alpha, real beta, real gamma, real zeta, real lambda,
                     real psi) {
    real p;
    p = 0;
    for(j in 1:J) {

      real prob1;
      real prob2;
      real p_delta;
      prob1 = inv_logit(alpha + beta*x1[j] + gamma*x2[j] + zeta*x3[j]);
      prob2 = inv_logit(lambda);
      p_delta = prob1^eff[j] * (1. - prob1)^(1. - eff[j]) * prob2^tox[j] *
        (1. - prob2)^(1. - tox[j]) +
        (-1.)^(eff[j] +tox[j]) * prob1 * prob2 * (1. - prob1) * (1. - prob2) *
        (exp(psi) - 1.) / (exp(psi) + 1.);
      p += log(p_delta);
    }
    return p;
  }
}

data {
  int<lower=0> J; // number of obs
  int<lower=0, upper=1> eff[J]; // event one
  int<lower=0, upper=1> tox[J]; // event two
  int<lower=0, upper=1> x1[J];  // x1=1 if pre-treated, else 0
  int<lower=0, upper=1> x2[J];  // x2=1 if PD-L1 is low
  int<lower=0, upper=1> x3[J];  // x3=1 if PD-L1 is medium
                                // Implicityly, x2=0 & x3=0 if PD-L1 is high

  // Hyperparameters
  real alpha_mean;
  real<lower=0> alpha_sd;
  real beta_mean;
  real<lower=0> beta_sd;
  real gamma_mean;
  real<lower=0> gamma_sd;
  real zeta_mean;
  real<lower=0> zeta_sd;
  real lambda_mean;
  real<lower=0> lambda_sd;
  real psi_mean;
  real<lower=0> psi_sd;
}

parameters {
  real alpha;
  real beta;
  real gamma;
  real zeta;
  real lambda;
  real psi;
}

transformed parameters {
  real<lower=0, upper=1> prob_eff[6];
  real<lower=0, upper=1> prob_tox[6];
  prob_eff[1] = inv_logit(alpha + beta*0 + gamma*1 + zeta*0);
  prob_eff[2] = inv_logit(alpha + beta*0 + gamma*0 + zeta*1);
  prob_eff[3] = inv_logit(alpha + beta*0 + gamma*0 + zeta*0);
  prob_eff[4] = inv_logit(alpha + beta*1 + gamma*1 + zeta*0);
  prob_eff[5] = inv_logit(alpha + beta*1 + gamma*0 + zeta*1);
  prob_eff[6] = inv_logit(alpha + beta*1 + gamma*0 + zeta*0);
  prob_tox[1] = inv_logit(lambda);
  prob_tox[2] = inv_logit(lambda);
  prob_tox[3] = inv_logit(lambda);
  prob_tox[4] = inv_logit(lambda);
  prob_tox[5] = inv_logit(lambda);
  prob_tox[6] = inv_logit(lambda);
}

model {
  target += normal_lpdf(alpha | alpha_mean, alpha_sd);
  target += normal_lpdf(beta | beta_mean, beta_sd);
  target += normal_lpdf(gamma | gamma_mean, gamma_sd);
  target += normal_lpdf(zeta | zeta_mean, zeta_sd);
  target += normal_lpdf(lambda | lambda_mean, lambda_sd);
  target += normal_lpdf(psi | psi_mean, psi_sd);
  target += log_joint_pdf(J, eff, tox, x1, x2, x3, alpha, beta, gamma, zeta, lambda, psi);
}
