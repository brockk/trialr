data {
  # Observed data and fixed trial parameters
  #   (N.b. not 'parameter' in the Stan sense)
  int<lower=1> m;                         // Number of cohorts
  int<lower=0> x[m];                      // Number of responses by cohort
  int<lower=0> n[m];                      // Number of patients by cohort
  real<lower=0.0, upper=1.0> target_resp; // Target response rate
  // Hyperparameters for mu
  real mu_mean;
  real<lower=0.0> mu_sd;
  // Hyperparameters for tau
  real<lower=0> tau_alpha;
  real<lower=0> tau_beta;
}

parameters {
  real mu;                // Mean of treatment effects
  real<lower=0.0> sigma;  // StDev of treatment effects
  real rho[m];            // Log-odds of response by cohort
}

transformed parameters {
  real<lower=0.0, upper=1.0> p[m];  // Prob of response by cohort
  real<lower=0, upper=1> pg[m];     // Prob that response rate > target_resp
  for (i in 1:m)
  {
    p[i] = inv_logit(rho[i]);
    pg[i] = int_step(p[i] - target_resp);
  }
}

model
{
  mu ~ normal(mu_mean, mu_sd);
  sigma ~ inv_gamma(tau_alpha, tau_beta);
  rho ~ normal(mu, sigma);
  x ~ binomial_logit(n, rho);
}
