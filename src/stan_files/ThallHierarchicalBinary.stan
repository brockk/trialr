// A common treatment is tested in num_groups related cohorts.
// Response rates in the cohorts are different yet exchangeable and correlated.
// Groups log-odds of response are assumed to come from a common normal distribution.
// See Hierarchical Bayesian approaches to phase II trials in diseases with
//    multiple subtypes, by Thall, Wathen, Bekele, Champlin, Baker & Benjamin.
// K Brock, April 2017

data {
  // Observed data and trial constants
  int<lower=1> num_groups;                   // Number of cohorts
  int<lower=0> group_responses[num_groups];  // Number of responses by cohort
  int<lower=0> group_sizes[num_groups];      // Number of patients by cohort

  // Hyperparameters for mu
  real mu_mean;
  real<lower=0.0> mu_sd;

  // Hyperparameters for tau
  real<lower=0> tau_alpha;
  real<lower=0> tau_beta;
}

parameters {
  real mu;                // Mean of the groupwise log-odds of response
  real<lower=0.0> sigma2;  // Variance of the groupwise log-odds of response
  real rho[num_groups];   // Log-odds of response by cohort
}

transformed parameters {
  real<lower=0.0> sigma;
  real<lower=0.0, upper=1.0> prob_response[num_groups]; // Probability of response by cohort
  sigma = sqrt(sigma2);
  for (i in 1:num_groups)
  {
    prob_response[i] = inv_logit(rho[i]);
  }
}

model
{
  mu ~ normal(mu_mean, mu_sd);  // TODO use non-cent
  sigma2 ~ inv_gamma(tau_alpha, tau_beta);
  rho ~ normal(mu, sigma);
  group_responses ~ binomial_logit(group_sizes, rho);
}
