
data {
  int<lower=1> N; // Number of patients
  // <lower=0> should be possible but it crashes when called with empty arrays
  vector[N] z0;  // Baseline tumour sizes
  vector[N] z1;  // Tumour sizes at assessment 1
  vector[N] z2;  // Tumour sizes at assessment 2

  int<lower=0, upper=1> d1[N]; // Non-shrinkage failures at assessment 1
  int<lower=0, upper=1> d2[N]; // Non-shrinkage failures at assessment 2

  // Hyperparameters
  // Tumour size model
  real alpha_mean;
  real<lower=0> alpha_sd;
  real beta_mean;
  real<lower=0> beta_sd;
  real gamma_mean;
  real<lower=0> gamma_sd;
  real sigma_mean;
  real<lower=0> sigma_sd;
  // LKJ prior parameter on the correlation matrix of log tumour sizes
  real omega_lkj_eta;
  // Non-shrinkage failure model
  real alpha_d1_mean;
  real<lower=0> alpha_d1_sd;
  real gamma_d1_mean;
  real<lower=0> gamma_d1_sd;
  real alpha_d2_mean;
  real<lower=0> alpha_d2_sd;
  real gamma_d2_mean;
  real<lower=0> gamma_d2_sd;
}

transformed data {
  vector[2] y[N]; // y, log tumour size ratios, has N rows and 2 columns
  for(i in 1:N) {
    y[i][1] = log(z1[i] / z0[i]);
    y[i][2] = log(z2[i] / z0[i]);
  }
}

parameters {
  real alpha;
  real beta;
  real gamma;
  corr_matrix[2] Omega;
  vector<lower=0>[2] sigma;
  real alphaD1;
  real gammaD1;
  real alphaD2;
  real gammaD2;
}

transformed parameters {
  vector[2] Mu[N];
  cov_matrix[2] Sigma;
  vector[2] ProbD[N];

  for(i in 1:N) {
    Mu[i][1] = alpha + gamma * z0[i];
    Mu[i][2] = beta + gamma * z0[i];
    ProbD[i][1] = inv_logit(alphaD1 + gammaD1 * z0[i]);
    ProbD[i][2] = inv_logit(alphaD2 + gammaD2 * z1[i]);
  }

  Sigma = quad_form_diag(Omega, sigma);
}

model {
  // Vector to accumulate patients non-shrinkage failure status through  time
  vector[N] d;

  // Priors on elements of tumour size model
  alpha ~ normal(alpha_mean, alpha_sd);
  beta ~ normal(beta_mean, beta_sd);
  gamma ~ normal(gamma_mean, gamma_sd);

  // Prior on the standard deviations of log-tumour size
  // sigma ~ cauchy(0, 1);
  sigma ~ normal(sigma_mean, sigma_sd);

  // LKJ prior on the correlation matrix  of log tumour sizes at times 1, 2,...
  Omega ~ lkj_corr(omega_lkj_eta);

  // Priors on elements in non-shrinkage failure model
  alphaD1 ~ normal(alpha_d1_mean, alpha_d1_sd);
  gammaD1 ~ normal(gamma_d1_mean, gamma_d1_sd);
  alphaD2 ~ normal(alpha_d2_mean, alpha_d2_sd);
  gammaD2 ~ normal(gamma_d2_mean, gamma_d2_sd);

  // Likelihood for tumour sizes
  y ~ multi_normal(Mu, Sigma);

  // Likelihood for non-shrinkage failures
  d = rep_vector(0, N);
  for(i in 1:N) {
    d1[i] ~ bernoulli(ProbD[i][1]);
    if(d1[i] == 1) d[i] = 1;
    if(d[i] == 0) d2[i] ~ bernoulli(ProbD[i][2]);
  }
}
