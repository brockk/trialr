

library("rstan")

# Prepare data
url <- "http://stat.columbia.edu/~gelman/arm/examples/arsenic/wells.dat"
wells <- read.table(url)
wells$dist100 <- with(wells, dist / 100)
X <- model.matrix(~ dist100 + arsenic, wells)
standata <- list(y = wells$switch, X = X, N = nrow(X), P = ncol(X))

model_code <- "data {
  int<lower=0> N;             // number of data points
int<lower=0> P;             // number of predictors (including intercept)
matrix[N,P] X;              // predictors (including 1s for intercept)
int<lower=0,upper=1> y[N];  // binary outcome
}
parameters {
vector[P] beta;
}
model {
beta ~ normal(0, 1);
y ~ bernoulli_logit(X * beta);
}
generated quantities {
vector[N] log_lik;
for (n in 1:N) {
// preferred Stan syntax as of version 2.10.0
log_lik[n] = bernoulli_logit_lpmf(y[n] | X[n] * beta);
}
}
"

fit_1 <- stan(model_code = model_code, data = standata)
print(fit_1, pars = "beta")
names(fit_1)

library(loo)
log_lik_1 <- extract_log_lik(fit_1)
loo_1 <- loo(log_lik_1)
print(loo_1)


