# Test stan_crm

# Operational checks ----
test_that('stan_crm passes ellipsis variables to rstan::sampling', {
  x <- stan_crm(outcome_str = '1NNN 2TNT',
                skeleton = c(0.1, 0.25, 0.4, 0.6),
                target = 0.25,
                model = 'empiric',
                beta_sd = sqrt(1.34),
                iter = 500, chains = 2, seed = 123)
  df <- as.data.frame(x$fit)
  # Expect 2 * 500 / 2 post-warmup samples
  expect_equal(nrow(df), 500)
})

# Accuracy checks ----
# Published example of logistic model:
# p.21 of Cheung's "Dose Finding by the Continual Reassessment Method",
# ISBN 9781420091519
test_that('stan_crm output for logistic model matches DFCRM by Cheung', {

  # The target is the result of this:
  # fooB <- dfcrm::crm(prior = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
  #                    tox = c(0, 0, 1, 0, 0), lev = c(3, 5, 5, 3, 4),
  #                    model = "logistic", intcpt = 3)

  epsilon <- 0.03

  skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
  target <- 0.25
  a0 <- 3
  beta_mean <- 0
  beta_sd <- sqrt(1.34)

  x <- stan_crm(outcome_str = '3N 5N 5T 3N 4N',
                skeleton = skeleton, target = target,
                model = 'logistic', a0 = a0,
                beta_mean = beta_mean, beta_sd = beta_sd,
                seed = 123)
  expect_equal(x$recommended_dose, 4)
  beta_samp <- as.data.frame(x, pars = 'beta')
  # fooB$estimate
  expect_lt(abs(mean(beta_samp$beta) - 0.28), epsilon)

  # fooB$doses
  coded_doses <- (gtools::logit(skeleton) - a0) / exp(beta_mean)
  expect_true(all(abs(coded_doses - c(-5.94, -4.99, -4.10, -3.41, -2.80)) <
                    epsilon))

  # Not given in the book but continuing the checking via the dfcrm package with
  # fooB$post.var
  expect_lt(abs(var(beta_samp$beta) - 0.11), epsilon)

  # fooB$ptox
  beta_mean_post <- mean(beta_samp$beta)
  prob_tox_post <- gtools::inv.logit(a0 + exp(beta_mean_post) * coded_doses)
  expect_true(all(abs(prob_tox_post - c(0.01, 0.03, 0.08, 0.18, 0.33)) <
                    epsilon))
  # Note that the estimated posterior probability of toxicity from dfcrm
  # is calculated by plugging the posterior estimate of beta into the assumed
  # dose-tox function. This differs from that calculated by trialr because
  # trialr draws samples from the posterior distribution of prob_tox and
  # calculates the mean from there.
})



# Invoke some expected errors ----
test_that("stan_crm demands monotonically increasing skeleton", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.4),
                        target = 0.25,
                        model = 'empiric',
                        beta_sd = sqrt(1.34),
                        seed = 123))
})

test_that("stan_crm requires beta_sd for 'empiric' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                skeleton = c(0.1, 0.25, 0.4, 0.6),
                target = 0.25,
                model = 'empiric',
                seed = 123))
})

test_that("stan_crm requires a0 for 'logistic' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic',
                        beta_mean = 0,
                        beta_sd = sqrt(1.34),
                        seed = 123))
})

test_that("stan_crm requires beta_mean for 'logistic' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic',
                        a0 = 1,
                        beta_sd = sqrt(1.34),
                        seed = 123))
})

test_that("stan_crm requires beta_sd for 'logistic' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic',
                        a0 = 1,
                        beta_mean = 0,
                        seed = 123))
})


test_that("stan_crm requires a0 for 'logistic_gamma' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic_gamma',
                        beta_shape = 1,
                        beta_inverse_scale = 1,
                        seed = 123))
})

test_that("stan_crm requires beta_shape for 'logistic_gamma' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic_gamma',
                        a0 = 1,
                        beta_inverse_scale = 1,
                        seed = 123))
})

test_that("stan_crm requires beta_inverse_scale for 'logistic_gamma' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic_gamma',
                        a0 = 1,
                        beta_shape = 1,
                        seed = 123))
})

test_that("stan_crm requires alpha_mean for 'logistic2' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic2',
                        alpha_sd = 1,
                        beta_mean = 0,
                        beta_sd = sqrt(1.34),
                        seed = 123))
})

test_that("stan_crm requires alpha_sd for 'logistic2' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic2',
                        alpha_mean = 0,
                        beta_mean = 0,
                        beta_sd = sqrt(1.34),
                        seed = 123))
})

test_that("stan_crm requires beta_mean for 'logistic2' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic2',
                        alpha_mean = 0,
                        alpha_sd = 1,
                        beta_sd = sqrt(1.34),
                        seed = 123))
})

test_that("stan_crm requires beta_sd for 'logistic2' model", {
  expect_error(stan_crm(outcome_str = '1NNN 2TNT',
                        skeleton = c(0.1, 0.25, 0.4, 0.6),
                        target = 0.25,
                        model = 'logistic2',
                        alpha_mean = 0,
                        alpha_sd = 1,
                        beta_mean = 0,
                        seed = 123))
})

