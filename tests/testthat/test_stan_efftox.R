# Test stan_efftox and stan_efftox_demo

# Operational checks ----
test_that('stan_efftox passes ellipsis variables to rstan::sampling', {
  l <- efftox_parameters_demo()
  x <- stan_efftox(outcome_str = '1NNN 2ENN',
                   real_doses = l$real_doses,
                   efficacy_hurdle = l$efficacy_hurdle,
                   toxicity_hurdle = l$toxicity_hurdle,
                   p_e = l$p_e,
                   p_t = l$p_t,
                   eff0 = l$eff0,
                   tox1 = l$tox1,
                   eff_star = l$eff_star,
                   tox_star = l$tox_star,
                   alpha_mean = l$alpha_mean, alpha_sd = l$alpha_sd,
                   beta_mean = l$beta_mean, beta_sd = l$beta_sd,
                   gamma_mean = l$gamma_mean, gamma_sd = l$gamma_sd,
                   zeta_mean = l$zeta_mean, zeta_sd = l$zeta_sd,
                   eta_mean = l$eta_mean, eta_sd = l$eta_sd,
                   psi_mean = l$psi_mean, psi_sd = l$psi_sd,
                   iter = 500, chains = 2, seed = 123)
  df <- as.data.frame(x$fit)
  # Expect 2 * 500 / 2 post-warmup samples
  expect_equal(nrow(df), 500)
})

# Accuracy checks ----
# Recreate some known (i.e. independently calculated) output:
test_that('Thall et al. (2014) model fits correctly to "1NNN 2ENN"', {
  mod <- stan_efftox_demo('1NNN 2ENN', seed = 123)

  # In each case, the expected values are taken from the MD Anderson app at
  # https://biostatistics.mdanderson.org/softwaredownload/SingleSoftware.aspx?Software_Id=2
  # N.b. their app uses numerical methods to approximate posterior parameters
  # (so there is variation) but does not allow a seed to be set for isolated
  # calculations (only simulations). Thus, the comparisons below are
  epsilon <- 0.03

  expect_true(all(abs(
    mod$prob_eff - c(0.05, 0.26, 0.72, 0.86, 0.91)) < epsilon))
  expect_true(all(abs(
    mod$prob_tox - c(0.01, 0.01, 0.02, 0.06, 0.12)) < epsilon))
  expect_true(all(abs(
    mod$utility - c(-0.91, -0.47, 0.42, 0.64, 0.64)) < epsilon))
  expect_true(all(abs(
    mod$prob_acc_eff - c(0.01, 0.13, 0.80, 0.91, 0.94)) < epsilon))
  expect_true(all(abs(
    mod$prob_acc_tox - c(1.00, 1.00, 0.98, 0.93, 0.86)) < epsilon))
})

test_that('Thall et al. (2014) model fits correctly to "1NNN 2ENN 3ETB"', {
  mod <- stan_efftox_demo('1NNN 2ENN 3ETB', seed = 123)

  # In each case, the expected values are taken from the MD Anderson app at
  # https://biostatistics.mdanderson.org/softwaredownload/SingleSoftware.aspx?Software_Id=2
  # N.b. their app uses numerical methods to approximate posterior parameters
  # (so there is variation) but does not allow a seed to be set for isolated
  # calculations (only simulations). Thus, the comparisons below are
  epsilon <- 0.03

  expect_true(all(abs(
    mod$prob_eff - c(0.06, 0.24, 0.71, 0.89, 0.94)) < epsilon))
  expect_true(all(abs(
    mod$prob_tox - c(0.02, 0.06, 0.41, 0.77, 0.87)) < epsilon))
  expect_true(all(abs(
    mod$utility - c(-0.92, -0.63, -0.24, -0.41, -0.47)) < epsilon))
  expect_true(all(abs(
    mod$prob_acc_eff - c(0.01, 0.07, 0.85, 0.97, 0.98)) < epsilon))
  expect_true(all(abs(
    mod$prob_acc_tox - c(1.00, 0.98, 0.36, 0.08, 0.05)) < epsilon))
})

test_that('EffTox fails on nonsense input string', {
  expect_error(stan_efftox_demo('1NZN 2ENN', seed = 123))
})


# Invoke some expected errors ----
# test_that('stan_efftox passes ellipsis variables to rstan::sampling', {
#   l <- efftox_parameters_demo()
#   expect_error(stan_efftox(outcome_str = '1NNN 2ENN',
#                    real_doses = l$real_doses,
#                    efficacy_hurdle = l$efficacy_hurdle,
#                    toxicity_hurdle = l$toxicity_hurdle,
#                    p_e = l$p_e,
#                    p_t = l$p_t,
#                    eff0 = l$eff0,
#                    tox1 = l$tox1,
#                    eff_star = l$eff_star,
#                    tox_star = l$tox_star,
#                    alpha_mean = l$alpha_mean, alpha_sd = l$alpha_sd,
#                    beta_mean = l$beta_mean, beta_sd = -10,
#                    gamma_mean = l$gamma_mean, gamma_sd = l$gamma_sd,
#                    zeta_mean = l$zeta_mean, zeta_sd = l$zeta_sd,
#                    eta_mean = l$eta_mean, eta_sd = l$eta_sd,
#                    psi_mean = l$psi_mean, psi_sd = l$psi_sd,
#                    seed = 123))
#
#   # This correctly throws an error but test_that stops for the error?!
# })
