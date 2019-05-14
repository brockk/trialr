



test_that('efftox_dtps fails when cohort_sizes is not vector of +ve integers', {

  expect_error(
    efftox_dtps(cohort_sizes = c(3, 0),
                next_dose = 1,
                real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                p_e = 0.1, p_t = 0.1,
                eff0 = 0.5, tox1 = 0.65,
                eff_star = 0.7, tox_star = 0.25,
                alpha_mean = -7.9593, alpha_sd = 3.5487,
                beta_mean = 1.5482, beta_sd = 3.5018,
                gamma_mean = 0.7367, gamma_sd = 2.5423,
                zeta_mean = 3.4181, zeta_sd = 2.4406,
                eta_mean = 0, eta_sd = 0.2,
                psi_mean = 0, psi_sd = 1, seed = 123)
  )

  expect_error(
    efftox_dtps(cohort_sizes = c(3, -1),
                next_dose = 1,
                real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                p_e = 0.1, p_t = 0.1,
                eff0 = 0.5, tox1 = 0.65,
                eff_star = 0.7, tox_star = 0.25,
                alpha_mean = -7.9593, alpha_sd = 3.5487,
                beta_mean = 1.5482, beta_sd = 3.5018,
                gamma_mean = 0.7367, gamma_sd = 2.5423,
                zeta_mean = 3.4181, zeta_sd = 2.4406,
                eta_mean = 0, eta_sd = 0.2,
                psi_mean = 0, psi_sd = 1, seed = 123)
  )

  expect_error(
    efftox_dtps(cohort_sizes = c(3, 2.3),
                next_dose = 1,
                real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                p_e = 0.1, p_t = 0.1,
                eff0 = 0.5, tox1 = 0.65,
                eff_star = 0.7, tox_star = 0.25,
                alpha_mean = -7.9593, alpha_sd = 3.5487,
                beta_mean = 1.5482, beta_sd = 3.5018,
                gamma_mean = 0.7367, gamma_sd = 2.5423,
                zeta_mean = 3.4181, zeta_sd = 2.4406,
                eta_mean = 0, eta_sd = 0.2,
                psi_mean = 0, psi_sd = 1, seed = 123)
  )

  expect_error(
    efftox_dtps(cohort_sizes = c(3, NA),
                next_dose = 1,
                real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                p_e = 0.1, p_t = 0.1,
                eff0 = 0.5, tox1 = 0.65,
                eff_star = 0.7, tox_star = 0.25,
                alpha_mean = -7.9593, alpha_sd = 3.5487,
                beta_mean = 1.5482, beta_sd = 3.5018,
                gamma_mean = 0.7367, gamma_sd = 2.5423,
                zeta_mean = 3.4181, zeta_sd = 2.4406,
                eta_mean = 0, eta_sd = 0.2,
                psi_mean = 0, psi_sd = 1, seed = 123)
  )
})
