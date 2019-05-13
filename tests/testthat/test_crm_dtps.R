

test_that('crm_dtps fails when cohort_sizes is not vector of +ve integers', {

  target <- 0.4
  skeleton <- seq(0.1, 0.6, 0.1)

  expect_error(
    crm_dtps(skeleton,
             target,
             model = 'empiric',
             cohort_sizes = c(3, 3, 0),
             previous_outcomes = '',
             beta_sd = 1)
  )

  expect_error(
    crm_dtps(skeleton,
             target,
             model = 'empiric',
             cohort_sizes = c(3, 3, -1),
             previous_outcomes = '',
             beta_sd = 1)
  )

  expect_error(
    crm_dtps(skeleton,
             target,
             model = 'empiric',
             cohort_sizes = c(3, 3, 2.3),
             previous_outcomes = '',
             beta_sd = 1)
  )

  expect_error(
    crm_dtps(skeleton,
             target,
             model = 'empiric',
             cohort_sizes = c(3, 3, NA),
             previous_outcomes = '',
             beta_sd = 1)
  )
})
