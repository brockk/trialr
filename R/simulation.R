

.run.sims <- function(model, num_sims, sample_data_func, summarise_func, ...) {
  sims <- list()
  for(i in 1:num_sims) {
    print(i)
    dat <- sample_data_func()
    fit <- rstan::sampling(model, data = dat, ...)
    sim <- summarise_func(dat, fit)
    sims[[i]] <- sim
  }
  return(sims)
}
