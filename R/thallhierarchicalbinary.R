
#' @title Get parameters to run the demo of Thall Hierarchical Binary model
#'
#' @description Get parameters to run the demo of Thall Hierarchical Binary model.
#' These match those used to demonstrate EffTox in Thall et al. 2003.
#'
#' @return a \code{list} of parameters
#' @export
#'
#' @examples
#' dat <- thallhierarchicalbinary_parameters_demo()
#' names(dat)
#' dat$target_resp == 0.3
#'
#' @seealso \code{\link{efftox_params}}
#' @references Thall et al. 2014, Effective sample size for computing prior
#' hyperparameters in Bayesian phase I-II dose-finding
thallhierarchicalbinary_parameters_demo <- function() {
  x <- list(
    m = 10,
    x = c(0, 0, 1, 3, 5, 0, 1, 2, 0, 0),
    n = c(0, 2 ,1, 7, 5, 0, 2, 3, 1, 0),
    target_resp = 0.3,
    mu_mean = -1.3863,
    mu_sd = sqrt(1 / 0.1),
    tau_alpha = 2,
    tau_beta = 20
  )
  return(x)
}
