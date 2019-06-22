# Ben Bolker's rLKJ function from emdbook package, copied from Richard
# McElreath's rethinking package on GitHub.


#' Sample LKJ correlation matrices.
#'
#' This function was copied from Richard McElreath's rethinking package hosted
#' at https://github.com/rmcelreath/rethinking. In turn, he appears to have
#' copied it from Ben Bolker's rLKJ function from the emdbook package, although
#' I cannot find it there (else I would have imported it).
#'
#' @param n Number of matrices to sample.
#' @param K dimenstion of matrix to sample.
#' @param eta Distribution parameter
#' @return matrix
#'
#' @importFrom stats rbeta rnorm
#' @export
rlkjcorr <- function (n, K, eta = 1) {

  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))

  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)

      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])

      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}
