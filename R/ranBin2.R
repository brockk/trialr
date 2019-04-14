# from binarySimCLF
# https://github.com/cran/binarySimCLF/blob/master/R/ranBin2.R

# This fails:
# install.packages('binarySimCLF')
# So I have taken the functions from archive tarball

.mardia <-
  function(a, b, c)
  {
    if (a == 0 )
    {
      return(-c / b);

    }
    k <- ifelse( b > 0, 1, ifelse(b < 0, -1, 0) );
    p <- -0.5 * (b + k*sqrt(b^2 - 4 * a * c));
    r1 <- p / a;
    r2 <- c / p;

    r <- ifelse( r2 > 0, r2, r1 );
    return(r);
  }

.solve2 <-
  function(mui, muj, psi)
  {
    if (psi == 1)
    {
      return(mui*muj);
    }
    else if (psi != 1)
    {
      a <- 1 - psi;
      b <- 1 - a*(mui + muj);
      c <- -psi*(mui * muj);
      muij <- .mardia(a, b, c);
    }
    return(muij);
  }

#' @title Sample pairs of correlated binary events
#'
#' @description This function is reproduced from the \code{binarySimCLF} package
#' on CRAN. The original package appears no longer to be maintained.
#' View the original source at:
#'  https://github.com/cran/binarySimCLF/blob/master/R/ranBin2.R
#'
#' @param nRep Number of simulated event pairs, positive integer.
#' @param u Mean event probabilities, expressed as a vector of length 2. E.g.
#' to simulate associated bivariate events with probabilities 80% and 30%, use
#' \code{u = c(0.8, 0.3)}.
#' @param psi Odds ratio, number. This parameter controls the strength of
#' association. Use \code{psi = 1} for no association. Values greater than 1
#' correspond to increasingly positive association between the two events,
#' and vice-versa.
#'
#' @return Matrix of events represented as 0s and 1s, with \code{nRep} rows
#' and 2 columns. The first column is the incidence of event 1.
#'
#' @export
#'
#' @examples
#' probs <- c(0.8, 0.3)
#' s <- ranBin2(1000, probs, psi=0.2)  # 1000 pairs of outcomes
#' cor(s)  # Negatively correlated because psi < 1
#' colMeans(s)  # Event rates as expected
ranBin2 <-
  function(nRep, u, psi)
  {
    if(nRep > 0) {
      u12 <- .solve2(u[1], u[2], psi);
      y <- matrix( rep(-1, 2 * nRep), nrow=nRep );
      y[,1] <- ifelse(stats::runif(nRep) <= u[1], 1, 0 );
      y[,2] <- y[,1] * (stats::runif(nRep) <= u12 / u[1]) + (1 - y[,1]) *
        (stats::runif(nRep) <= (u[2] - u12) / (1 - u[1]));
      return(y)
    } else {
      return(matrix(ncol = 2, nrow = 0))
    }
  }
