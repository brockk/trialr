# from binarySimCLF
# https://github.com/cran/binarySimCLF/blob/master/R/ranBin2.R

# This fails:
# install.packages('binarySimCLF')
# So take the functions from archive tarball

`mardia` <-
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

`solve2` <-
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
      muij <- mardia(a, b, c);
    }
    return(muij);
  }

`ranBin2` <-
  function(nRep, u, psi)
  {
    if(nRep > 0) {
      u12 <- solve2(u[1], u[2], psi);
      y <- matrix( rep(-1, 2 * nRep), nrow=nRep );
      y[,1] <- ifelse(runif(nRep) <= u[1], 1, 0 );
      y[,2] <- y[,1] * (runif(nRep) <= u12 / u[1]) + (1 - y[,1]) *
        (runif(nRep) <= (u[2] - u12) / (1 - u[1]));
      return(y)
    } else {
      return(matrix(ncol = 2, nrow = 0))
    }
  }

# E.g.
# probs <- c(0.3, 0.8)
# s <- ranBin2(1000, probs, psi=0.2)
# cor(s)
# plot(jitter(s))
# colMeans(s)
