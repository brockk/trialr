\name{eight_schools_get_data}
\alias{eight_schools_get_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provide data for eight schools example
}
\description{
BLAH BLAH
BLAH BLAH
BLAH BLAH
}
\usage{
eight_schools_get_data()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    schools_dat <- list(J = 8, y = c(28, 8, -3, 7, -1, 1, 18, 
        12), sigma = c(15, 10, 16, 11, 9, 11, 10, 18))
    return(schools_dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
