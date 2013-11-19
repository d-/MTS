\name{ECMxmtx}
\alias{ECMxmtx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix of Explanatory Variables of an ECM Model
%%  ~~function to do ... ~~
}
\description{Constructs the matrix of explanatory variables of an ECM model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ECMxmtx(x, p, m, include.con)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Matrix of multivariate time series data
%%     ~~Describe \code{x} here~~
}
  \item{p}{VAR order
%%     ~~Describe \code{p} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{include.con}{A logical switch for including the mean vector
%%     ~~Describe \code{include.con} here~~
}
}
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
\author{Ruey S Tsay
%%  ~~who you are~~
}
\note{An internal routine, not intended to be a general R command
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
function (x, p, m, include.con) 
{
    nT = dim(x)[1]
    k = dim(x)[2]
    dx = x[2:nT, ] - x[1:(nT - 1), ]
    dx = rbind(rep(0, k), dx)
    ist = p
    if (ist < 2) 
        ist = 2
    xm = x[(ist - 1):(nT - 1), ]
    ne = nT - ist + 1
    if (include.con) 
        xm = cbind(xm, rep(1, ne))
    if (p > 1) {
        for (ii in 1:(p - 1)) {
            xm = cbind(xm, dx[(ist - ii):(nT - ii), ])
        }
    }
    ECMxmtx <- list(xm = xm, y = dx[ist:nT, ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
