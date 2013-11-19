\name{WLS}
\alias{WLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Least Squares 
%%  ~~function to do ... ~~
}
\description{Perform weighted least squares estimation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WLS(y, x, size, xpxi = NULL, xpy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Data vector of dependent variable
%%     ~~Describe \code{y} here~~
}
  \item{x}{Data matrix of explanatory variables
%%     ~~Describe \code{x} here~~
}
  \item{size}{Window size used to estimate the weights
%%     ~~Describe \code{size} here~~
}
  \item{xpxi}{Inverse of the X'X matrix
%%     ~~Describe \code{xpxi} here~~
}
  \item{xpy}{X'y matrix
%%     ~~Describe \code{xpy} here~~
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
\author{Ruey S. Tsay
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
function (y, x, size, xpxi = NULL, xpy = NULL) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nT = dim(x)[1]
    k = dim(x)[2]
    ist = size
    if (size <= k) 
        size = k + 1
    if (is.null(xpxi)) {
        ist = size
        xpx0 = t(x[1:ist, ]) \%*\% x[1:ist, ]
        xpxi = solve(xpx0)
        xpy0 = t(x[1:ist, ]) \%*\% matrix(y[1:ist], ist, 1)
    }
    ist = ist + 1
    resi = NULL
    beta = matrix(c(xpxi \%*\% xpy0), 1, k)
    for (t in ist:nT) {
        jdx = t - ist + 1
        x1 = matrix(x[t, ], 1, k)
        wk = x1 \%*\% xpxi
        k1 = 1 + sum(wk * x1)
        wk1 = t(wk) \%*\% wk/k1
        xpxt = xpxi - wk1
        res = sum(beta[jdx, ] * x1) - y[t]
        tmpt = beta[jdx, ] - (wk/k1) * res
        x2 = matrix(x[(t - size), ], 1, k)
        wk2 = x2 \%*\% xpxt
        k2 = sum(wk2 * x2) - 1
        r2 = sum(x2 * tmpt) - y[t - size]
        tmp = tmpt - (wk2/k2) * r2
        beta = rbind(beta, tmp)
        xpxi = xpxt - t(wk2) \%*\% wk2/k2
        res = y[t] - tmp \%*\% t(x1)
        resi = c(resi, res)
    }
    neb = dim(beta)[1]
    beta = beta[2:neb, ]
    WLS <- list(beta = beta, resi = resi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
