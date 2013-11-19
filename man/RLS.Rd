\name{RLS}
\alias{RLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursive Least squares
%%  ~~function to do ... ~~
}
\description{Compute recursive least squares estimation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RLS(y, x, ist = 30, xpxi = NULL, xpy0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data of dependent variable
%%     ~~Describe \code{y} here~~
}
  \item{x}{data matrix of regressors
%%     ~~Describe \code{x} here~~
}
  \item{ist}{initial number of data points used to start the estimation
%%     ~~Describe \code{ist} here~~
}
  \item{xpxi}{Inverse of the X'X matrix
%%     ~~Describe \code{xpxi} here~~
}
  \item{xpy0}{Initial value of X'y.
%%     ~~Describe \code{xpy0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Recursive least squares estimates
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
\note{This function is used internally, but can be used as a command too.
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
function (y, x, ist = 30, xpxi = NULL, xpy0 = NULL) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nT = dim(x)[1]
    k = dim(x)[2]
    if (is.null(xpxi)) {
        if (ist <= k) 
            ist = k + 1
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
        xpxi = xpxi - wk1
        Kmtx = xpxi \%*\% t(x1)
        tmp = beta[jdx, ] - t(Kmtx) * (sum(x1 * beta[jdx, ]) - 
            y[t])
        beta = rbind(beta, tmp)
        res = y[t] - tmp \%*\% t(x1)
        resi = c(resi, res)
    }
    beta = beta[-1, ]
    RLS <- list(beta = beta, resi = resi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
