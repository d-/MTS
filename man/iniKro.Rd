\name{iniKro}
\alias{iniKro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial Estimates for Hronecker Index Approach
%%  ~~function to do ... ~~
}
\description{Obtain initial parameter estimates of a VARMA model specified 
by the Kronecker index approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
iniKro(da, at, ARid, MAid, include.mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{at}{A T-by-k preliminary estimates of the innovations
%%     ~~Describe \code{at} here~~
}
  \item{ARid}{AR polynomial structure
%%     ~~Describe \code{ARid} here~~
}
  \item{MAid}{MA polynomial structure
%%     ~~Describe \code{MAid} here~~
}
  \item{include.mean}{A logical switch to include the mean vector
%%     ~~Describe \code{include.mean} here~~
}
}
\details{Initial function, not intended as a general command
%%  ~~ If necessary, more details than the description above ~~
}
\value{Initial parameter estimates
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 4)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{Kronfit
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, at, ARid, MAid, include.mean) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    if (!is.matrix(at)) 
        at = as.matrix(at)
    nT = dim(da)[1]
    k = dim(da)[2]
    p = floor(dim(ARid)[2]/k) - 1
    if (p <= 0) 
        p = 1
    ist = p + 1
    est = NULL
    estse = NULL
    for (i in 1:k) {
        X = NULL
        Y = da[ist:nT, i]
        if (include.mean) 
            X = rep(1, (nT - p))
        if (i > 1) {
            for (j in 1:(i - 1)) {
                if (ARid[i, j] > 1) {
                  tmp = at[ist:nT, j] - da[ist:nT, j]
                  X = cbind(X, tmp)
                }
            }
        }
        for (lag in 1:p) {
            jst = lag * k
            for (j in 1:k) {
                if (ARid[i, jst + j] > 1) {
                  tmp = da[(ist - lag):(nT - lag), j]
                  X = cbind(X, tmp)
                }
            }
        }
        for (lag in 1:p) {
            jst = lag * k
            for (j in 1:k) {
                if (MAid[i, jst + j] > 1) {
                  tmp = at[(ist - lag):(nT - lag), j]
                  X = cbind(X, tmp)
                }
            }
        }
        XPX = crossprod(X, X)/nT
        XPXinv = solve(XPX)
        XPY = crossprod(X, Y)/nT
        beta = XPXinv \%*\% XPY
        l1 = dim(XPX)[1]
        resi = Y - X \%*\% matrix(beta, l1, 1)
        evar = crossprod(resi, resi)/(nT - p)
        est = c(est, beta)
        estse = c(estse, sqrt(diag(XPXinv) * evar/nT))
    }
    iniKro <- list(par = est, se = estse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
