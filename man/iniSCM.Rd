\name{iniSCM}
\alias{iniSCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial Estimates of SCM Approach
%%  ~~function to do ... ~~
}
\description{Obtain initial parameter estimates of a model specified by the 
SCM approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
iniSCM(da, at, locTmtx, locAR, locMA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{at}{A T-by-k initial estimates of the innovations
%%     ~~Describe \code{at} here~~
}
  \item{locTmtx}{Locator of 1 for the transformation matrix
%%     ~~Describe \code{locTmtx} here~~
}
  \item{locAR}{Locator of AR polynomial structure
%%     ~~Describe \code{locAR} here~~
}
  \item{locMA}{Locator of MA polynomial structure
%%     ~~Describe \code{locMA} here~~
}
}
\details{Initial function used by SCMfit, not intended as a general 
command
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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

\seealso{SCMfit
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, at, locTmtx, locAR, locMA) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    if (!is.matrix(at)) 
        at = as.matrix(at)
    inc.mean <- inc.meanSCM
    nT = dim(da)[1]
    k = dim(da)[2]
    p = 0
    q = 0
    if (length(locAR) > 0) 
        p = floor(dim(locAR)[2]/k)
    if (length(locMA) > 0) 
        q = floor(dim(locMA)[2]/k)
    cat("iniSCM p and q: ", c(p, q), "\n")
    pq = p + q
    if (pq < 1) {
        cat("The series is white noise. No estimation is needed", 
            "\n")
    }
    est = NULL
    estse = NULL
    ist = pq + 1
    for (i in 1:k) {
        X = NULL
        loki = Tdx[i]
        Y = da[ist:nT, loki]
        if (inc.mean) 
            X = rep(1, (nT - pq))
        for (j in 1:k) {
            if (j != loki) {
                if (locTmtx[i, j] > 1) 
                  X = cbind(X, -da[ist:nT, j])
            }
        }
        pi = scms[i, 1]
        qi = scms[i, 2]
        if (pi > 0) {
            for (lag in 1:pi) {
                jst = (lag - 1) * k
                for (j in 1:k) {
                  if (locAR[i, jst + j] > 1) {
                    tmp = da[(ist - lag):(nT - lag), j]
                    X = cbind(X, tmp)
                  }
                }
            }
        }
        if (qi > 0) {
            for (lag in 1:qi) {
                jst = (lag - 1) * k
                for (j in 1:k) {
                  if (locMA[i, jst + j] > 1) {
                    tmp = at[(ist - lag):(nT - lag), j]
                    X = cbind(X, tmp)
                  }
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
    iniSCM <- list(par = est, se = estse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
