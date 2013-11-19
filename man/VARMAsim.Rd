\name{VARMAsim}
\alias{VARMAsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generating a VARMA Process
%%  ~~function to do ... ~~
}
\description{Performs simulation of a given VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMAsim(nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL, theta = NULL, skip = 200, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nobs}{Sample size
%%     ~~Describe \code{nobs} here~~
}
  \item{arlags}{The exact lags of the VAR matrix polynomial.
%%     ~~Describe \code{arlags} here~~
}
  \item{malags}{The exact lags of the VMA matrix polynomial.
%%     ~~Describe \code{malags} here~~
}
  \item{cnst}{Constant vector, Phi0
%%     ~~Describe \code{cnst} here~~
}
  \item{phi}{Matrix of VAR coefficient matrices in the order of the given arlags.
%%     ~~Describe \code{phi} here~~
}
  \item{theta}{Matrix of VMA coefficient matrices in the order of the given 
malags.
%%     ~~Describe \code{theta} here~~
}
  \item{skip}{The number of initial data to be omitted. Deafult is 200.
%%     ~~Describe \code{skip} here~~
}
  \item{sigma}{Covariance matrix (k-by-k, positive definite) 
of the innovations
%%     ~~Describe \code{sigma} here~~
}
}
\details{Use multivariate Gaussian distribution to generate random shocks. 
Then, generate a given VARMA model
%%  ~~ If necessary, more details than the description above ~~
}
\value{The vector time series and the innovations used
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 3)
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
p1=matrix(c(0.2,-0.6,0.3,1,1),2,2)
sig=matrix(c(4,0.8,0.8,1),2,2)
th1=matrix(c(-0.5,0,0,-0.6),2,2)
m1=VARMAsim(300,arlags=c(1),malags=c(1),phi=p1,theta=th1,sigma=sig)
zt=m1$series
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL, 
    theta = NULL, skip = 200, sigma) 
{
    if (!is.matrix(sigma)) 
        sigma = as.matrix(sigma)
    k = nrow(sigma)
    nT = nobs + skip
    at = rmvnorm(nT, rep(0, k), sigma)
    nar = length(arlags)
    p = 0
    if (nar > 0) {
        arlags = sort(arlags)
        p = arlags[nar]
    }
    q = 0
    nma = length(malags)
    if (nma > 0) {
        malags = sort(malags)
        q = malags[nma]
    }
    ist = max(p, q) + 1
    zt = matrix(0, nT, k)
    if (length(cnst) == 0) 
        cnst = rep(0, k)
    for (it in ist:nT) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
            for (j in 1:nma) {
                jdx = (j - 1) * k
                thej = theta[, (jdx + 1):(jdx + k)]
                atm = matrix(at[it - malags[j], ], 1, k)
                tmp = tmp - atm \%*\% t(thej)
            }
        }
        if (nar > 0) {
            for (i in 1:nar) {
                idx = (i - 1) * k
                phj = phi[, (idx + 1):(idx + k)]
                ztm = matrix(zt[it - arlags[i], ], 1, k)
                tmp = tmp + ztm \%*\% t(phj)
            }
        }
        zt[it, ] = cnst + tmp
    }
    zt = zt[(1 + skip):nT, ]
    at = at[(1 + skip):nT, ]
    VARMAsim <- list(series = zt, noises = at)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
