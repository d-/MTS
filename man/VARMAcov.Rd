\name{VARMAcov}
\alias{VARMAcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Autocovariance Matrices of a VARMA Model
%%  ~~function to do ... ~~
}
\description{Uses psi-weights to compute the autocovariance matrices of 
a VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMAcov(Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, trun = 120)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{A k-by-kp matrix consisting of VAR coefficient matrices, 
Phi = [Phi1, Phi2, ..., Phip].
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{A k-by-kq matrix consisting of VMA coefficient magrices,
Theta = [Theta1, Theta2, ..., Thetaq]
%%     ~~Describe \code{Theta} here~~
}
  \item{Sigma}{Covariance matrix of the innovations (k-by-k).
%%     ~~Describe \code{Sigma} here~~
}
  \item{lag}{Number of cross-covariance matrices to be computed. 
Default is 12.
%%     ~~Describe \code{lag} here~~
}
  \item{trun}{The lags of pis-weights used in calculation. Default is 120.
%%     ~~Describe \code{trun} here~~
}
}
\details{Use psi-weight matrices to compute approximate autocovariance 
matrices of a VARMA model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The theoretical auto-covariance matrices of a given VARMA model.
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
Phi=matrix(c(0.2,-0.6,0.3,1.1),2,2)
Sig=matrix(c(4,1,1,1),2,2)
VARMAcov(Phi=Phi,Sigma=Sig)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, trun = 120) 
{
    m1 = PSIwgt(Phi = Phi, Theta = Theta, lag = trun, plot = FALSE)
    Psi = m1$psi.weight
    nc = dim(Psi)[2]
    k = dim(Psi)[1]
    if (is.null(Sigma)) {
        wk = PSI
    }
    else {
        wk = NULL
        for (i in 0:trun) {
            ist = i * k
            wk = cbind(wk, Psi[, (ist + 1):(ist + k)] \%*\% Sigma)
        }
    }
    Gam0 = wk \%*\% t(Psi)
    SE = diag(1/sqrt(diag(Gam0)))
    covmtx = Gam0
    cormtx = SE \%*\% Gam0 \%*\% SE
    for (i in 1:lag) {
        ist = i * k
        Gami = wk[, (ist + 1):nc] \%*\% t(Psi[, 1:(nc - ist)])
        covmtx = cbind(covmtx, Gami)
        cormtx = cbind(cormtx, SE \%*\% Gami \%*\% SE)
    }
    for (i in 0:lag) {
        ist = i * k
        cat("Auto-Covariance matrix of lag: ", i, "\n")
        print(round(covmtx[, (ist + 1):(ist + k)], 5))
    }
    for (i in 0:lag) {
        ist = i * k
        cat("cross correlation matrix of lag: ", i, "\n")
        print(round(cormtx[, (ist + 1):(ist + k)], 4))
    }
    VARMAcov <- list(autocov = covmtx, ccm = cormtx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
