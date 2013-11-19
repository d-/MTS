\name{FEVdec}
\alias{FEVdec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forecast Error Variance Decomposition
%%  ~~function to do ... ~~
}
\description{Computes the forecast error variance decomposition of a 
VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FEVdec(Phi, Theta, Sig, lag = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{VAR coefficient matrices in the form Phi=[Phi1, Phi2, ..., Phip], 
a k-by-kp matrix.
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{VMA coefficient matrices in form form Theta=[Theta1, Theta2, 
..., Thetaq], a k-by-kq matrix.
%%     ~~Describe \code{Theta} here~~
}
  \item{Sig}{The residual covariance matrix Sigma, a k-by-k positive definite 
matrix.
%%     ~~Describe \code{Sig} here~~
}
  \item{lag}{The number of lags of forecast errors variance to be computed. 
Default is 4.
%%     ~~Describe \code{lag} here~~
}
}
\details{Use the psi-weight matrices to compute the forecast error covariance 
and use Cholesky decomposition to perform the decomposition
%%  ~~ If necessary, more details than the description above ~~
}
\value{irf: impulse response matrices, orthirf: orthogonal impulse response matrices, Omega: forecast error variance matrices, and OmegaR: forecast error 
variance decomposition
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
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
theta1=matrix(c(-0.5,0,0,-0.6),2,2)
Sig=matrix(c(3,1,1,1),2,2)
m1=FEVdec(p1,theta1,Sig)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Phi, Theta, Sig, lag = 4) 
{
    if (length(Phi) > 0) {
        if (!is.matrix(Phi)) 
            Phi = as.matrix(Phi)
    }
    if (length(Theta) > 0) {
        if (!is.matrix(Theta)) 
            Theta = as.matrix(Theta)
    }
    if (!is.matrix(Sig)) 
        Sig = as.matrix(Sig)
    if (lag < 1) 
        lag = 1
    p = 0
    if (length(Phi) > 0) {
        k = nrow(Phi)
        m = ncol(Phi)
        p = floor(m/k)
    }
    q = 0
    if (length(Theta) > 0) {
        k = dim(Theta)[1]
        m = dim(Theta)[2]
        q = floor(m/k)
    }
    cat("Order of the ARMA mdoel: ", "\n")
    print(c(p, q))
    Si = diag(rep(1, k))
    if (q > 0) {
        Si = cbind(Si, -Theta)
    }
    m = (lag + 1) * k
    m1 = (q + 1) * k
    if (m > m1) {
        Si = cbind(Si, matrix(0, k, (m - m1)))
    }
    if (p > 0) {
        for (i in 1:lag) {
            if (i <= p) {
                idx = (i - 1) * k
                tmp = Phi[, (idx + 1):(idx + k)]
            }
            else {
                tmp = matrix(0, k, k)
            }
            jj = i - 1
            jp = min(jj, p)
            if (jp > 0) {
                for (j in 1:jp) {
                  jdx = (j - 1) * k
                  idx = (i - j) * k
                  w1 = Phi[, (jdx + 1):(jdx + k)]
                  w2 = Si[, (idx + 1):(idx + k)]
                  tmp = tmp + w1 \%*\% w2
                }
            }
            kdx = i * k
            Si[, (kdx + 1):(kdx + k)] = tmp
        }
    }
    orSi = NULL
    m1 = chol(Sig)
    P = t(m1)
    orSi = P
    for (i in 1:lag) {
        idx = i * k
        w1 = Si[, (idx + 1):(idx + k)]
        w2 = w1 \%*\% P
        orSi = cbind(orSi, w2)
    }
    orSi2 = orSi^2
    Ome = orSi2[, 1:k]
    wk = Ome
    for (i in 1:lag) {
        idx = i * k
        wk = wk + orSi2[, (idx + 1):(idx + k)]
        Ome = cbind(Ome, wk)
    }
    FeV = NULL
    OmeRa = Ome[, 1:k]
    FeV = cbind(FeV, apply(OmeRa, 1, sum))
    OmeRa = OmeRa/FeV[, 1]
    for (i in 1:lag) {
        idx = i * k
        wk = Ome[, (idx + 1):(idx + k)]
        FeV = cbind(FeV, apply(wk, 1, sum))
        OmeRa = cbind(OmeRa, wk/FeV[, (i + 1)])
    }
    cat("Standard deviation of forecast error: ", "\n")
    print(sqrt(FeV))
    cat("Forecast-Error-Variance Decomposition", "\n")
    for (i in 1:(lag + 1)) {
        idx = (i - 1) * k
        cat("Forecast horizon: ", i, "\n")
        Ratio = OmeRa[, (idx + 1):(idx + k)]
        print(Ratio)
    }
    FEVdec <- list(irf = Si, orthirf = orSi, Omega = Ome, OmegaR = OmeRa)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
