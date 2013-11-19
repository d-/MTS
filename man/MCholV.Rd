\name{MCholV}
\alias{MCholV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Cholesky Volatility Model
%%  ~~function to do ... ~~
}
\description{Use Cholesky decomposition to obtain multivariate volatility models
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MCholV(rtn, size = 36, lambda = 0.96, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns.
%%     ~~Describe \code{rtn} here~~
}
  \item{size}{The initial sample size used to start recursive least squares 
estimation
%%     ~~Describe \code{size} here~~
}
  \item{lambda}{The exponential smoothing parameter. Default is 0.96.
%%     ~~Describe \code{lambda} here~~
}
  \item{p}{VAR order for the mean equation. Default is 0.
%%     ~~Describe \code{p} here~~
}
}
\details{Use recursive least squares to perform the time-varying 
Cholesky decomposition. The least squares estimates are then smoothed 
via the exponentially weighted moving-average method with decaying rate 0.96.
University GARCH(1,1) model is used for the innovations of each 
linear regression. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{Sigma.t: Volatility matrices, betat: Time-varying beta's.
bt = innovation series, and Vol: volatility series of innovations.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 7)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{fGarch
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-ibmspko-6111.txt",header=T)
rtn=log(da[,2:4]+1)
require(fGarch)
m3=MCholV(rtn)
names(m3)
Sample means: 0.007728
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rtn, size = 36, lambda = 0.96, p = 0) 
{
    if (!is.matrix(rtn)) 
        rtn = as.matrix(rtn)
    if (size <= 0) 
        size = 36
    if (p < 0) 
        p = 0
    mu = apply(rtn, 2, mean)
    if (p == 0) {
        RTN = scale(rtn, center = T, scale = F)
        cat("Sample means: ", round(mu, 6), "\n")
    }
    else {
        m1 = VAR(rtn, p)
        RTN = m1$residuals
    }
    nT = dim(RTN)[1]
    k = dim(RTN)[2]
    Mtxcoef = matrix(0, k, 3)
    Mtxse = matrix(0, k, 3)
    Mtxtval = matrix(0, k, 3)
    bt = RTN[(size + 1):nT, 1]
    m1 = garchFit(~garch(1, 1), data = bt, include.mean = F, 
        trace = F)
    Mtxcoef[1, ] = m1@fit$par
    Mtxse[1, ] = c(m1@fit$se.coef)
    Mtxtval[1, ] = m1@fit$tval
    cat("Estimation of the first component", "\n")
    cat("Estimate (alpha0, alpha1, beta1): ", round(Mtxcoef[1, 
        ], 6), "\n")
    cat("s.e.                            : ", round(Mtxse[1, 
        ], 6), "\n")
    cat("t-value                         : ", round(Mtxtval[1, 
        ], 6), "\n")
    VOL = volatility(m1)^2
    BetaU = NULL
    for (i in 2:k) {
        betat = NULL
        y = RTN[, i]
        x = RTN[, 1:(i - 1)]
        m2 = RLS(y, x, ist = size)
        betat = cbind(betat, -m2$beta)
        bt = cbind(bt, m2$resi)
        betat = as.matrix(betat)
        BetaM = apply(betat, 2, mean)
        betat = scale(betat, center = T, scale = F)
        RTN1 <- RTN[(size + 1):nT, 1:i]
        k1 = dim(betat)[2]
        nT1 = dim(betat)[1]
        Beta = NULL
        for (j in 1:k1) {
            b1 = c(betat[1, j], (1 - lambda) * betat[-nT1, j])
            bb = filter(b1, lambda, "r", init = 0) + BetaM[j]
            Beta = cbind(Beta, bb)
        }
        BetaU = cbind(BetaU, Beta)
        wk = cbind(Beta, rep(1, nT1))
        bit = apply(RTN1 * wk, 1, sum)
        bt[, i] = bit
        m3 = garchFit(~garch(1, 1), data = bit, include.mean = F, 
            trace = F)
        cat("Component ", i, " Estimation Results (residual series):", 
            "\n")
        Mtxcoef[i, ] = m3@fit$par
        Mtxse[i, ] = c(m3@fit$se.coef)
        Mtxtval[i, ] = m3@fit$tval
        cat("Estimate (alpha0, alpha1, beta1): ", round(Mtxcoef[i, 
            ], 6), "\n")
        cat("s.e.                            : ", round(Mtxse[i, 
            ], 6), "\n")
        cat("t-value                         : ", round(Mtxtval[i, 
            ], 6), "\n")
        VOL = cbind(VOL, volatility(m3)^2)
    }
    BetaU = as.matrix(BetaU)
    Sigma.t = NULL
    for (it in 1:nT1) {
        A = diag(k)
        icnt = 0
        for (j in 2:k) {
            A[j, 1:(j - 1)] = BetaU[it, (icnt + 1):(icnt + j - 
                1)]
            icnt = icnt + j - 1
        }
        B = LinV(A)
        for (i in 1:k) {
            A[, i] = VOL[it, i] * B[, i]
        }
        Sigma.t = rbind(Sigma.t, c(A \%*\% t(B)))
    }
    MCholV <- list(betat = -BetaU, bt = bt, Vol = VOL, Sigma.t = Sigma.t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
