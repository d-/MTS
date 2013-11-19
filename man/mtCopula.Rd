\name{mtCopula}
\alias{mtCopula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mulivariate t-Copula for Volatility
%%  ~~function to do ... ~~
}
\description{Fits a t-copula to k-dimensional standardized return series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mtCopula(rt, g1, g2, grp = NULL, th0 = NULL, m = 0, include.th0 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{A T-by-k data matrix of k standardized time series (after 
univariate volatility modeling)
%%     ~~Describe \code{rt} here~~
}
  \item{g1}{lamda1 parameter, nononegative and less than 1
%%     ~~Describe \code{g1} here~~
}
  \item{g2}{lambda2 parameter, nonnegative and satisfying lambda1+lambda2 < 1.
%%     ~~Describe \code{g2} here~~
}
  \item{grp}{a vector to indicate the number of assets divided into groups.
Default means each individual asset forms a group.
%%     ~~Describe \code{grp} here~~
}
  \item{th0}{initial estimate of theta0
%%     ~~Describe \code{th0} here~~
}
  \item{m}{number of lags used to estimate the local theta-angles
%%     ~~Describe \code{m} here~~
}
  \item{include.th0}{A logical switch to include theta0 in estimation. 
Default is to inlcude.
%%     ~~Describe \code{include.th0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{estimates: parameter estimates. Hessian: Hessian matrix.
rho.t: correlation matrices. theta.t: angel mtrices
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-ibmspko-6111.txt",header=T)
rtn=log(da[,-1]+1)
m1=dccPre(rtn,cond.dist="std")
Vol=m1$marVol; eta=m1$sresi
m2=mtCopula(eta,0.8,0.04)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rt, g1, g2, grp = NULL, th0 = NULL, m = 0, include.th0 = TRUE) 
{
    if (!is.matrix(rt)) 
        rt = as.matrix(rt)
    k = dim(rt)[2]
    nT = dim(rt)[1]
    mtCopuladata <<- rt
    if (is.null(grp)) 
        grp = rep(1, k)
    if (length(th0) == 0) {
        mi = SCCor(rt, nT, nT, grp)$conCor
        THmi = CorTheta(mi)$THmtx
        th0 = TH2idth(THmi, grp)
    }
    if (!include.th0) {
        cat("Value of angles: ", "\n")
        print(th0)
    }
    g = length(grp)
    ncor = g * (g - 1)/2
    for (i in 1:g) {
        if (grp[i] > 1) 
            ncor = ncor + 1
    }
    mtCopulancor <<- ncor
    if (m <= 0) 
        m = ncor + 1
    ist = m + 1
    mgsize = max(grp)
    Theta = matrix(1, m, 1) \%*\% matrix(th0, 1, length(th0))
    for (t in ist:nT) {
        V1 = SCCor(rt, t - 1, m, grp)$conCor
        Tmp = CorTheta(V1)$THmtx
        if (mgsize == 1) {
            thet = Vech(t(Tmp))
        }
        else {
            thet = TH2idth(Tmp, grp)
        }
        Theta = rbind(Theta, thet)
    }
    mtCopulaTheta <<- Theta
    mtCopulath0 <<- th0
    mtCopulagrp <<- grp
    mtCopulam <<- m
    mtCopulaincludeth0 <<- include.th0
    if (include.th0) {
        par = c(7, g1, g2, th0)
        c1 = c(5.1, 0.2, 1e-04, th0 * 0.8)
        c2 = c(20, 0.95, 0.04999999, th0 * 1.1)
    }
    else {
        par = c(7, g1, g2)
        c1 = c(5.1, 0.2, 1e-05)
        c2 = c(20, 0.95, 0.0499999)
    }
    cat("Lower limits: ", c1, "\n")
    cat("Upper limits: ", c2, "\n")
    m1 = optim(par, mtlikeC, method = "L-BFGS-B", lower = c1, 
        upper = c2, hessian = T)
    est = m1$par
    H = m1$hessian
    Hi = solve(H)
    se.coef = sqrt(diag(Hi))
    cat("estimates:  ", est, "\n")
    cat("std.errors: ", se.coef, "\n")
    cat("t-values:   ", est/se.coef, "\n")
    est = m1$par
    npar = length(par)
    epsilon = 1e-04 * est
    epsilon[1] = 0.3 * est[1]
    Hessian = matrix(0, ncol = npar, nrow = npar)
    for (i in 1:npar) {
        for (j in 1:npar) {
            x1 = x2 = x3 = x4 = m1$par
            x1[i] = x1[i] + epsilon[i]
            x1[j] = x1[j] + epsilon[j]
            x2[i] = x2[i] + epsilon[i]
            x2[j] = x2[j] - epsilon[j]
            x3[i] = x3[i] - epsilon[i]
            x3[j] = x3[j] + epsilon[j]
            x4[i] = x4[i] - epsilon[i]
            x4[j] = x4[j] - epsilon[j]
            Hessian[i, j] = (mtlikeC(x1) - mtlikeC(x2) - mtlikeC(x3) + 
                mtlikeC(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    d1 = det(Hessian)
    if (d1 < 1e-13) {
        se.coef = rep(1, npar)
    }
    else {
        se.coef = sqrt(diag(solve(Hessian)))
    }
    cat("Alternative numerical estimates of se:", "\n")
    cat("st.errors: ", se.coef, "\n")
    cat("t-values:  ", est/se.coef, "\n")
    mf = mtCopulaVol(est)
    rho.t = mf$rho.t
    angles = mf$angles
    mtCopula <- list(estimates = est, Hessian = H, rho.t = rho.t, 
        theta.t = angles)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
