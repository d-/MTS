\name{VARX}
\alias{VARX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR Model with Exogenous Variables
%%  ~~function to do ... ~~
}
\description{Estimation of a VARX model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARX(zt, p, xt = NULL, m = 0, include.mean = T, fixed = NULL, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{p}{The VAR order
%%     ~~Describe \code{p} here~~
}
  \item{xt}{A T-by-m data matrix of m exogenous variables
%%     ~~Describe \code{xt} here~~
}
  \item{m}{The number of exogenous variables
%%     ~~Describe \code{m} here~~
}
  \item{include.mean}{A logical switch to include the constant vector. 
Default is to include the constant.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix for setting parameters to zero.
%%     ~~Describe \code{fixed} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{Performs least squares estimation of a VARX(p,s) model
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a VARX model, including model and data
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 6)
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
da=read.table("m-gasoil.txt",header=T)
zt=da[,3:4]; xt=da[,5:6]
m1=VARX(zt,2,xt,1)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (zt, p, xt = NULL, m = 0, include.mean = T, fixed = NULL, 
    output = T) 
{
    zt = as.matrix(zt)
    if (length(xt) < 1) {
        m = -1
        kx = 0
    }
    else {
        xt = as.matrix(xt)
        kx = dim(xt)[2]
    }
    if (p < 0) 
        p = 0
    ist = max(p, m) + 1
    nT = dim(zt)[1]
    k = dim(zt)[2]
    yt = zt[ist:nT, ]
    xmtx = NULL
    if (include.mean) 
        xmtx = rep(1, (nT - ist + 1))
    if (p > 0) {
        for (i in 1:p) {
            xmtx = cbind(xmtx, zt[(ist - i):(nT - i), ])
        }
    }
    if (m > -1) {
        for (j in 0:m) {
            xmtx = cbind(xmtx, xt[(ist - j):(nT - j), ])
        }
    }
    p1 = dim(xmtx)[2]
    nobe = dim(xmtx)[1]
    if (length(fixed) < 1) {
        xpx = t(xmtx) \%*\% xmtx
        xpy = t(xmtx) \%*\% yt
        xpxi = solve(xpx)
        beta = xpxi \%*\% xpy
        resi = as.matrix(yt - xmtx \%*\% beta)
        sig = crossprod(resi, resi)/nobe
        co = kronecker(sig, xpxi)
        se = sqrt(diag(co))
        se.beta = matrix(se, nrow(beta), k)
        npar = nrow(beta) * k
        d1 = log(det(sig))
        aic = d1 + 2 * npar/nobe
        bic = d1 + (npar * log(nobe))/nobe
    }
    else {
        beta = matrix(0, p1, k)
        se.beta = matrix(1, p1, k)
        resi = yt
        npar = 0
        for (i in 1:k) {
            idx = c(1:p1)[fixed[, i] > 0]
            npar = npar + length(idx)
            if (length(idx) > 0) {
                xm = as.matrix(xmtx[, idx])
                y1 = matrix(yt[, i], nobe, 1)
                xpx = t(xm) \%*\% xm
                xpy = t(xm) \%*\% y1
                xpxi = solve(xpx)
                beta1 = xpxi \%*\% xpy
                res = y1 - xm \%*\% beta1
                sig1 = sum(res^2)/nobe
                se = sqrt(diag(xpxi) * sig1)
                beta[idx, i] = beta1
                se.beta[idx, i] = se
                resi[, i] = res
            }
        }
        sig = crossprod(resi, resi)/nobe
        d1 = log(det(sig))
        aic = d1 + 2 * npar/nobe
        bic = d1 + log(nobe) * npar/nobe
    }
    Ph0 = NULL
    icnt = 0
    if (include.mean) {
        Ph0 = beta[1, ]
        icnt = icnt + 1
        cat("constant term: ", "\n")
        cat("est: ", round(Ph0, 4), "\n")
        cat(" se: ", round(se.beta[1, ], 4), "\n")
    }
    Phi = NULL
    if (p > 0) {
        Phi = t(beta[(icnt + 1):(icnt + k * p), ])
        sePhi = t(se.beta[(icnt + 1):(icnt + k * p), ])
        for (j in 1:p) {
            cat("AR(", j, ") matrix", "\n")
            jcnt = (j - 1) * k
            print(round(Phi[, (jcnt + 1):(jcnt + k)], 3))
            cat("standard errors", "\n")
            print(round(sePhi[, (jcnt + 1):(jcnt + k)], 3))
        }
        icnt = icnt + k * p
    }
    if (m > -1) {
        cat("Coefficients of exogenous", "\n")
        Beta = t(beta[(icnt + 1):(icnt + (m + 1) * kx), ])
        seBeta = t(se.beta[(icnt + 1):(icnt + (m + 1) * kx), 
            ])
        for (i in 0:m) {
            jdx = i * kx
            cat("lag-", i, " coefficient matrix", "\n")
            print(round(Beta[, (jdx + 1):(jdx + kx)], 3))
            cat("standard errors", "\n")
            print(round(seBeta[, (jdx + 1):(jdx + kx)], 3))
        }
    }
    cat("Residual Covariance Matrix", "\n")
    print(round(sig, 5))
    cat("===========", "\n")
    cat("Information criteria: ", "\n")
    cat("AIC: ", aic, "\n")
    cat("BIC: ", bic, "\n")
    VARX <- list(data = zt, xt = xt, aror = p, m = m, Ph0 = Ph0, 
        Phi = Phi, beta = Beta, residuals = resi, Sigma = sig, 
        coef = beta, se.coef = se.beta, include.mean = include.mean)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
