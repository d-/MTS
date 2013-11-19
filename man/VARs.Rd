\name{VARs}
\alias{VARs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR Model of Selected Lags
%%  ~~function to do ... ~~
}
\description{This is a modified version of VAR command by allowing the 
users to specify which AR lags to be included in the model. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARs(x, lags, include.mean = T, output = T, fixed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of k-dimensional time series with T observations
%%     ~~Describe \code{x} here~~
}
  \item{lags}{A vector of non-zero AR lags. For instance, lags=c(1,3) 
denotes a VAR(3) model with Phi2 = 0.
%%     ~~Describe \code{lags} here~~
}
  \item{include.mean}{A logical switch to include the mean vector
%%     ~~Describe \code{include.mean} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
  \item{fixed}{A logical matrix to fix parameters to zero.
%%     ~~Describe \code{fixed} here~~
}
}
\details{Performs VAR estimation by allowing certain lag 
coefficient matrices being zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{VAR estimation
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ruey S. Tsay (Chapter 2)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{VAR command
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("q-gdp-ukcaus.txt",header=T)
zt=log(da[,3:5])
m1=VARs(zt,lags=c(1,2,4))
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, lags, include.mean = T, output = T, fixed = NULL) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nT = nrow(x)
    k = ncol(x)
    nlags = length(lags)
    if (nlags < 0) {
        lags = c(1)
        nlags = 1
    }
    lags = sort(lags)
    idm = k * nlags + 1
    p = lags[nlags]
    if (p < 1) 
        p = 1
    ne = nT - p
    ist = p + 1
    y = x[ist:nT, ]
    jj = lags[1]
    if (include.mean) {
        xmtx = cbind(rep(1, ne), x[(ist - jj):(nT - jj), ])
    }
    else {
        xmtx = cbind(x[(ist - jj):(T - jj), ])
    }
    if (nlags > 1) {
        for (i in 2:nlags) {
            jj = lags[i]
            xmtx = cbind(xmtx, x[(ist - jj):(nT - jj), ])
        }
    }
    xmtx = as.matrix(xmtx)
    ndim = dim(xmtx)[2]
    if (length(fixed) == 0) {
        paridx = matrix(1, ndim, k)
    }
    else {
        paridx = fixed
    }
    res = NULL
    beta = matrix(0, ndim, k)
    sdbeta = matrix(0, ndim, k)
    for (i in 1:k) {
        idx = c(1:ndim)[paridx[, i] == 1]
        resi = y[, i]
        if (length(idx) > 0) {
            xm = as.matrix(xmtx[, idx])
            xpx = t(xm) \%*\% xm
            xpxinv = solve(xpx)
            xpy = t(xm) \%*\% as.matrix(y[, i], ne, 1)
            betai = xpxinv \%*\% xpy
            beta[idx, i] = betai
            resi = y[, i] - xm \%*\% betai
            sse = sum(resi * resi)/nT
            dd = diag(xpxinv)
            sdbeta[idx, i] = sqrt(sse * dd)
        }
        res = cbind(res, resi)
    }
    sse = t(res) \%*\% res/nT
    Phi = NULL
    Ph0 = rep(0, k)
    if (output) {
        jst = 0
        if (include.mean) {
            Ph0 = beta[1, ]
            se = sdbeta[1, ]
            cat("Constant term:", "\n")
            print(Ph0, digits = 4)
            cat("Std error:", "\n")
            print(se, digits = 4)
            jst = 1
        }
        cat("AR coefficient matrix:", "\n")
        Phi = matrix(0, k, p * k)
        for (i in 1:nlags) {
            ord = lags[i]
            cat("AR(", ord, ")-matrix:", "\n")
            phi = t(beta[(jst + 1):(jst + k), ])
            se = t(sdbeta[(jst + 1):(jst + k), ])
            print(phi, digits = 3)
            cat("Standard error:", "\n")
            print(se, digits = 3)
            jst = jst + k
            cat("      ", "\n")
            kdx = (ord - 1) * k
            Phi[, (kdx + 1):(kdx + k)] = phi
        }
        cat("Residuals cov-mtx:", "\n")
        print(sse)
        cat("       ", "\n")
        dd = det(sse)
        cat("det(SSE) = ", dd, "\n")
        d1 = log(dd)
        aic = d1 + (2 * nlags * k * k)/nT
        bic = d1 + log(nT) * nlags * k * k/nT
        cat("AIC = ", aic, "\n")
        cat("BIC = ", bic, "\n")
    }
    VARs <- list(data = x, lags = lags, order = p, cnst = include.mean, 
        coef = beta, aic = aic, bic = bic, residuals = res, secoef = sdbeta, 
        Sigma = sse, Phi = Phi, Ph0 = Ph0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
