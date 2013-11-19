\name{REGts}
\alias{REGts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regression Model with Time Series Errors
%%  ~~function to do ... ~~
}
\description{Perform maximum likelihood estimation of a multivariate 
linear regression model with time-series errors
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
REGts(zt, p, xt, include.mean = T, fixed = NULL, par = NULL, se.par = NULL, details = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{p}{The VAR order
%%     ~~Describe \code{p} here~~
}
  \item{xt}{A T-by-v data matrix of independent variables, where 
v denotes the number of independent variables (excluding constant 1).
%%     ~~Describe \code{xt} here~~
}
  \item{include.mean}{A logical switch to include the constant term
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used to set parameters to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{par}{Initial parameter estimates of the beta coefficients
%%     ~~Describe \code{par} here~~
}
  \item{se.par}{Standard errors of the parameters in par.
%%     ~~Describe \code{se.par} here~~
}
  \item{details}{A logical switch to control the output
%%     ~~Describe \code{details} here~~
}
}
\details{Perform maximum likelihood estimation of a multivariate linear 
regression model with time series errors. Use multivariate linear regression 
to obtain initial estimates of regression coefficients if not provided
%%  ~~ If necessary, more details than the description above ~~
}
\value{Estimation results of a regression model with time series errors.
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
m1=REGts(zt,2,xt)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (zt, p, xt, include.mean = T, fixed = NULL, par = NULL, 
    se.par = NULL, details = F) 
{
    if (!is.matrix(zt)) 
        zt = as.matrix(zt)
    if (!is.matrix(xt)) 
        xt = as.matrix(xt)
    nT = dim(zt)[1]
    k = dim(zt)[2]
    kx = dim(xt)[2]
    if (p < 0) 
        p = 0
    if (length(par) < 1) {
        m1 = Mlm(zt, xt, constant = include.mean, output = F)
        par = c(m1$beta)
        r1 = dim(m1$beta)[1]
        se.par = c(m1$se.beta)
        res = m1$residuals
        if (p > 0) {
            m2 = VAR(res, p, include.mean = F, output = F)
            par1 = c(m2$coef)
            par = c(par, par1)
            se1 = c(m2$secoef)
            se.par = c(se.par, se1)
            r2 = dim(m2$coef)[1]
        }
        if (length(fixed) < 1) {
            fixed = matrix(1, (r1 + r2), k)
        }
    }
    r1 = kx
    if (include.mean) 
        r1 = r1 + 1
    r2 = p * k
    RegZt <<- zt
    RegXt <<- xt
    RegCnt <<- include.mean
    RegP <<- p
    Regfixed <<- fixed
    cat("Number of parameters: ", length(par), "\n")
    cat("initial estimates: ", par, "\n")
    lowerBounds = par
    upperBounds = par
    npar = length(par)
    mult = 1.5
    if (npar > 10) 
        mult = 1.3
    for (j in 1:npar) {
        lowerBounds[j] = par[j] - mult * se.par[j]
        upperBounds[j] = par[j] + mult * se.par[j]
    }
    cat("Par. Lower-bounds: ", lowerBounds, "\n")
    cat("Par. Upper-bounds: ", upperBounds, "\n")
    if (details) {
        fit = nlminb(start = par, objective = lRegts, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = par, objective = lRegts, control = list(step.min = 0.2, 
            step.max = 0.5), lower = lowerBounds, upper = upperBounds)
    }
    epsilon = 1e-04 * fit$par
    Hessian = matrix(0, ncol = npar, nrow = npar)
    for (i in 1:npar) {
        for (j in 1:npar) {
            x1 = x2 = x3 = x4 = fit$par
            x1[i] = x1[i] + epsilon[i]
            x1[j] = x1[j] + epsilon[j]
            x2[i] = x2[i] + epsilon[i]
            x2[j] = x2[j] - epsilon[j]
            x3[i] = x3[i] - epsilon[i]
            x3[j] = x3[j] + epsilon[j]
            x4[i] = x4[i] - epsilon[i]
            x4[j] = x4[j] - epsilon[j]
            Hessian[i, j] = (lRegts(x1) - lRegts(x2) - lRegts(x3) + 
                lRegts(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    est = fit$par
    cat("Final   Estimates: ", est, "\n")
    se.coef = sqrt(diag(solve(Hessian)))
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", 
        " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
    beta = NULL
    se.beta = NULL
    Phi = NULL
    se.Phi = NULL
    icnt = 0
    if (r1 > 0) {
        beta = matrix(0, r1, k)
        se.beta = matrix(1, r1, k)
        for (i in 1:k) {
            idx = c(1:r1)[fixed[1:r1, i] > 0]
            ii = length(idx)
            if (ii > 0) {
                beta[idx, i] = est[(icnt + 1):(icnt + ii)]
                se.beta[idx, i] = se.coef[(icnt + 1):(icnt + 
                  ii)]
                icnt = icnt + ii
            }
        }
        cat("======= ", "\n")
        cat("Coefficient matrix for constant + exogenous variable", 
            "\n")
        cat("Estimates: ", "\n")
        print(round(t(beta), 3))
        cat("Standard errors: ", "\n")
        print(round(t(se.beta), 3))
    }
    res = zt - RegXmtx \%*\% beta
    if (r2 > 0) {
        Phi = matrix(0, r2, k)
        se.Phi = matrix(1, r2, k)
        fix1 = fixed[(r1 + 1):(r1 + r2), ]
        for (i in 1:k) {
            idx = c(1:r2)[fix1[1:r2, i] > 0]
            ii = length(idx)
            if (ii > 0) {
                Phi[idx, i] = est[(icnt + 1):(icnt + ii)]
                se.Phi[idx, i] = se.coef[(icnt + 1):(icnt + ii)]
                icnt = icnt + ii
            }
        }
        cat("VAR coefficient matrices: ", "\n")
        for (i in 1:p) {
            kdx = (i - 1) * k
            cat("AR(", i, ") coefficient: ", "\n")
            phi = t(Phi[(kdx + 1):(kdx + k), ])
            print(round(phi, 3))
            cat("standard errors:", "\n")
            sephi = t(se.Phi[(kdx + 1):(kdx + k), ])
            print(round(sephi, 3))
        }
    }
    ist = p + 1
    resi = res
    if (p > 0) {
        tsxm = res[(ist - 1):(nT - 1), ]
        if (p > 1) {
            for (j in 2:p) {
                tsxm = cbind(tsxm, res[(ist - j):(nT - j), ])
            }
        }
        resi = res[ist:nT, ] - as.matrix(tsxm) \%*\% Phi
    }
    sig = t(resi) \%*\% resi/(nT - p)
    cat("Residual Covariance matrix: ", "\n")
    print(sig, digits = 4)
    d1 = log(det(sig))
    aic = d1 + 2 * npar/(nT - p)
    bic = d1 + log(nT - p) * npar/(nT - p)
    cat("============", "\n")
    cat("Information criteria: ", "\n")
    cat("AIC: ", aic, "\n")
    cat("BIC: ", bic, "\n")
    coef = rbind(beta, Phi)
    se.coef = rbind(se.beta, se.Phi)
    REGts <- list(data = zt, xt = xt, aror = p, include.mean = include.mean, 
        Phi = t(Phi), se.Phi = t(se.Phi), beta = t(beta), se.beta = t(se.beta), 
        residuals = resi, Sigma = sig, coef = coef, se.coef = se.coef)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
