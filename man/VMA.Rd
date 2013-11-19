\name{VMA}
\alias{VMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vector Moving Averge Model
%%  ~~function to do ... ~~
}
\description{Performs VMA estimation using the conditional multivariate 
Gaussian likelihood function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMA(da, q = 1, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix of a k-dimensional VMA process with each column 
containing one time series
%%     ~~Describe \code{da} here~~
}
  \item{q}{The order of VMA model
%%     ~~Describe \code{q} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. 
The default is to include the mean vector in estimation.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used to fix parameter to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch to select parameters to be included in estimation
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control the amount of output
%%     ~~Describe \code{details} here~~
}
  \item{thres}{Threshold for t-ratio used to fix parameeter to zero. Deault is 2.
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 3).
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
\examples{da=read.table("m-dec15678-6111.txt",header=T)
rtn=cbind(da$dec5,da$dec8)
zt=log(rtn+1)*100
m1=VMA(zt,q=1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, q = 1, include.mean = T, fixed = NULL, prelim = F, 
    details = F, thres = 2) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    nT = dim(da)[1]
    k = dim(da)[2]
    if (q < 1) 
        q = 1
    kq = k * q
    VMAdata <<- da
    m1 = VARorder(da, q + 12, output = FALSE)
    porder = m1$aicor
    if (porder < 1) 
        porder = 1
    m2 = VAR(da, porder, output = FALSE)
    y = da[(porder + 1):nT, ]
    x = m2$residuals
    m3 = THini(y, x, q, include.mean)
    beta = m3$estimates
    sebeta = m3$se
    nr = dim(beta)[1]
    if (prelim) {
        fixed = matrix(0, nr, k)
        for (j in 1:k) {
            tt = beta[, j]/sebeta[, j]
            idx = c(1:nr)[abs(tt) >= thres]
            fixed[idx, j] = 1
        }
    }
    if (length(fixed) == 0) {
        fixed = matrix(1, nr, k)
    }
    fix1 <<- fixed
    inc.mean <<- include.mean
    MAq <<- q
    par = NULL
    separ = NULL
    VMAcnt = 0
    ist = 0
    if (include.mean) {
        jdx = c(1:k)[fix1[1, ] == 1]
        VMAcnt = length(jdx)
        if (VMAcnt > 0) {
            par = beta[1, jdx]
            separ = sebeta[1, jdx]
        }
        TH = -beta[2:(kq + 1), ]
        seTH = sebeta[2:(kq + 1), ]
        ist = 1
    }
    else {
        TH = -beta
        seTH = sebeta
    }
    for (j in 1:k) {
        idx = c(1:(nr - ist))[fix1[(ist + 1):nr, j] == 1]
        if (length(idx) > 0) {
            par = c(par, TH[idx, j])
            separ = c(separ, seTH[idx, j])
        }
    }
    ParMA <<- par
    VMAcnt <<- VMAcnt
    cat("Number of parameters: ", length(par), "\n")
    cat("initial estimates: ", round(par, 4), "\n")
    lowerBounds = par
    upperBounds = par
    npar = length(par)
    mult = 2
    if ((npar > 10) || (MAq > 2)) 
        mult = 1.2
    for (j in 1:npar) {
        lowerBounds[j] = par[j] - mult * separ[j]
        upperBounds[j] = par[j] + mult * separ[j]
    }
    cat("Par. Lower-bounds: ", round(lowerBounds, 4), "\n")
    cat("Par. Upper-bounds: ", round(upperBounds, 4), "\n")
    if (details) {
        fit = nlminb(start = ParMA, objective = LLKvma, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = ParMA, objective = LLKvma, lower = lowerBounds, 
            upper = upperBounds)
    }
    epsilon = 1e-04 * fit$par
    npar = length(par)
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
            Hessian[i, j] = (LLKvma(x1) - LLKvma(x2) - LLKvma(x3) + 
                LLKvma(x4))/(4 * epsilon[i] * epsilon[j])
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
    cat("---", "\n")
    cat("Estimates in matrix form:", "\n")
    icnt = 0
    ist = 0
    cnt = NULL
    if (include.mean) {
        ist = 1
        cnt = rep(0, k)
        secnt = rep(1, k)
        jdx = c(1:k)[fix1[1, ] == 1]
        icnt = length(jdx)
        if (icnt > 0) {
            cnt[jdx] = est[1:icnt]
            secnt[jdx] = se.coef[1:icnt]
            cat("Constant term: ", "\n")
            cat("Estimates: ", cnt, "\n")
        }
    }
    cat("MA coefficient matrix", "\n")
    TH = matrix(0, kq, k)
    seTH = matrix(1, kq, k)
    for (j in 1:k) {
        idx = c(1:kq)[fix1[(ist + 1):nr, j] == 1]
        jcnt = length(idx)
        if (jcnt > 0) {
            TH[idx, j] = est[(icnt + 1):(icnt + jcnt)]
            seTH[idx, j] = se.coef[(icnt + 1):(icnt + jcnt)]
            icnt = icnt + jcnt
        }
    }
    icnt = 0
    for (i in 1:q) {
        cat("MA(", i, ")-matrix", "\n")
        theta = t(TH[(icnt + 1):(icnt + k), ])
        print(theta, digits = 3)
        icnt = icnt + k
    }
    zt = VMAdata
    if (include.mean) {
        for (i in 1:k) {
            zt[, i] = zt[, i] - cnt[i]
        }
    }
    at = mFilter(zt, t(TH))
    sig = t(at) \%*\% at/nT
    cat(" ", "\n")
    cat("Residuals cov-matrix:", "\n")
    print(sig)
    dd = det(sig)
    d1 = log(dd)
    aic = d1 + 2 * npar/nT
    bic = d1 + log(nT) * npar/nT
    cat("----", "\n")
    cat("aic= ", aic, "\n")
    cat("bic= ", bic, "\n")
    Theta = t(TH)
    if (include.mean) {
        TH = rbind(cnt, TH)
        seTH = rbind(secnt, seTH)
    }
    VMA <- list(data = da, MAorder = q, cnst = include.mean, 
        coef = TH, secoef = seTH, residuals = at, Sigma = sig, 
        Theta = Theta, mu = cnt, aic = aic, bic = bic)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
