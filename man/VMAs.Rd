\name{VMAs}
\alias{VMAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VMA Model with Selected Lags
%%  ~~function to do ... ~~
}
\description{Performs conditional maximum likelihood estimation 
of a VMA model with selected lags in the model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMAs(da, malags, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k matrix of a k-dimensional time series with T observations
%%     ~~Describe \code{da} here~~
}
  \item{malags}{A vector consisting of non-zero MA lags
%%     ~~Describe \code{malags} here~~
}
  \item{include.mean}{A logical switch to include the mean vector
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to fix coefficients to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch concerning initial estimation
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control output level
%%     ~~Describe \code{details} here~~
}
  \item{thres}{A threshold value for setting coefficient estimates to zero
%%     ~~Describe \code{thres} here~~
}
}
\details{A modified version of VMA model by allowing the users to select 
non-zero MA lags 
%%  ~~ If necessary, more details than the description above ~~
}
\value{VMA estimation
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

\seealso{VMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-ibmsp-6111.txt",header=T)
rtn=log(da[,2:3]+1)
m1=VMAs(rtn,malags=c(1,3))
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, malags, include.mean = T, fixed = NULL, prelim = F, 
    details = F, thres = 2) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    nT = dim(da)[1]
    k = dim(da)[2]
    nlags = length(malags)
    if (nlags < 1) {
        malags = c(1)
        nlags = 1
    }
    MAlag <<- sort(malags)
    kq = k * nlags
    q = MAlag[nlags]
    VMAdata <<- da
    m1 = VARorder(da, q + 10, output = FALSE)
    porder = m1$aicor
    m2 = VAR(da, porder, output = FALSE)
    y = da[(porder + 1):nT, ]
    x = m2$residuals
    m3 = THinis(y, x, MAlag, include.mean)
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
    par = NULL
    separ = NULL
    ist = 0
    if (include.mean) {
        jdx = c(1:k)[fix1[1, ] == 1]
        if (length(jdx) > 0) {
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
    cat("Initial estimates: ", round(par, 4), "\n")
    lowerBounds = par
    upperBounds = par
    for (j in 1:length(par)) {
        lowerBounds[j] = par[j] - 2 * separ[j]
        upperBounds[j] = par[j] + 2 * separ[j]
    }
    cat("Par. lower-bounds: ", round(lowerBounds, 4), "\n")
    cat("Par. upper-bounds: ", round(upperBounds, 4), "\n")
    if (details) {
        fit = nlminb(start = par, objective = LLKvmas, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = par, objective = LLKvmas, lower = lowerBounds, 
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
            Hessian[i, j] = (LLKvmas(x1) - LLKvmas(x2) - LLKvmas(x3) + 
                LLKvmas(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    est = fit$par
    cat("Final    Estimates: ", est, "\n")
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
    cnt = rep(0, k)
    secnt = rep(1, k)
    if (inc.mean) {
        ist = 1
        jdx = c(1:k)[fix1[1, ] == 1]
        icnt = length(jdx)
        if (icnt > 0) {
            cnt[jdx] = est[1:icnt]
            secnt = se.coef[1:icnt]
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
    for (i in 1:nlags) {
        ii = MAlag[i]
        cat("MA(", ii, ")-matrix", "\n")
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
    at = zt[1:MAlag[1], ]
    if (MAlag[1] == 1) 
        at = matrix(at, 1, k)
    if (q >= (MAlag[1] + 1)) {
        for (t in (MAlag[1] + 1):q) {
            Past = NULL
            for (ii in 1:nlags) {
                jj = MAlag[ii]
                if ((t - jj) > 0) {
                  Past = c(Past, at[t - jj, ])
                }
                else {
                  Past = c(Past, rep(0, k))
                }
            }
            tmp = zt[t, ] + matrix(Past, 1, kq) \%*\% TH
            at = rbind(at, tmp)
        }
    }
    for (t in (q + 1):nT) {
        Past = NULL
        for (ii in 1:nlags) {
            jj = MAlag[ii]
            Past = c(Past, at[t - jj, ])
        }
        tmp = zt[t, ] + matrix(Past, 1, kq) \%*\% TH
        at = rbind(at, tmp)
    }
    sig = t(at) \%*\% at/nT
    cat(" ", "\n")
    cat("Residuals cov-matrix:", "\n")
    print(sig)
    dd = det(sig)
    d1 = log(dd)
    aic = d1 + 2 * npar/nT
    bic = d1 + log(nT) * npar/nT
    cat("---", "\n")
    cat("aic = ", aic, "\n")
    cat("bic = ", bic, "\n")
    Theta = t(TH)
    if (include.mean) {
        TH = rbind(cnt, TH)
        seTH = rbind(secnt, seTH)
    }
    VMAs <- list(data = da, MAlags = MAlag, cnst = include.mean, 
        coef = TH, secoef = seTH, residuals = at, aic = aic, 
        bic = bic, Sigma = sig, Theta = Theta, mu = cnt, MAorder = q)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
