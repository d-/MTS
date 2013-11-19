\name{VMAe}
\alias{VMAe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of VMA Model with Exact likelihood
%%  ~~function to do ... ~~
}
\description{Estimation of a VMA(q) model using the exact likelihood method. 
Multivariate Gaussian likelihood function is used.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMAe(da, q = 1, include.mean = T, coef0 = NULL, secoef0 = NULL, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix (T-by-k) for a k-dimensional VMA process
%%     ~~Describe \code{da} here~~
}
  \item{q}{The order of a VMA model
%%     ~~Describe \code{q} here~~
}
  \item{include.mean}{A logicak switch to include the mean vector in estimation. Default is to include the mean vector.
%%     ~~Describe \code{include.mean} here~~
}
  \item{coef0}{
%%     ~~Describe \code{coef0} here~~
}
  \item{secoef0}{
%%     ~~Describe \code{secoef0} here~~
}
  \item{fixed}{
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{
%%     ~~Describe \code{details} here~~
}
  \item{thres}{
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, q = 1, include.mean = T, coef0 = NULL, secoef0 = NULL, 
    fixed = NULL, prelim = F, details = F, thres = 2) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    nT = dim(da)[1]
    k = dim(da)[2]
    if (q < 1) 
        q = 1
    kq = k * q
    VMAdata <<- da
    EstStep <<- TRUE
    if (length(coef0) < 1) {
        m1 = VARorder(da, q + 12, output = FALSE)
        porder = m1$aicor
        if (porder < 1) 
            porder = 1
        m2 = VAR(da, porder, output = FALSE)
        y = da[(porder + 1):nT, ]
        x = m2$residuals
        m3 = THini(y, x, q, include.mean)
        beta = -m3$estimates
        sebeta = m3$se
        nr = dim(beta)[1]
        if (include.mean) {
            beta[1, ] = -beta[1, ]
        }
        if (prelim) {
            fixed = matrix(0, nr, k)
            for (j in 1:k) {
                tt = beta[, j]/sebeta[, j]
                idx = c(1:nr)[abs(tt) >= thres]
                fixed[idx, j] = 1
            }
        }
    }
    else {
        beta = coef0
        sebeta = secoef0
        nr = dim(beta)[1]
    }
    if (length(fixed) == 0) {
        fixed = matrix(1, nr, k)
    }
    fix1 <<- fixed
    inc.mean <<- include.mean
    MAq <<- q
    par = NULL
    separ = NULL
    VMAecnt = 0
    ist = 0
    if (include.mean) {
        jdx = c(1:k)[fix1[1, ] == 1]
        VMAecnt = length(jdx)
        if (VMAecnt > 0) {
            par = beta[1, jdx]
            separ = sebeta[1, jdx]
        }
        TH = beta[2:(kq + 1), ]
        seTH = sebeta[2:(kq + 1), ]
        ist = 1
    }
    else {
        TH = beta
        seTH = sebeta
    }
    for (j in 1:k) {
        idx = c(1:(nr - ist))[fix1[(ist + 1):nr, j] == 1]
        if (length(idx) > 0) {
            par = c(par, TH[idx, j])
            separ = c(separ, seTH[idx, j])
        }
    }
    ParE <<- par
    VMAecnt <<- VMAecnt
    cat("Number of parameters: ", length(par), "\n")
    cat("initial estimates: ", par, "\n")
    lowerBounds = par
    upperBounds = par
    npar = length(par)
    mult = 2
    if ((npar > 10) || (MAq > 2)) 
        mult = 1.5
    if (length(coef0) > 0) {
        mult = 1
    }
    for (j in 1:npar) {
        lowerBounds[j] = par[j] - mult * separ[j]
        upperBounds[j] = par[j] + mult * separ[j]
    }
    cat("Par. Lower-bounds: ", lowerBounds, "\n")
    cat("Par. Upper-bounds: ", upperBounds, "\n")
    if (details) {
        fit = nlminb(start = ParE, objective = EVMAq, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = ParE, objective = EVMAq, control = list(step.min = 0.4, 
            step.max = 0.8), lower = lowerBounds, upper = upperBounds)
    }
    est = fit$par
    zt <- VMAdata
    ist = 0
    mu = rep(0, k)
    icnt = 0
    if (include.mean) {
        ist = 1
        jdx = c(1:k)[fix1[1, ] == 1]
        icnt = length(jdx)
        if (icnt > 0) 
            mu[jdx] = est[1:icnt]
    }
    for (j in 1:k) {
        zt[, j] = zt[, j] - mu[j]
    }
    theta = matrix(0, kq, k)
    for (j in 1:k) {
        idx = c(1:kq)[fix1[(ist + 1):(ist + kq), j] == 1]
        jcnt = length(idx)
        if (jcnt > 0) {
            theta[idx, j] = est[(icnt + 1):(icnt + jcnt)]
            icnt = icnt + jcnt
        }
    }
    theta = t(theta)
    k1 = dim(theta)[2]
    Theta = theta[, 1:k]
    Zt = zt
    if (q > 1) {
        z0 = cbind(diag(1, k * (q - 1)), matrix(0, k * (q - 1), 
            k))
        Theta = rbind(theta, z0)
        Zt = cbind(zt, matrix(0, nT, k * (q - 1)))
    }
    m1 = eigen(Theta)
    V1 = m1$values
    M1 = Mod(V1)
    ich = 0
    for (i in 1:k1) {
        if (M1[i] > 1) {
            V1[i] = 1/V1[i]
            ich = 1
        }
    }
    est1 = NULL
    if (VMAecnt > 0) 
        est1 = est[1:VMAecnt]
    if (ich > 0) {
        P1 = m1$vectors
        P1i = solve(P1)
        Theta = Re(P1 \%*\% diag(V1) \%*\% P1i)
        beta = t(Theta[1:k, ])
        ist = 0
        nr = kq
        if (include.mean) {
            ist = 1
            nr = kq + 1
        }
        for (j in 1:k) {
            idx = c(1:k1)[fix1[(ist + 1):nr, j] == 1]
            if (length(idx) > 0) {
                est1 = c(est1, beta[idx, j])
            }
        }
        est = est1
    }
    cat("Final   Estimates: ", est, "\n")
    EstStep <- FALSE
    epsilon = 1e-04 * est
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
            Hessian[i, j] = (EVMAq(x1) - EVMAq(x2) - EVMAq(x3) + 
                EVMAq(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    se.coef = sqrt(diag(solve(Hessian)))
    tval = est/se.coef
    matcoef = cbind(est, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
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
    Past = matrix(0, 1, kq)
    at = NULL
    for (t in 1:nT) {
        tmp = zt[t, ] + Past \%*\% TH
        at = rbind(at, tmp)
        if (q == 1) {
            Past = tmp
        }
        else {
            Past = c(tmp, Past[1:(kq - k)])
        }
    }
    sig = t(at) \%*\% at/nT
    cat(" ", "\n")
    cat("Residuals cov-matrix:", "\n")
    print(sig)
    dd = det(sig)
    d1 = log(dd)
    aic = d1 + 2 * npar/nT
    bic = d1 + log(T) * npar/nT
    cat("----", "\n")
    cat("aic= ", aic, "\n")
    cat("bic= ", bic, "\n")
    Theta = t(TH)
    if (include.mean) {
        TH = rbind(cnt, TH)
        seTH = rbind(secnt, seTH)
    }
    VMAe <- list(data = da, MAorder = q, cnst = include.mean, 
        coef = TH, secoef = seTH, residuals = at, Sigma = sig, 
        Theta = Theta, mu = cnt, aic = aic, bic = bic)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
