\name{Kronfit}
\alias{Kronfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of VARMA Model via Kronecker Index
%%  ~~function to do ... ~~
}
\description{Perform estimation of a VARMA model specified by Kronecker 
indices
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Kronfit(da, kidx, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix (T-by-k) of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{kidx}{The set of Kronecker indices
%%     ~~Describe \code{kidx} here~~
}
  \item{include.mean}{A logical switch for including the mean vector in 
estimation. Default is to include the mean vector.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used to set zero parameter constraints. 
This is used mainly in the command refKronfit.
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A lofical switch for a preliminary estimation.
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control output.
%%     ~~Describe \code{details} here~~
}
  \item{thres}{A threshold for t-ratio in setting parameter to zero. 
Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of estimation, including details of model specification and 
residuals
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 4)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{refKronfit
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{zt=read.table("flourc.txt")
kdx=c(1,1,1)
m1=Kronfit(da,kdx)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, kidx, include.mean = T, fixed = NULL, prelim = F, 
    details = F, thres = 1) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    nT = dim(da)[1]
    k = dim(da)[2]
    if (length(fixed) < 1) {
        k1 = length(kidx)
        if (k1 <= 0) {
            k1 = 1
            kidx = c(1)
        }
        maxk = kidx[1]
        for (i in 1:k1) {
            if (maxk < kidx[i]) 
                maxk = kidx[i]
        }
        Vtsdata <<- da
        maxK <<- maxk
        m0 = Kronspec(kidx, output = F)
        ARid = m0$PhiID
        MAid = m0$ThetaID
        print(ARid)
        print(MAid)
        m1 = VARorder(da, maxk + 9, output = FALSE)
        porder = m1$aicor
        if (porder < 1) 
            porder = 1
        m2 = VAR(da, porder, output = FALSE)
        y = da[(porder + 1):nT, ]
        x = m2$residuals
        m3 = iniKro(y, x, ARid, MAid, include.mean)
        Kpar <<- m3$par
        seKpar = m3$se
        nr = length(Kpar)
        if (prelim) {
            fixed = rep(0, nr)
            for (j in 1:nr) {
                tt = Kpar[j]/seKpar[j]
                if (abs(tt) >= thres) {
                  fixed[j] = 1
                }
                else {
                  Kpar[j] = 0
                }
            }
        }
        if (length(fixed) < 1) {
            fixed = rep(1, nr)
        }
        fixK <<- fixed
        inc.meanK <<- include.mean
        ARid <<- ARid
        MAid <<- MAid
    }
    maxk <- maxK
    nr = length(Kpar)
    JJdx = c(1:nr)[fixed == 1]
    par = Kpar[JJdx]
    separ = seKpar[JJdx]
    JJdx <<- JJdx
    cat("Number of parameters: ", length(par), "\n")
    cat("initial estimates: ", round(par, 4), "\n")
    lowerBounds = par
    upperBounds = par
    for (j in 1:length(par)) {
        lowerBounds[j] = par[j] - 2 * separ[j]
        upperBounds[j] = par[j] + 2 * separ[j]
    }
    cat("Upper-bound: ", round(upperBounds, 4), "\n")
    cat("Lower-bound: ", round(lowerBounds, 4), "\n")
    if (details) {
        fit = nlminb(start = par, objective = LLKron, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = par, objective = LLKron, lower = lowerBounds, 
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
            Hessian[i, j] = (LLKron(x1) - LLKron(x2) - LLKron(x3) + 
                LLKron(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    d1 = det(Hessian)
    if (d1 < 1e-10) {
        se.coef = rep(1, npar)
    }
    else {
        se.coef = sqrt(diag(solve(Hessian)))
    }
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", 
        " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
    Kpar[JJdx] = fit$par
    seKpar[JJdx] = se.coef
    Cnt = rep(0, k)
    seCnt = rep(0, k)
    Ph0 = diag(rep(1, k))
    sePh0 = diag(rep(1, k))
    kp1 = dim(ARid)[2]
    kp = kp1 - k
    PH = matrix(0, k, kp)
    sePH = matrix(0, k, kp)
    TH = matrix(0, k, kp)
    seTH = matrix(0, k, kp)
    icnt = 0
    for (i in 1:k) {
        idx = c(1:kp1)[ARid[i, ] > 1]
        jdx = c(1:kp1)[MAid[i, ] > 1]
        kdx = c(1:k)[ARid[i, 1:k] > 1]
        if (length(kdx) > 0) {
            idx = idx[-kdx]
            jdx = jdx[-kdx]
        }
        iend = length(idx)
        jend = length(jdx)
        kend = length(kdx)
        if (inc.meanK) {
            icnt = icnt + 1
            Cnt[i] = Kpar[icnt]
            seCnt[i] = seKpar[icnt]
        }
        if (kend > 0) {
            Ph0[i, kdx] = Kpar[(icnt + 1):(icnt + kend)]
            sePh0[i, kdx] = seKpar[(icnt + 1):(icnt + kend)]
            icnt = icnt + kend
        }
        if (iend > 0) {
            PH[i, idx - k] = Kpar[(icnt + 1):(icnt + iend)]
            sePH[i, idx - k] = seKpar[(icnt + 1):(icnt + iend)]
            icnt = icnt + iend
        }
        if (jend > 0) {
            TH[i, jdx - k] = Kpar[(icnt + 1):(icnt + jend)]
            seTH[i, jdx - k] = seKpar[(icnt + 1):(icnt + jend)]
            icnt = icnt + jend
        }
    }
    cat("---", "\n")
    cat("Estimates in matrix form:", "\n")
    if (include.mean) {
        cat("Constant term: ", "\n")
        cat("Estimates: ", round(Cnt, 3), "\n")
    }
    cat("AR and MA lag-0 coefficient matrix", "\n")
    print(round(Ph0, 3))
    cat("AR coefficient matrix", "\n")
    jcnt = 0
    for (i in 1:maxk) {
        cat("AR(", i, ")-matrix", "\n")
        ph = PH[, (jcnt + 1):(jcnt + k)]
        print(round(ph, 3))
        jcnt = jcnt + k
    }
    cat("MA coefficient matrix", "\n")
    icnt = 0
    for (i in 1:maxk) {
        cat("MA(", i, ")-matrix", "\n")
        theta = -TH[, (icnt + 1):(icnt + k)]
        print(round(theta, 3))
        icnt = icnt + k
    }
    Ph0i = solve(Ph0)
    ARc = Ph0i \%*\% PH
    MAc = Ph0i \%*\% TH
    Cntc = Ph0i \%*\% as.matrix(Cnt, k, 1)
    zt = Vtsdata
    ist = maxk + 1
    at = matrix((zt[1, ] - Cntc), 1, k)
    if (maxk > 1) {
        for (t in 2:maxk) {
            tmp = matrix((zt[t, ] - Cntc), 1, k)
            for (j in 1:maxk) {
                if ((t - j) > 0) {
                  jdx = (j - 1) * k
                  tmp1 = matrix(zt[(t - j), ], 1, k) \%*\% t(as.matrix(ARc[, 
                    (jdx + 1):(jdx + k)]))
                  tmp = tmp - tmp1
                }
            }
            for (j in 1:maxk) {
                jdx = (j - 1) * k
                if ((t - j) > 0) {
                  tmp2 = matrix(at[(t - j), ], 1, k) \%*\% t(as.matrix(MAc[, 
                    (jdx + 1):(jdx + k)]))
                  tmp = tmp - tmp2
                }
            }
            at = rbind(at, tmp)
        }
    }
    ist = maxk + 1
    Pcnt = NULL
    beta = NULL
    if (inc.meanK) {
        beta = matrix(Cntc, 1, k)
        Pcnt = c(1)
    }
    beta = rbind(beta, t(ARc), t(MAc))
    idim = k * maxk * 2
    if (inc.meanK) {
        Pcnt = c(1)
        idim = idim + 1
    }
    for (t in (maxk + 1):nT) {
        Past = NULL
        for (j in 1:maxk) {
            Past = c(Past, zt[(t - j), ])
        }
        for (j in 1:maxk) {
            Past = c(Past, at[(t - j), ])
        }
        tmp = matrix(c(Pcnt, Past), 1, idim) \%*\% beta
        tmp3 = zt[t, ] - tmp
        at = rbind(at, tmp3)
    }
    at = at[(ist:nT), ]
    sig = crossprod(at, at)/(nT - maxk)
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
    Kronfit <- list(data = da, Kindex = kidx, ARid = ARid, MAid = MAid, 
        cnst = include.mean, coef = Kpar, secoef = seKpar, residuals = at, 
        Sigma = sig, aic = aic, bic = bic, Ph0 = Ph0, Phi = PH, 
        Theta = -TH)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
