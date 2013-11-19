\name{SCMfit}
\alias{SCMfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting a VARMA Model via SCM Approch
%%  ~~function to do ... ~~
}
\description{Perform estimation of a VARMA model specified via 
the SCM approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMfit(da, scms, Tdx, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 1, ref = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{The T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{scms}{A k-by-2 matrix of the orders of SCMs
%%     ~~Describe \code{scms} here~~
}
  \item{Tdx}{A k-dimensional vector for locating "1" of each row in the 
transformation matrix.
%%     ~~Describe \code{Tdx} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. 
Default is to include mean vector.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to set parameters to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch for preliminary estimation. 
Default is false. 
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control details of output
%%     ~~Describe \code{details} here~~
}
  \item{thres}{Threshold for individual t-ratio when setting parameters 
to zero. Default is 1.
%%     ~~Describe \code{thres} here~~
}
  \item{ref}{A switch to use SCMmod in model specification.
%%     ~~Describe \code{ref} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of estimation, including model specification and residuals
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{zt=read.table("flourc.txt")
ord=matrix(c(1,1,2,0,0,0),3,2)
Tdx=c(2,3,1)
m1=SCMfit(zt,ord,Tdx)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, scms, Tdx, include.mean = T, fixed = NULL, prelim = F, 
    details = F, thres = 1, ref = 0) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    nT = dim(da)[1]
    k = dim(da)[2]
    p = max(scms[, 1])
    q = max(scms[, 2])
    pq = max(p, q)
    kp = k * p
    kq = k * q
    pq = max(p, q)
    cat("Maximum VARMA order: (", p, ",", q, ")", "\n")
    if (ref < 1) {
        mm1 = SCMmod(scms, Tdx, FALSE)
        locTmtx <<- mm1$Tmtx
        locAR <<- mm1$ARpar
        locMA <<- mm1$MApar
        Vtsdata <<- da
        Tdx <<- Tdx
        inc.meanSCM <<- include.mean
        scms <<- scms
        cat("Locations of estimable parameters: Transformation Matrix", 
            "\n")
        print(locTmtx)
        cat("AR parameters", "\n")
        print(locAR)
        cat("MA parameters", "\n")
        print(locMA)
        m1 = VARorder(da, pq + 9, output = FALSE)
        porder = m1$aicor
        if (porder < 1) 
            porder = 1
        m2 = VAR(da, porder, output = FALSE)
        y = da[(porder + 1):nT, ]
        x = m2$residuals
        m3 = iniSCM(y, x, locTmtx, locAR, locMA)
        SCMpar <<- m3$par
        seSCMpar = m3$se
        nr = length(SCMpar)
        if (prelim) {
            fixed = rep(0, nr)
            for (j in 1:nr) {
                tt = SCMpar[j]/seSCMpar[j]
                if (abs(tt) >= thres) {
                  fixed[j] = 1
                }
                else {
                  SCMpar[j] = 0
                }
            }
        }
        if (length(fixed) < 1) 
            fixed = rep(1, nr)
    }
    fixSCM <<- fixed
    nr = length(SCMpar)
    JJdx = c(1:nr)[fixed == 1]
    par = SCMpar[JJdx]
    separ = seSCMpar[JJdx]
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
        fit = nlminb(start = par, objective = LLSCM, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = par, objective = LLSCM, lower = lowerBounds, 
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
            Hessian[i, j] = (LLSCM(x1) - LLSCM(x2) - LLSCM(x3) + 
                LLSCM(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    d1 = det(Hessian)
    if (d1 < 1e-13) {
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
    SCMpar[JJdx] = fit$par
    seSCMpar[JJdx] = se.coef
    Cnt = rep(0, k)
    seCnt = rep(0, k)
    Ph0 = locTmtx
    sePh0 = matrix(0, k, k)
    PH = NULL
    TH = NULL
    if (p > 0) {
        PH = matrix(0, k, kp)
        sePH = matrix(0, k, kp)
    }
    else {
        PH = NULL
        sePH = NULL
    }
    if (q > 0) {
        TH = matrix(0, k, kq)
        seTH = matrix(0, k, kq)
    }
    else {
        TH = NULL
        seTH = NULL
    }
    icnt = 0
    for (i in 1:k) {
        idx = NULL
        jdx = NULL
        kdx = NULL
        if (p > 0) 
            idx = c(1:kp)[locAR[i, ] > 1]
        if (q > 0) 
            jdx = c(1:kq)[locMA[i, ] > 1]
        kdx = c(1:k)[locTmtx[i, ] > 1]
        iend = length(idx)
        jend = length(jdx)
        kend = length(kdx)
        if (inc.meanSCM) {
            icnt = icnt + 1
            Cnt[i] = SCMpar[icnt]
            seCnt[i] = seSCMpar[icnt]
        }
        if (kend > 0) {
            Ph0[i, kdx] = SCMpar[(icnt + 1):(icnt + kend)]
            sePh0[i, kdx] = seSCMpar[(icnt + 1):(icnt + kend)]
            icnt = icnt + kend
        }
        if (iend > 0) {
            PH[i, idx] = SCMpar[(icnt + 1):(icnt + iend)]
            sePH[i, idx] = seSCMpar[(icnt + 1):(icnt + iend)]
            icnt = icnt + iend
        }
        if (jend > 0) {
            TH[i, jdx] = SCMpar[(icnt + 1):(icnt + jend)]
            seTH[i, jdx] = seSCMpar[(icnt + 1):(icnt + jend)]
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
    jcnt = 0
    if (p > 0) {
        cat("AR coefficient matrix", "\n")
        for (i in 1:p) {
            cat("AR(", i, ")-matrix", "\n")
            ph = PH[, (jcnt + 1):(jcnt + k)]
            print(round(ph, 3))
            jcnt = jcnt + k
        }
    }
    if (q > 0) {
        cat("MA coefficient matrix", "\n")
        icnt = 0
        for (i in 1:q) {
            cat("MA(", i, ")-matrix", "\n")
            theta = -TH[, (icnt + 1):(icnt + k)]
            print(round(theta, 3))
            icnt = icnt + k
        }
    }
    Ph0i = solve(Ph0)
    ARc = NULL
    MAc = NULL
    if (p > 0) 
        ARc = Ph0i \%*\% PH
    if (q > 0) 
        MAc = Ph0i \%*\% TH
    Cntc = Ph0i \%*\% as.matrix(Cnt, k, 1)
    zt = Vtsdata
    at = matrix((zt[1, ] - Cntc), 1, k)
    if (pq > 1) {
        for (t in 2:pq) {
            tmp = matrix((zt[t, ] - Cntc), 1, k)
            if (p > 0) {
                for (j in 1:p) {
                  if ((t - j) > 0) {
                    jdx = (j - 1) * k
                    tmp1 = matrix(zt[(t - j), ], 1, k) \%*\% t(as.matrix(ARc[, 
                      (jdx + 1):(jdx + k)]))
                    tmp = tmp - tmp1
                  }
                }
            }
            if (q > 0) {
                for (j in 1:q) {
                  jdx = (j - 1) * k
                  if ((t - j) > 0) {
                    tmp2 = matrix(at[(t - j), ], 1, k) \%*\% t(as.matrix(MAc[, 
                      (jdx + 1):(jdx + k)]))
                    tmp = tmp - tmp2
                  }
                }
            }
            at = rbind(at, tmp)
        }
    }
    ist = pq + 1
    Pcnt = NULL
    beta = NULL
    if (inc.meanSCM) {
        beta = matrix(Cntc, 1, k)
        Pcnt = c(1)
    }
    if (length(ARc) > 0) 
        beta = rbind(beta, t(ARc))
    if (length(MAc) > 0) 
        beta = rbind(beta, t(MAc))
    idim = k * (p + q)
    if (inc.meanSCM) {
        Pcnt = c(1)
        idim = idim + 1
    }
    for (t in ist:nT) {
        Past = NULL
        if (p > 0) {
            for (j in 1:p) {
                Past = c(Past, zt[(t - j), ])
            }
        }
        if (q > 0) {
            for (j in 1:q) {
                Past = c(Past, at[(t - j), ])
            }
        }
        tmp = matrix(c(Pcnt, Past), 1, idim) \%*\% beta
        tmp3 = zt[t, ] - tmp
        at = rbind(at, tmp3)
    }
    at = as.matrix(at[(ist:nT), ])
    sig = crossprod(at, at)/(nT - pq)
    cat(" ", "\n")
    cat("Residuals cov-matrix:", "\n")
    print(sig)
    dd = det(sig)
    d1 = log(dd)
    jj = 0
    for (i in 1:k) {
        kdx = c(1:k)[locTmtx[i, ] > 1]
        jj = jj + length(kdx)
    }
    aic = d1 + 2 * (npar - jj)/nT
    bic = d1 + log(nT) * (npar - jj)/nT
    cat("----", "\n")
    cat("aic= ", aic, "\n")
    cat("bic= ", bic, "\n")
    TH1 = NULL
    if (length(TH) > 0) 
        TH1 = -TH
    SCMfit <- list(data = da, SCMs = scms, Tdx = Tdx, locTmtx = locTmtx, 
        locAR = locAR, locMA = locMA, cnst = include.mean, coef = SCMpar, 
        secoef = seSCMpar, residuals = at, Sigma = sig, aic = aic, 
        bic = bic, Ph0 = Ph0, Phi = PH, Theta = TH1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
