\name{sVARMA}
\alias{sVARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Seasonal VARMA Model
%%  ~~function to do ... ~~
}
\description{Performs conditional maximum likelihood estimation of 
a seasonal VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sVARMA(da, order, sorder, s, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 0.8, switch = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k data matrix of a k-dimensional seasonal time series
%%     ~~Describe \code{da} here~~
}
  \item{order}{Regular order (p,d,q) of the model
%%     ~~Describe \code{order} here~~
}
  \item{sorder}{Seasonal order (P,D,Q) of the model
%%     ~~Describe \code{sorder} here~~
}
  \item{s}{Seasonality. s=4 for quarterly data and s=12 for monthly series
%%     ~~Describe \code{s} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. Deafult is 
to include the mean
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to set zero parameter constraints
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A lofical switch for preliminary estimation. 
Default is none.
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch for output
%%     ~~Describe \code{details} here~~
}
  \item{thres}{Threshold value for individual t-ratio used to set zero parameter constraints
%%     ~~Describe \code{thres} here~~
}
  \item{switch}{A logical switch to exchange the ordering of the regular 
and seasonal VMA factors. Default is theta(B)*Theta(B).
%%     ~~Describe \code{switch} here~~
}
}
\details{Estimation of a sesonal VARMA model
%%  ~~ If necessary, more details than the description above ~~
}
\value{Estimation results of a seasonal VARMA model, including model specification
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
da=read.table("m-hsoldhst6312.txt",header=T)
zt=log(da[,3:4])
m3=sVARMA(zt,order=c(0,1,2),sorder=c(0,1,1),s=12)
names(m3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, order, sorder, s, include.mean = T, fixed = NULL, 
    prelim = F, details = F, thres = 0.8, switch = F) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    p = order[1]
    d = order[2]
    q = order[3]
    P = sorder[1]
    D = sorder[2]
    Q = sorder[3]
    nT = dim(da)[1]
    k = dim(da)[2]
    if (p < 0) 
        p = 0
    if (q < 0) 
        q = 0
    if (P < 0) 
        P = 0
    if (Q < 0) 
        Q = 0
    if (s < 0) 
        s = -s
    if (d > 1) {
        cat("Regular difference is adjusted to d=1", "\n")
        d = 1
    }
    if (D > 1) {
        cat("Seasonal difference is adjusted to D=1", "\n")
        D = 1
    }
    kp = k * p
    kq = k * q
    kP = k * P
    kQ = k * Q
    if (d == 1) {
        X = NULL
        MEAN = rep(0, k)
        for (j in 1:k) {
            X = cbind(X, diff(da[, j]))
            t1 = t.test(X[, j])
            if (t1$p.value < 0.05) 
                MEAN[j] = 1
        }
        if (sum(MEAN) < 1) 
            include.mean = F
    }
    else {
        X = da
    }
    if (D == 1) {
        DX = NULL
        Smean = rep(0, k)
        for (j in 1:k) {
            DX = cbind(DX, diff(X[, j], s))
            t1 = t.test(DX[, j])
            if (t1$p.value < 0.05) 
                Smean[j] = 1
        }
        if (sum(Smean) < 1) 
            include.mean = F
    }
    else {
        DX = X
    }
    nT = dim(DX)[1]
    arlags = NULL
    if (p > 0) {
        arlags = c(1:p)
        if (P > 0) 
            arlags = c(arlags, c(1:P) * s, c(1:P) * s + c(1:p))
    }
    else {
        if (P > 0) 
            arlags = c(1:P) * s
    }
    malags = NULL
    if (q > 0) {
        malags = c(1:q)
        if (Q > 0) 
            malags = c(malags, c(1:Q) * s, c(1:Q) * s + c(1:q))
    }
    else {
        if (Q > 0) 
            malags = c(1:Q) * s
    }
    nar = length(arlags)
    nma = length(malags)
    idim = k * (nar + nma)
    if (include.mean) 
        idim = idim + 1
    if (length(fixed) == 0) 
        fixed = matrix(1, idim, k)
    Mtdata <<- DX
    ARlags <<- arlags
    MAlags <<- malags
    Order <<- c(order, sorder)
    inc.mean <<- include.mean
    fix1 <<- fixed
    swi <<- switch
    phi = NULL
    sphi = NULL
    sephi = NULL
    sesphi = NULL
    if (p > 0) 
        phi = matrix(0, k, k * p)
    sephi = phi
    if (P > 0) 
        sphi = matrix(0, k, k * P)
    sesphi = sphi
    theta = NULL
    stheta = NULL
    setheta = NULL
    sestheta = NULL
    if (q > 0) 
        theta = matrix(0, k, k * q)
    setheta = theta
    if (Q > 0) 
        stheta = matrix(0, k, k * Q)
    sestheta = stheta
    resi = NULL
    for (j in 1:k) {
        m1 = arima(DX[, j], order = c(p, 0, q), seasona = list(order = c(P, 
            0, Q), period = s))
        resi = cbind(resi, m1$residuals)
        seest = sqrt(diag(m1$var.coef))
        icnt = 0
        if (p > 0) {
            for (i in 1:p) {
                icnt = icnt + 1
                ii = (i - 1) * k
                phi[j, (ii + j)] = m1$coef[icnt]
                sephi[j, (ii + j)] = seest[icnt]
            }
        }
        if (q > 0) {
            for (i in 1:q) {
                ii = (i - 1) * k
                icnt = icnt + 1
                theta[j, (ii + j)] = -m1$coef[icnt]
                setheta[j, (ii + j)] = seest[icnt]
            }
        }
        if (P > 0) {
            for (i in 1:P) {
                icnt = icnt + 1
                ii = (i - 1) * k
                sphi[j, (ii + j)] = m1$coef[icnt]
                sesphi[j, (ii + j)] = seest[icnt]
            }
        }
        if (Q > 0) {
            for (i in 1:Q) {
                ii = (i - 1) * k
                icnt = icnt + 1
                stheta[j, (ii + j)] = -m1$coef[icnt]
                sestheta[j, (ii + j)] = seest[icnt]
            }
        }
    }
    m2 = siniEST(DX, resi, arlags, malags, include.mean)
    beta = t(m2$estimates)
    sebeta = t(m2$se)
    icnst = 0
    if (include.mean) 
        icnst = 1
    if (nar > 0) {
        if (p > 0) {
            for (i in 1:p) {
                idx = (i - 1) * k
                for (ii in 1:k) {
                  jdx = c(1:k)[-ii]
                  phi[jdx, (idx + ii)] = beta[jdx, (icnst + idx + 
                    ii)]
                  sephi[jdx, (idx + ii)] = sebeta[jdx, (icnst + 
                    idx + ii)]
                }
            }
        }
        if (P > 0) {
            for (i in 1:P) {
                kdx = (i - 1) * k
                idx = k * p + kdx
                for (ii in 1:k) {
                  jdx = c(1:k)[-ii]
                  sphi[jdx, (kdx + ii)] = beta[jdx, (icnst + 
                    idx + ii)]
                  sesphi[jdx, (kdx + ii)] = sebeta[jdx, (icnst + 
                    idx + ii)]
                }
            }
        }
    }
    if (nma > 0) {
        if (q > 0) {
            for (i in 1:q) {
                kdx = (i - 1) * k
                idx = nar * k + kdx
                for (ii in 1:k) {
                  jdx = c(1:k)[-ii]
                  theta[jdx, (kdx + ii)] = -beta[jdx, (icnst + 
                    idx + ii)]
                  setheta[jdx, (kdx + ii)] = sebeta[jdx, (icnst + 
                    idx + ii)]
                }
            }
        }
        if (Q > 0) {
            for (i in 1:Q) {
                kdx = (i - 1) * k
                idx = (nar + q) * k + kdx
                for (ii in 1:k) {
                  jdx = c(1:k)[-ii]
                  stheta[jdx, (kdx + ii)] = -beta[jdx, (icnst + 
                    idx + ii)]
                  sestheta[jdx, (kdx + ii)] = sebeta[jdx, (icnst + 
                    idx + ii)]
                }
            }
        }
    }
    par = NULL
    separ = NULL
    ist = 0
    if (include.mean) {
        jdx = c(1:k)[fix1[1, ] == 1]
        if (length(jdx) > 0) {
            par = beta[jdx, 1]
            separ = sebeta[jdx, 1]
        }
        ist = 1
    }
    if (nar > 0) {
        if (p > 0) {
            for (j in 1:k) {
                idx = c(1:kp)[fix1[(ist + 1):(ist + kp), j] == 
                  1]
                if (length(idx) > 0) {
                  par = c(par, phi[j, idx])
                  separ = c(separ, sephi[j, idx])
                }
            }
            ist = ist + kp
        }
        if (P > 0) {
            for (j in 1:k) {
                idx = c(1:kP)[fix1[(ist + 1):(ist + kP), j] == 
                  1]
                if (length(idx) > 0) {
                  par = c(par, sphi[j, idx])
                  separ = c(separ, sesphi[j, idx])
                }
            }
            ist = ist + kP
        }
    }
    if (nma > 0) {
        if (q > 0) {
            for (j in 1:k) {
                idx = c(1:kq)[fix1[(ist + 1):(ist + kq), j] == 
                  1]
                if (length(idx) > 0) {
                  par = c(par, theta[j, idx])
                  separ = c(separ, setheta[j, idx])
                }
            }
            ist = ist + kq
        }
        if (Q > 0) {
            for (j in 1:k) {
                idx = c(1:kQ)[fix1[(ist + 1):(ist + kQ), j] == 
                  1]
                if (length(idx) > 0) {
                  par = c(par, stheta[j, idx])
                  separ = c(separ, sestheta[j, idx])
                }
            }
        }
    }
    jst = max(arlags[nar], malags[nma])
    Sresi <<- resi[1:jst, ]
    cat("Number of parameters: ", length(par), "\n")
    cat("initial estimates: ", par, "\n")
    lowerBounds = par
    upperBounds = par
    for (j in 1:length(par)) {
        lowerBounds[j] = par[j] - 2 * separ[j]
        upperBounds[j] = par[j] + 2 * separ[j]
    }
    if (details) {
        fit = nlminb(start = par, objective = LLKsvarma, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = par, objective = LLKsvarma, lower = lowerBounds, 
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
            Hessian[i, j] = (LLKsvarma(x1) - LLKsvarma(x2) - 
                LLKsvarma(x3) + LLKsvarma(x4))/(4 * epsilon[i] * 
                epsilon[j])
        }
    }
    se.coef = sqrt(diag(solve(Hessian)))
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", 
        " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
    est = fit$par
    ist = 0
    icnt = 0
    Ph0 = rep(0, k)
    sePh0 = rep(0, k)
    beta = NULL
    sebeta = NULL
    if (include.mean) {
        idx = c(1:k)[fix1[1, ] == 1]
        icnt = length(idx)
        if (icnt > 0) {
            Ph0[idx] = est[1:icnt]
            sePh0[idx] = se.coef[1:icnt]
        }
        ist = 1
        beta = rbind(beta, Ph0)
        sebeta = rbind(sebeta, sePh0)
    }
    PH = NULL
    sePH = NULL
    sPH = NULL
    sesPH = NULL
    if (p > 0) {
        PH = matrix(0, kp, k)
        sePH = matrix(0, kp, k)
        for (j in 1:k) {
            idx = c(1:kp)[fix1[(ist + 1):(ist + kp), j] == 1]
            jdx = length(idx)
            if (jdx > 0) {
                PH[idx, j] = est[(icnt + 1):(icnt + jdx)]
                sePH[idx, j] = se.coef[(icnt + 1):(icnt + jdx)]
                icnt = icnt + jdx
            }
        }
        ist = ist + kp
        beta = rbind(beta, PH)
        sebeta = rbind(sebeta, sePH)
    }
    if (P > 0) {
        sPH = matrix(0, kP, k)
        sesPH = matrix(0, kP, k)
        for (j in 1:k) {
            idx = c(1:kP)[fix1[(ist + 1):(ist + kP), j] == 1]
            jdx = length(idx)
            if (jdx > 0) {
                sPH[idx, j] = est[(icnt + 1):(icnt + jdx)]
                sesPH[idx, j] = se.coef[(icnt + 1):(icnt + jdx)]
                icnt = icnt + jdx
            }
        }
        ist = ist + kP
        beta = rbind(beta, sPH)
        sebeta = rbind(sebeta, sesPH)
    }
    TH = NULL
    seTH = NULL
    sTH = NULL
    sesTH = NULL
    if (q > 0) {
        TH = matrix(0, kq, k)
        seTH = matrix(0, kq, k)
        for (j in 1:k) {
            idx = c(1:kq)[fix1[(ist + 1):(ist + kq), j] == 1]
            jdx = length(idx)
            if (jdx > 0) {
                TH[idx, j] = est[(icnt + 1):(icnt + jdx)]
                seTH[idx, j] = se.coef[(icnt + 1):(icnt + jdx)]
                icnt = icnt + jdx
            }
        }
        ist = ist + kq
        beta = rbind(beta, -TH)
        sebeta = rbind(sebeta, seTH)
    }
    if (Q > 0) {
        sTH = matrix(0, kQ, k)
        sesTH = matrix(0, kQ, k)
        for (j in 1:k) {
            idx = c(1:kQ)[fix1[(ist + 1):(ist + kQ), j] == 1]
            jdx = length(idx)
            if (jdx > 0) {
                sTH[idx, j] = est[(icnt + 1):(icnt + jdx)]
                sesTH[idx, j] = se.coef[(icnt + 1):(icnt + jdx)]
                icnt = icnt + jdx
            }
        }
        beta = rbind(beta, -sTH)
        sebeta = rbind(sebeta, sesTH)
    }
    cat("---", "\n")
    cat("Estimates in matrix form:", "\n")
    if (include.mean) {
        cat("Constant term: ", "\n")
        cat("Estimates: ", Ph0, "\n")
    }
    if (p > 0) {
        cat("Regular AR coefficient matrix", "\n")
        jcnt = 0
        for (i in 1:p) {
            cat("AR(", i, ")-matrix", "\n")
            ph = t(PH[(jcnt + 1):(jcnt + k), ])
            print(ph, digits = 3)
            jcnt = jcnt + k
        }
    }
    if (P > 0) {
        cat("Seasonal AR coefficient matrix", "\n")
        jcnt = 0
        for (i in 1:P) {
            cat("AR(", i * s, ")-matrix", "\n")
            ph = t(sPH[(jcnt + 1):(jcnt + k), ])
            print(ph, digits = 3)
            jcnt = jcnt + k
        }
    }
    if (q > 0) {
        cat("Regular MA coefficient matrix", "\n")
        icnt = 0
        for (i in 1:q) {
            cat("MA(", i, ")-matrix", "\n")
            the = t(TH[(icnt + 1):(icnt + k), ])
            print(the, digits = 3)
            icnt = icnt + k
        }
    }
    if (Q > 0) {
        cat("Seasonal MA coefficient matrix", "\n")
        icnt = 0
        for (i in 1:Q) {
            cat("MA(", i * s, ")-matrix", "\n")
            the = t(sTH[(icnt + 1):(icnt + k), ])
            print(the, digists = 3)
            icnt = icnt + k
        }
    }
    if ((p > 0) && (P > 0)) {
        if (switch) {
            Phi = t(Mtxprod1(t(PH), t(sPH), p, P))
        }
        else {
            Phi = t(Mtxprod(t(PH), t(sPH), p, P))
        }
    }
    if ((p > 0) && (P == 0)) 
        Phi = PH
    if ((p == 0) && (P > 0)) 
        Phi = sPH
    if ((q > 0) && (Q > 0)) {
        if (switch) {
            Theta = t(Mtxprod1(t(TH), t(sTH), q, Q))
        }
        else {
            Theta = t(Mtxprod(t(TH), t(sTH), q, Q))
        }
    }
    if ((q > 0) && (Q == 0)) 
        Theta = TH
    if ((q == 0) && (Q > 0)) 
        Theta = sTH
    zt = Mtdata
    pqmax = max(ARlags[nar], MAlags[nma])
    ist = pqmax + 1
    at = Sresi[1:pqmax, ]
    for (t in ist:nT) {
        tmp = zt[t, ] - Ph0
        if (nar > 0) {
            for (j in 1:nar) {
                jj = ARlags[j]
                jdx = (j - 1) * k
                ph = Phi[(jdx + 1):(jdx + k), ]
                tmp = tmp - matrix(zt[(t - jj), ], 1, k) \%*\% 
                  ph
            }
        }
        if (nma > 0) {
            for (j in 1:nma) {
                jj = MAlags[j]
                jdx = (j - 1) * k
                th = Theta[(jdx + 1):(jdx + k), ]
                tmp = tmp + matrix(at[(t - jj), ], 1, k) \%*\% 
                  th
            }
        }
        at = rbind(at, tmp)
    }
    at = at[(ist:nT), ]
    c1 = rep("resi", k)
    colnames(at) <- c1
    sig = t(at) \%*\% at/(nT - pqmax)
    cat(" ", "\n")
    cat("Residuals cov-matrix:", "\n")
    print(sig)
    dd = det(sig)
    d1 = log(dd)
    aic = d1 + 2 * npar/nT
    bic = d1 + log(nT) * npar/nT
    cat("----", "\n")
    cat("aic= ", round(aic, 4), "\n")
    cat("bic= ", round(bic, 4), "\n")
    if (length(PH) > 0) 
        PH = t(PH)
    if (length(sPH) > 0) 
        sPH = t(sPH)
    if (length(TH) > 0) 
        TH = t(TH)
    if (length(sTH) > 0) 
        sTH = t(sTH)
    sVARMA <- list(data = da, order = order, sorder = sorder, 
        period = s, cnst = include.mean, coef = beta, secoef = sebeta, 
        residuals = at, Sigma = sig, aic = aic, bic = bic, regPhi = PH, 
        seaPhi = sPH, regTheta = TH, seaTheta = sTH, Ph0 = Ph0, 
        switch = switch)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
