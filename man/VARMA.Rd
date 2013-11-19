\name{VARMA}
\alias{VARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vector Autoregressive Moving-Average Models
%%  ~~function to do ... ~~
}
\description{Performs conditional maximum likelihood estimation of 
a VARMA model. Multivariate Gaussian likelihood function is used.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMA(da, p = 0, q = 0, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix (T-by-k) of a k-dimensional time series with sample 
size T.
%%     ~~Describe \code{da} here~~
}
  \item{p}{AR order
%%     ~~Describe \code{p} here~~
}
  \item{q}{MA order
%%     ~~Describe \code{q} here~~
}
  \item{include.mean}{A logical switch to control estimation of the mean vector.
Dafault is to include the mean in estimation.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to control zero coefficients in estimation. 
It is mainly used by the command refVARMA.
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch to control preliminary estimation. 
Deafult is none.
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control the amount of output.
%%     ~~Describe \code{details} here~~
}
  \item{thres}{A threshold used to set zero parameter constraints based on 
individual t-ratio.
Default is 2.
%%     ~~Describe \code{thres} here~~
}
}
\details{Coefficient estimates, model specification, and residuals are returned.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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

\seealso{refVARMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
d1=read.table("m-pce.txt",header=T); d2=read.table("m-dspi.txt",header=T)
z1=diff(log(d1$pce)); z2=diff(log(d2$dspi))
zt=cbind(z1,z2)
colnames(zt) <- c("pce","dspi")
m2=VARMA(zt,p=3,q=1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, p = 0, q = 0, include.mean = T, fixed = NULL, prelim = F, 
    details = F, thres = 2) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    nT = dim(da)[1]
    k = dim(da)[2]
    if (p < 0) 
        p = 0
    if (q < 0) 
        q = 0
    if ((p + q) < 1) 
        p = 1
    pqmax = max(p, q)
    kq = k * q
    kp = k * p
    Vtsdata <<- da
    m1 = VARorder(da, p + q + 9, output = FALSE)
    porder = m1$aicor
    if (porder < 1) 
        porder = 1
    m2 = VAR(da, porder, output = FALSE)
    y = da[(porder + 1):nT, ]
    x = m2$residuals
    m3 = iniEST(y, x, p, q, include.mean)
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
    ARp <<- p
    MAq <<- q
    par = NULL
    separ = NULL
    ist = 0
    if (include.mean) {
        jdx = c(1:k)[fix1[1, ] == 1]
        if (length(jdx) > 0) {
            par = beta[1, jdx]
            separ = sebeta[1, jdx]
        }
        ist = 1
    }
    if (p > 0) {
        for (j in 1:k) {
            idx = c(1:kp)[fix1[(ist + 1):(ist + kp), j] == 1]
            if (length(idx) > 0) {
                tmp = beta[(ist + 1):(ist + kp), j]
                setmp = sebeta[(ist + 1):(ist + kp), j]
                par = c(par, tmp[idx])
                separ = c(separ, setmp[idx])
            }
        }
        ist = ist + kp
    }
    if (q > 0) {
        for (j in 1:k) {
            idx = c(1:kq)[fix1[(ist + 1):(ist + kq), j] == 1]
            if (length(idx) > 0) {
                tmp = beta[(ist + 1):(ist + kq), j]
                setmp = sebeta[(ist + 1):(ist + kq), j]
                par = c(par, tmp[idx])
                separ = c(separ, setmp[idx])
            }
        }
    }
    cat("Number of parameters: ", length(par), "\n")
    cat("initial estimates: ", round(par, 4), "\n")
    lowerBounds = par
    upperBounds = par
    for (j in 1:length(par)) {
        lowerBounds[j] = par[j] - 2 * separ[j]
        upperBounds[j] = par[j] + 2 * separ[j]
    }
    cat("Par. lower-bounds: ", round(lowerBounds, 4), "\n")
    cat("Par. upper-bounds: ", round(upperBounds, 4), "\n")
    if (details) {
        fit = nlminb(start = par, objective = LLKvarma, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = par, objective = LLKvarma, lower = lowerBounds, 
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
            Hessian[i, j] = (LLKvarma(x1) - LLKvarma(x2) - LLKvarma(x3) + 
                LLKvarma(x4))/(4 * epsilon[i] * epsilon[j])
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
    TH = NULL
    seTH = NULL
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
        beta = rbind(beta, TH)
        sebeta = rbind(sebeta, seTH)
    }
    cat("---", "\n")
    cat("Estimates in matrix form:", "\n")
    if (include.mean) {
        cat("Constant term: ", "\n")
        cat("Estimates: ", Ph0, "\n")
    }
    if (p > 0) {
        cat("AR coefficient matrix", "\n")
        jcnt = 0
        for (i in 1:p) {
            cat("AR(", i, ")-matrix", "\n")
            ph = t(PH[(jcnt + 1):(jcnt + k), ])
            print(ph, digits = 3)
            jcnt = jcnt + k
        }
    }
    if (q > 0) {
        cat("MA coefficient matrix", "\n")
        icnt = 0
        for (i in 1:q) {
            cat("MA(", i, ")-matrix", "\n")
            theta = -t(TH[(icnt + 1):(icnt + k), ])
            print(theta, digits = 3)
            icnt = icnt + k
        }
    }
    zt = Vtsdata
    ist = pqmax + 1
    at = matrix((zt[1, ] - Ph0), 1, k)
    if (pqmax > 1) {
        for (t in 2:pqmax) {
            tmp = matrix((zt[t, ] - Ph0), 1, k)
            if (p > 0) {
                for (j in 1:p) {
                  if ((t - j) > 0) {
                    jdx = (j - 1) * k
                    tmp1 = matrix(zt[(t - j), ], 1, k) \%*\% as.matrix(PH[(jdx + 
                      1):(jdx + k), ])
                    tmp = tmp - tmp1
                  }
                }
            }
            if (q > 0) {
                for (j in 1:q) {
                  jdx = (j - 1) * k
                  if ((t - j) > 0) {
                    tmp2 = matrix(at[(t - j), ], 1, k) \%*\% as.matrix(TH[(jdx + 
                      1):(jdx + k), ])
                    tmp = tmp - tmp2
                  }
                }
            }
            at = rbind(at, tmp)
        }
    }
    Pcnt = NULL
    idim = kp + kq
    if (inc.mean) {
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
    at = at[(ist:nT), ]
    sig = t(at) \%*\% at/(nT - pqmax)
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
    if (length(PH) > 0) 
        PH = t(PH)
    if (length(TH) > 0) 
        TH = -t(TH)
    VARMA <- list(data = da, ARorder = p, MAorder = q, cnst = include.mean, 
        coef = beta, secoef = sebeta, residuals = at, Sigma = sig, 
        aic = aic, bic = bic, Phi = PH, Theta = TH, Ph0 = Ph0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
