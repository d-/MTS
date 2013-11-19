\name{SCMid2}
\alias{SCMid2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scalar Component Model Identification
%%  ~~function to do ... ~~
}
\description{Provides detailed analysis of scalar component models 
for a specified VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMid2(zt, maxp = 2, maxq = 2, h = 0, crit = 0.05, sseq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{The T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{maxp}{Maximum AR order specified. Default is 2.
%%     ~~Describe \code{maxp} here~~
}
  \item{maxq}{Maximum MA order specified. Default is 2.
%%     ~~Describe \code{maxq} here~~
}
  \item{h}{The additional past lags used in canonical correlation analysis. 
Default is zero.
%%     ~~Describe \code{h} here~~
}
  \item{crit}{Type-I error used in testing
%%     ~~Describe \code{crit} here~~
}
  \item{sseq}{The search sequence for SCM components
%%     ~~Describe \code{sseq} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Tmatrix: The transformation matrix T. SCMorder: The orders of 
SCM.
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

\seealso{SCMid
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{zt=read.table("flourc.txt")
m2=SCMid2(zt)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (zt, maxp = 2, maxq = 2, h = 0, crit = 0.05, sseq = NULL) 
{
    if (!is.matrix(zt)) 
        zt = as.matrix(zt)
    nT = dim(zt)[1]
    k = dim(zt)[2]
    nar = maxp + 1
    nma = maxq + 1
    maxpq = maxp + maxq
    if (length(sseq) < 1) {
        sseq = c(0, 0)
        Order = matrix(0, k, 2)
        for (ell in 1:maxpq) {
            for (jj in 0:ell) {
                m = ell - jj
                j = jj
                if ((m <= maxp) && (j <= maxq)) {
                  sseq = rbind(sseq, c(m, j))
                }
            }
        }
    }
    tcases = dim(sseq)[1]
    cmx = (maxp + 1) * k * k * k
    wkspace = matrix(0, tcases, cmx)
    Nscm = rep(0, tcases)
    Tmx = NULL
    Jcnt = 0
    nc = 1
    while (Jcnt < k) {
        m = sseq[nc, 1]
        j = sseq[nc, 2]
        ist = m + j + h + 2
        Ymt = zt[ist:nT, ]
        if (m > 0) {
            for (i in 1:m) {
                Ymt = cbind(Ymt, zt[(ist - i):(nT - i), ])
            }
        }
        Ymt = as.matrix(Ymt)
        k1 = dim(Ymt)[2]
        Pt = zt[(ist - j - 1):(nT - j - 1), ]
        if (m > 0) {
            for (i in 1:m) {
                Pt = cbind(Pt, zt[(ist - 1 - i - j):(nT - 1 - 
                  i - j), ])
            }
        }
        if (h > 0) {
            for (i in 1:h) {
                Pt = cbind(Pt, zt[(ist - 1 - j - m - i):(nT - 
                  1 - j - m - i), ])
            }
        }
        Pt = as.matrix(Pt)
        mcan = cancor(Ymt, Pt)
        corsq = mcan$cor^2
        cat("For (pi,qi) = (", m, ",", j, ")", "\n")
        dsq = rep(1, k1)
        xM = as.matrix(mcan$xcoef)
        yM = as.matrix(mcan$ycoef)
        for (kk in 1:k1) {
            xM[, kk] = xM[, kk]/sqrt(sum(xM[, kk]^2))
            yM[, kk] = yM[, kk]/sqrt(sum(yM[, kk]^2))
        }
        xT = Ymt \%*\% xM
        yT = Pt \%*\% yM[, 1:k1]
        for (jj in 1:k1) {
            d1 = 1
            if (j > 0) {
                m1a = acf(xT[, jj], lag.max = j, plot = F)
                m1b = acf(yT[, jj], lag.max = j, plot = F)
                for (ij in 2:(j + 1)) {
                  d1 = d1 + 2 * m1a$acf[ij] * m1b$acf[ij]
                }
            }
            dsq[jj] = d1
        }
        chk = qnorm(1 - crit/2)
        chk1 = chk^2/nT
        idx = c(1:k1)[corsq > chk1]
        dsq[idx] = 1
        cat("Tests:", "\n")
        re = NULL
        icnt = 0
        tst = 0
        k2 = (j + 1) * k
        k2 = min(k1, k2)
        for (kk in 1:k2) {
            ik = k1 - kk + 1
            tmp = corsq[ik]/dsq[ik]
            if (tmp >= 1) 
                tmp = 0.999
            tmp = log(1 - tmp)
            tst = tst - (nT - m - j) * tmp
            deg = kk * (kk + h * k)
            pv = 1 - pchisq(tst, deg)
            re = rbind(re, c(corsq[ik], dsq[ik], tst, deg, pv))
            if (pv >= crit) 
                icnt = icnt + 1
        }
        re1 = round(re, 3)
        colnames(re1) <- c("Eigvalue", "St.dev", "Test", "deg", 
            "p-value")
        print(re1)
        cat("Summary:", "\n")
        cat("Number of SCMs detected: ", icnt, "\n")
        Nscm[nc] = icnt
        if (icnt > 0) {
            ldx = rev(c((k1 - icnt + 1):k1))
            wk1 = as.matrix(xM[, ldx])
            wkspace[nc, 1:(k1 * icnt)] = c(wk1)
            if (nc > 1) {
                wk2 = NULL
                if (Jcnt > 0) {
                  wk2 = NULL
                  for (j1 in 1:(nc - 1)) {
                    wk3 = NULL
                    ndup = min(m - sseq[j1, 1], j - sseq[j1, 
                      2])
                    m1 = sseq[j1, 1]
                    wcnt = Nscm[j1]
                    if ((ndup >= 0) && (wcnt > 0)) {
                      dcnt = (m1 + 1) * k * wcnt
                      wk3 = matrix(wkspace[j1, 1:dcnt], (m1 + 
                        1) * k, wcnt)
                      n1 = dim(wk1)[1]
                      n2 = dim(wk3)[1]
                      if (n1 > n2) {
                        wk3 = rbind(wk3, matrix(0, (n1 - n2), 
                          wcnt))
                      }
                    }
                    else {
                      mj1 = sseq[j1, 2]
                      if (((m1 + mj1) == (m + j)) && (wcnt > 
                        0)) {
                        dcnt = (m1 + 1) * k * wcnt
                        wk3 = matrix(wkspace[j1, 1:dcnt], (m1 + 
                          1) * k, wcnt)
                        n1 = dim(wk1)[1]
                        n2 = dim(wk3)[1]
                        if (n1 > n2) {
                          wk3 = rbind(wk3, matrix(0, (n1 - n2), 
                            wcnt))
                        }
                        if (n1 < n2) {
                          wk3 = wk3[1:n1, ]
                        }
                      }
                    }
                    wk2 = cbind(wk2, wk3)
                  }
                }
                if (length(wk2) < 1) {
                  Tmx = cbind(Tmx, wk1[1:k, 1:icnt])
                  cat("Found ", icnt, " new SCMs", "\n")
                  cat("Updated transformation matrix (columwise): ", 
                    "\n")
                  print(round(Tmx, 3))
                  Order[(Jcnt + 1):(Jcnt + icnt), ] = c(m, j)
                  Jcnt = Jcnt + icnt
                }
                else {
                  kold = dim(wk2)[2]
                  if (icnt <= kold) {
                    cat("No new SCM found", "\n")
                    Nscm[nc] = 0
                  }
                  if (icnt > kold) {
                    nn = min(icnt - kold, k - Jcnt)
                    cat("Found ", nn, " new SCMs", "\n")
                    xpx = t(wk2) \%*\% wk2
                    xpxinv = solve(xpx)
                    if (nn < (icnt - kold)) {
                      nover = icnt - kold - nn
                      wk1 = wk1[, 1:(icnt - nover)]
                    }
                    xpy = t(wk2) \%*\% wk1
                    beta = xpxinv \%*\% xpy
                    resi = wk1 - wk2 \%*\% beta
                    cvr = t(resi) \%*\% resi
                    mm1 = eigen(cvr)
                    vv1 = mm1$vector
                    newv = resi \%*\% vv1[, 1:nn]
                    Tmx = cbind(Tmx, newv[1:k, ])
                    cat("Transpose of Transformation-matrix: ", 
                      "\n")
                    print(round(Tmx, 3))
                    Order[(Jcnt + 1):(Jcnt + nn), ] = c(m, j)
                    Jcnt = Jcnt + nn
                    Nscm[nc] = nn
                    wkspace[nc, 1:(k1 * nn)] = c(newv)
                  }
                }
            }
            else {
                for (ij in 1:icnt) {
                  ic = k - ij + 1
                  ldx = c(1:k)[abs(xM[, ic]) == max(abs(xM[, 
                    ic]))]
                  xM[, ij] = xM[, ic]/xM[ldx[1], ic]
                }
                Tmx = cbind(Tmx, xM[, 1:icnt])
                cat("Newly detected SCMs: ", "\n")
                print(round(Tmx, 3))
                Jcnt = icnt
                Order[1:icnt, ] = c(m, j)
            }
        }
        cat("Cumulative SCMs found", Jcnt, "\n")
        nc = nc + 1
    }
    SCMid2 <- list(Tmatrix = t(Tmx), SCMorder = Order)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
