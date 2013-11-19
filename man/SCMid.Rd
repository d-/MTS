\name{SCMid}
\alias{SCMid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scalar Component Identification
%%  ~~function to do ... ~~
}
\description{Find the overall order of a VARMA process via the 
scalar component model approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMid(zt, maxp = 5, maxq = 5, h = 0, crit = 0.05, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{The T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{maxp}{Maximum AR order entertained. Default is 5.
%%     ~~Describe \code{maxp} here~~
}
  \item{maxq}{Maximum MA order entertained. Default is 5.
%%     ~~Describe \code{maxq} here~~
}
  \item{h}{The additional past lags used in canonical correlation analysis.
Default is 0.
%%     ~~Describe \code{h} here~~
}
  \item{crit}{Type-I error of the test
%%     ~~Describe \code{crit} here~~
}
  \item{output}{A logical switch to control the output.
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Nmtx: The table of the numbers of zero canonical correlations, and 
DDmtx: 
the table of diagonal differences.
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
m1=SCMid(zt)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (zt, maxp = 5, maxq = 5, h = 0, crit = 0.05, output = FALSE) 
{
    if (!is.matrix(zt)) 
        zt = as.matrix(zt)
    nT = dim(zt)[1]
    k = dim(zt)[2]
    nar = maxp + 1
    nma = maxq + 1
    zeroTbl = matrix(0, nar, nma)
    diagDif = matrix(0, nar, nma)
    for (m in 0:maxp) {
        for (j in 0:maxq) {
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
            m1 = cancor(Ymt, Pt)
            corsq = m1$cor^2
            if (output) {
                cat("For (m,j) = ", c(m, j), "\n")
                cat("Squares of canonical correlations: ", "\n")
                print(round(sort(corsq), 3))
            }
            dsq = rep(1, k1)
            if (j > 0) {
                xM = as.matrix(m1$xcoef)
                yM = as.matrix(m1$ycoef)
                for (kk in 1:k1) {
                  xM[, kk] = xM[, kk]/sqrt(sum(xM[, kk]^2))
                  yM[, kk] = yM[, kk]/sqrt(sum(yM[, kk]^2))
                }
                xT = Ymt \%*\% xM
                yT = Pt \%*\% yM[, 1:k1]
                for (jj in 1:k1) {
                  d1 = 1
                  m1a = acf(xT[, jj], lag.max = j, plot = F)
                  m1b = acf(yT[, jj], lag.max = j, plot = F)
                  for (ij in 2:(j + 1)) {
                    d1 = d1 + 2 * m1a$acf[ij] * m1b$acf[ij]
                  }
                  dsq[jj] = d1
                }
            }
            chk = qnorm(1 - crit/2)
            chk1 = chk^2/nT
            idx = c(1:k1)[corsq > chk1]
            dsq[idx] = 1
            if (output) {
                cat("Variance of canonical correlations", "\n")
                print(round(dsq, 3))
                cat("Test results: ", "\n")
            }
            icnt = 0
            tst = 0
            k2 = (j + 1) * k
            k2 = min(k1, k2)
            for (kk in 1:k2) {
                tmp = corsq[k1 - kk + 1]/dsq[k1 - kk + 1]
                if (tmp >= 1) 
                  tmp = 0.999
                tmp = log(1 - tmp)
                tst = tst - (nT - m - j) * tmp
                deg = kk * (kk + h * k)
                pv = 1 - pchisq(tst, deg)
                if (output) 
                  print(c(tst, deg, pv), digits = 3)
                if (pv >= crit) 
                  icnt = icnt + 1
            }
            zeroTbl[(m + 1), (j + 1)] = icnt
        }
    }
    cat("Column: MA order", "\n")
    cat("Row   : AR order", "\n")
    cat("Number of zero canonical correlations", "\n")
    colnames(zeroTbl) <- c(c(0:maxq))
    rownames(zeroTbl) <- c(c(0:maxp))
    printCoefmat(zeroTbl)
    diagDif = zeroTbl
    for (i in 1:maxp) {
        for (j in 1:maxq) {
            diagDif[(i + 1), (j + 1)] = min(zeroTbl[(i + 1), 
                (j + 1)] - zeroTbl[i, j], k)
        }
    }
    cat("Diagonal Differences: ", "\n")
    colnames(diagDif) <- c(c(0:maxq))
    rownames(diagDif) <- c(c(0:maxp))
    printCoefmat(diagDif)
    SCMid <- list(Nmtx = zeroTbl, DDmtx = diagDif)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
