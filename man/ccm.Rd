\name{ccm}
\alias{ccm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-Correlation Matrices
%%  ~~function to do ... ~~
}
\description{Computes sample cross-correlation matrices of 
a multivariate time series, including simplified ccm matrix 
and p-value plot of Ljung-Box statistics.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ccm(x, lags = 12, level = FALSE, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix of vector time series, each column represents a series.
%%     ~~Describe \code{x} here~~
}
  \item{lags}{
The number of lags of CCM to be computed. Default is 12.
%%     ~~Describe \code{lags} here~~
}
  \item{level}{
A logical switch. When level=T, numerical values of CCM is printed. 
Default is no print.
%%     ~~Describe \code{level} here~~
}
  \item{output}{
A logical switch. If ouput=F, no output is given. Default is with output.
%%     ~~Describe \code{output} here~~
}
}
\details{
The p-value of Ljung-Box statistics does not include any adjustment 
in degrees of freedom.
%%  ~~ If necessary, more details than the description above ~~
}
\value{For simplified CCM
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
xt=matrix(rnorm(1500),500,3)
ccm(xt)
ccm(xt,lag=20)

## The function is currently defined as
function (x, lags = 12, level = FALSE, output = T) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nT = dim(x)[1]
    k = dim(x)[2]
    if (lags < 1) 
        lags = 1
    av = apply(x, 2, mean)
    y = x
    for (j in 1:k) {
        y[, j] = x[, j] - av[j]
    }
    V1 = cov(y)
    if (output) {
        print("Covariance matrix:")
        print(V1, digits = 3)
    }
    se = sqrt(diag(V1))
    SD = diag(1/se)
    S0 = SD \%*\% V1 \%*\% SD
    ksq = k * k
    wk = matrix(0, ksq, (lags + 1))
    wk[, 1] = c(S0)
    j = 0
    if (output) {
        cat("CCM at lag: ", j, "\n")
        print(S0, digits = 3)
        cat("Simplified matrix:", "\n")
    }
    y = y \%*\% SD
    crit = 2/sqrt(nT)
    for (j in 1:lags) {
        y1 = y[1:(nT - j), ]
        y2 = y[(j + 1):nT, ]
        Sj = t(y2) \%*\% y1/nT
        Smtx = matrix(".", k, k)
        for (ii in 1:k) {
            for (jj in 1:k) {
                if (Sj[ii, jj] > crit) 
                  Smtx[ii, jj] = "+"
                if (Sj[ii, jj] < -crit) 
                  Smtx[ii, jj] = "-"
            }
        }
        if (output) {
            cat("CCM at lag: ", j, "\n")
            for (ii in 1:k) {
                cat(Smtx[ii, ], "\n")
            }
            if (level) {
                cat("Correlations:", "\n")
                print(Sj, digits = 3)
            }
        }
        wk[, (j + 1)] = c(Sj)
    }
    if (output) {
        par(mfcol = c(k, k))
        k0 = 4
        if (k > k0) 
            par(mfcol = c(k0, k0))
        tdx = c(0, 1:lags)
        jcnt = 0
        if (k > 10) {
            print("Skip the plots due to high dimension!")
        }
        else {
            for (j in 1:ksq) {
                plot(tdx, wk[j, ], type = "h", xlab = "lag", 
                  ylab = "ccf", ylim = c(-1, 1))
                abline(h = c(0))
                crit = 2/sqrt(nT)
                abline(h = c(crit), lty = 2)
                abline(h = c(-crit), lty = 2)
                jcnt = jcnt + 1
                if ((jcnt == k0^2) && (k > k0)) {
                  jcnt = 0
                  cat("Hit Enter for more plots:", "\n")
                  readline()
                }
            }
        }
        par(mfcol = c(1, 1))
        cat("Hit Enter for p-value plot of individual ccm: ", 
            "\n")
        readline()
    }
    r0i = solve(S0)
    R0 = kronecker(r0i, r0i)
    pv = rep(0, lags)
    for (i in 1:lags) {
        tmp = matrix(wk[, (i + 1)], ksq, 1)
        tmp1 = R0 \%*\% tmp
        ci = crossprod(tmp, tmp1) * nT * nT/(nT - i)
        pv[i] = 1 - pchisq(ci, ksq)
    }
    if (output) {
        plot(pv, xlab = "lag", ylab = "p-value", ylim = c(0, 
            1))
        abline(h = c(0))
        abline(h = c(0.05), col = "blue")
        title(main = "Signficance plot of CCM")
    }
    ccm <- list(ccm = wk, pvalue = pv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
