\name{Eccm}
\alias{Eccm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extended Cross-Correlation Matrices
%%  ~~function to do ... ~~
}
\description{Compute the extended cross-correlation matrices and the 
associated two-way table of p-values of multivariate Ljung-Box 
statistics for specifying the order of a VARMA process.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Eccm(zt, maxp = 5, maxq = 6, include.mean = FALSE, rev = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{Data matrix (T-by-k) of a vector time series, where T is the sample size and 
k is the dimension. 
%%     ~~Describe \code{zt} here~~
}
  \item{maxp}{Maximum AR order entertained. Default is 5.
%%     ~~Describe \code{maxp} here~~
}
  \item{maxq}{Maximum MA order entertained. Default is 6.
%%     ~~Describe \code{maxq} here~~
}
  \item{include.mean}{A logical switch controling the mean vector in estimation. 
Default assumes zero mean.
%%     ~~Describe \code{include.mean} here~~
}
  \item{rev}{A logical switch to control the cross-correlation matrices used 
to compute the multivariate Ljung-Box statistics. Traditional way is to 
compute test statistics from lag-1 to lag-m. If rev = TRUE, then the 
test statistics are compute from lag-(m-1) to lag-m, from lag-(m-2) to lage m, etc.
%%     ~~Describe \code{rev} here~~
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
zt=matrix(rnorm(900),300,3)
m1=Eccm(zt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (zt, maxp = 5, maxq = 6, include.mean = FALSE, rev = TRUE) 
{
    if (!is.matrix(zt)) 
        zt = as.matrix(zt)
    x = zt
    nT = dim(x)[1]
    k = dim(x)[2]
    if (include.mean) {
        av = apply(x, 2, mean)
        for (i in 1:k) {
            x[, i] = x[, i] - av[i]
        }
    }
    ksq = k * k
    if (rev) {
        maxq1 = maxq + 1
        m1 = revmq(x, maxq1, output = F)
        vEccm = m1$ccm[, 1:maxq1]
        pEccm = m1$pvalue
        ARcoef = NULL
        for (p in 1:maxp) {
            Phi = NULL
            m1 = VAR(x, p, include.mean = F, output = F)
            Phi = rbind(Phi, m1$Phi)
            resi = m1$residuals
            m2 = revmq(resi, maxq1, output = F)
            pv1 = m2$pvalue[1]
            Eccmit = m2$ccm[, 1]
            y = x[(p + 1):nT, ]
            xreg = NULL
            for (j in 1:p) {
                xreg = cbind(xreg, x[(p + 1 - j):(nT - j), ])
            }
            kx = dim(xreg)[2]
            for (it in 1:maxq) {
                y = y[-1, ]
                yT = dim(y)[1]
                if (it == 1) {
                  xreg = xreg[-1, ]
                }
                else {
                  nx = dim(xreg)[2]
                  xT = dim(xreg)[1]
                  xreg = cbind(xreg[-1, 1:kx], xreg[-xT, (kx + 
                    1):nx])
                }
                TT = dim(resi)[1]
                xreg = cbind(xreg, resi[-TT, ])
                xpx = crossprod(xreg, xreg)/nT
                xpy = crossprod(xreg, y)/nT
                xpxinv = solve(xpx)
                beta = xpxinv \%*\% xpy
                wt = y - xreg[, 1:kx] \%*\% beta[1:kx, ]
                resi = y - xreg \%*\% beta
                Phi = rbind(Phi, t(beta[1:kx, ]))
                m3 = revmq(wt, maxq1, output = F)
                Eccmit = cbind(Eccmit, m3$ccm[, (it + 1)])
                pv1 = c(pv1, m3$pvalue[it + 1])
            }
            ARcoef = cbind(ARcoef, Phi)
            vEccm = rbind(vEccm, Eccmit)
            pEccm = rbind(pEccm, pv1)
        }
    }
    else {
        m1 = ccm(x, (maxq + 1), output = F)
        vEccm = m1$ccm[, 2:(maxq + 2)]
        pEccm = m1$pvalue
        ARcoef = NULL
        for (p in 1:maxp) {
            Phi = NULL
            m1 = VAR(x, p, include.mean = F, output = F)
            Phi = rbind(Phi, m1$Phi)
            resi = m1$residuals
            m2 = ccm(resi, 1, output = F)
            pv1 = m2$pvalue
            Eccmit = matrix(m2$ccm[, 2], ksq, 1)
            y = x[(p + 1):nT, ]
            xreg = NULL
            for (j in 1:p) {
                xreg = cbind(xreg, x[(p + 1 - j):(nT - j), ])
            }
            kx = dim(xreg)[2]
            for (it in 1:maxq) {
                y = y[-1, ]
                yT = dim(y)[1]
                if (it == 1) {
                  xreg = xreg[-1, ]
                }
                else {
                  nx = dim(xreg)[2]
                  xT = dim(xreg)[1]
                  xreg = cbind(xreg[-1, 1:kx], xreg[-xT, (kx + 
                    1):nx])
                }
                TT = dim(resi)[1]
                xreg = cbind(xreg, resi[-TT, ])
                xpx = crossprod(xreg, xreg)/nT
                xpy = crossprod(xreg, y)/nT
                xpxinv = solve(xpx)
                beta = xpxinv \%*\% xpy
                wt = y - xreg[, 1:kx] \%*\% beta[1:kx, ]
                resi = y - xreg \%*\% beta
                Phi = rbind(Phi, t(beta[1:kx, ]))
                m3 = ccm(wt, it + 1, output = F)
                Eccmit = cbind(Eccmit, m3$ccm[, (it + 2)])
                pv1 = c(pv1, m3$pvalue[it + 1])
            }
            ARcoef = cbind(ARcoef, Phi)
            vEccm = rbind(vEccm, Eccmit)
            pEccm = rbind(pEccm, pv1)
        }
    }
    cat("p-values table of Extended Cross-correlation Matrices:", 
        "\n")
    cat("Column: MA order", "\n")
    cat("Row   : AR order", "\n")
    colnames(pEccm) <- c(c(0:maxq))
    rownames(pEccm) <- c(c(0:maxp))
    tmp = round(pEccm, 4)
    printCoefmat(tmp)
    Eccm <- list(pEccm = pEccm, vEccm = vEccm, ARcoef = ARcoef)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
