\name{VARXorder}
\alias{VARXorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VARX Order Specification
%%  ~~function to do ... ~~
}
\description{Specifies the orders of a VARX model, including AR order 
and the number of lags of exogenous variables
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARXorder(x, exog, maxp = 13, maxm = 3, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{x} here~~
}
  \item{exog}{A T-by-m data matrix of exogenous variables
%%     ~~Describe \code{exog} here~~
}
  \item{maxp}{The maximum VAR order detertained
%%     ~~Describe \code{maxp} here~~
}
  \item{maxm}{The maximum lags of exogenous variables entertained
%%     ~~Describe \code{maxm} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{Computes the information criteria of a VARX process
%%  ~~ If necessary, more details than the description above ~~
}
\value{Information criteria
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
da=read.table("m-gasoil.txt",header=T)
zt=da[,3:4]; xt=da[,5:6]
m1=VARXorder(zt,xt,maxp=11,maxm=2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, exog, maxp = 13, maxm = 3, output = T) 
{
    x1 = as.matrix(x)
    exog = as.matrix(exog)
    nT = dim(x1)[1]
    k = dim(x1)[2]
    ksq = k * k
    maxp = max(maxp, 2 * floor(log(nT)))
    nT1 = dim(exog)[1]
    m = dim(exog)[2]
    if (nT1 > nT) {
        cat("Different nobs:", c(nT, nT1), "\n")
    }
    if (nT > nT1) {
        cat("Different nobs:", c(nT, nT1), "\n")
        nT = nT1
    }
    aic = matrix(0, maxp + 1, maxm + 1)
    bic = aic
    hq = aic
    for (mm in 0:maxm) {
        isto = mm + 1
        y = x1[isto:nT, ]
        xm = rep(1, rep(nT - mm))
        for (i in 0:mm) {
            xm = cbind(xm, exog[(isto - i):(nT - i), ])
        }
        xm = as.matrix(xm)
        xpx = t(xm) \%*\% xm
        xpxi = solve(xpx)
        xpy = t(xm) \%*\% y
        beta = xpxi \%*\% xpy
        y = y - xm \%*\% beta
        s = t(y) \%*\% y/(nT - mm)
        enob = nT - mm
        c1 = log(det(s))
        aic[1, mm + 1] = c1 + 2 * k * mm/enob
        bic[1, mm + 1] = c1 + log(enob) * k * m/enob
        hq[1, mm + 1] = c1 + 2 * log(log(enob)) * k * m/enob
        for (p in 1:maxp) {
            ist = max(mm + 1, p + 1)
            enob = nT - ist + 1
            y = as.matrix(x1[ist:nT, ])
            xmtx = rep(1, enob)
            for (i in 0:mm) {
                xmtx = cbind(xmtx, exog[(ist - i):(nT - i), ])
            }
            for (j in 1:p) {
                xmtx = cbind(xmtx, x1[(ist - j):(nT - j), ])
            }
            xm1 = as.matrix(xmtx)
            xpx = t(xm1) \%*\% xm1
            xpxinv = solve(xpx)
            xpy = t(xm1) \%*\% y
            beta = xpxinv \%*\% xpy
            resi = y - xm1 \%*\% beta
            sse = (t(resi) \%*\% resi)/enob
            d1 = log(det(sse))
            npar = p * ksq + k * (mm + 1)
            aic[p + 1, mm + 1] = d1 + (2 * npar)/enob
            bic[p + 1, mm + 1] = d1 + (log(enob) * npar)/enob
            hq[p + 1, mm + 1] = d1 + (2 * log(log(enob)) * npar)/enob
        }
    }
    aicor = ind.min(aic) - 1
    bicor = ind.min(bic) - 1
    hqor = ind.min(hq) - 1
    if (output) {
        cat("selected order(p,s): aic = ", aicor, "\n")
        cat("selected order(p,s): bic = ", bicor, "\n")
        cat("selected order(p,s): hq = ", hqor, "\n")
    }
    VARXorder <- list(aic = aic, aicor = aicor, bic = bic, bicor = bicor, 
        hq = hq, hqor = hqor)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
