\name{VARorderI}
\alias{VARorderI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR order specification I
%%  ~~function to do ... ~~
}
\description{This program is similar to VARorder, but it uses 
observations from t=p+1 to T to compute the information criteria for 
a given VAR(p) model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARorderI(x, maxp = 13, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of vector time series
%%     ~~Describe \code{x} here~~
}
  \item{maxp}{The maximum VAR order entertained
%%     ~~Describe \code{maxp} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{For a given VAR(p) model, the program uses observations 
from t=p+1 to T to compute the information criteria.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2013)
%% ~put references to the literature/web site here ~
}
\author{Ruey S Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{VARorder
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, maxp = 13, output = T) 
{
    x1 = as.matrix(x)
    nT = nrow(x1)
    k = ncol(x1)
    ksq = k * k
    maxp = max(maxp, 2 * floor(log(nT)))
    chidet = rep(0, (maxp + 1))
    s = cov(x1) * (nT - 1)/nT
    chidet[1] = log(det(s))
    aic = chidet
    bic = aic
    hq = aic
    for (p in 1:maxp) {
        idm = k * p + 1
        ist = p + 1
        enob = nT - p
        y = as.matrix(x1[ist:nT, ])
        xmtx = rep(1, enob)
        for (j in 1:p) {
            xmtx = cbind(xmtx, x1[(ist - j):(nT - j), ])
        }
        xm = as.matrix(xmtx)
        xpx = t(xm) \%*\% xm
        xpxinv = solve(xpx)
        xpy = t(xm) \%*\% y
        beta = xpxinv \%*\% xpy
        yhat = xm \%*\% beta
        resi = y - yhat
        sse = (t(resi) \%*\% resi)/enob
        d1 = log(det(sse))
        aic[p + 1] = d1 + (2 * p * ksq)/enob
        bic[p + 1] = d1 + (log(enob) * p * ksq)/enob
        hq[p + 1] = d1 + (2 * log(log(enob)) * p * ksq)/enob
        chidet[p + 1] = d1
    }
    maic = min(aic)
    aicor = c(1:(maxp + 1))[aic == maic] - 1
    mbic = min(bic)
    bicor = c(1:(maxp + 1))[bic == mbic] - 1
    mhq = min(hq)
    hqor = c(1:(maxp + 1))[hq == mhq] - 1
    Mstat = rep(0, maxp)
    pv = rep(0, maxp)
    for (j in 1:maxp) {
        Mstat[j] = (nT - maxp - k * j - 1.5) * (chidet[j] - chidet[j + 
            1])
        pv[j] = 1 - pchisq(Mstat[j], ksq)
    }
    if (output) {
        cat("selected order: aic = ", aicor, "\n")
        cat("selected order: bic = ", bicor, "\n")
        cat("selected order: hq = ", hqor, "\n")
        cat("M statistic and its p-value", "\n")
        tmp = cbind(Mstat, pv)
        print(tmp, digits = 4)
        n1 = length(aic) - 1
        cri = cbind(c(0:n1), aic, bic, hq, c(0, Mstat), c(0, 
            pv))
        colnames(cri) <- c("p", "AIC", "BIC", "HQ", "M(p)", "p-value")
        cat("Summary table: ", "\n")
        print(cri, digits = 5)
    }
    VARorderI <- list(aic = aic, aicor = aicor, bic = bic, bicor = bicor, 
        hq = hq, hqor = hqor, Mstat = Mstat, Mpv = pv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
