\name{VARorder}
\alias{VARorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR Order Specification
%%  ~~function to do ... ~~
}
\description{Computes information criteria and the sequential Chi-square 
statistics for a vector autoregressive process
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARorder(x, maxp = 13, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data matrix of dimension T-by-k with T being the sample size and 
k the number of time series
%%     ~~Describe \code{x} here~~
}
  \item{maxp}{The maximum VAR order entertained. Default is 13.
%%     ~~Describe \code{maxp} here~~
}
  \item{output}{A logical switch to control the output. Default is to provide output
%%     ~~Describe \code{output} here~~
}
}
\details{For a given maxp, the command computes Akaike, Bayesin and 
Hannan-Quinn informtation criteria for various VAR models using the 
data from t=maxp+1 to T. It also computes the Tiao-Box sequential Chi-square 
statistics and their p-values.
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
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{VARorderI
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{da=read.table("q-gdp-ukcaus.txt",header=T)
zt=diffM(log(da[,3:5]))
VARorder(zt)
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
    enob = nT - maxp
    y = x1[(maxp + 1):nT, ]
    ist = maxp + 1
    xmtx = cbind(rep(1, enob), x1[maxp:(nT - 1), ])
    if (maxp > 1) {
        for (i in 2:maxp) {
            xmtx = cbind(xmtx, x1[(ist - i):(nT - i), ])
        }
    }
    chidet = rep(0, (maxp + 1))
    s = cov(y) * (enob - 1)/enob
    chidet[1] = log(det(s))
    aic = rep(0, (maxp + 1))
    aic[1] = chidet[1]
    bic = aic
    hq = aic
    y = as.matrix(y)
    for (p in 1:maxp) {
        idm = k * p + 1
        xm = xmtx[, 1:idm]
        xm = as.matrix(xm)
        xpx = t(xm) \%*\% xm
        xpxinv = solve(xpx)
        xpy = t(xm) \%*\% y
        beta = xpxinv \%*\% xpy
        yhat = xm \%*\% beta
        resi = y - yhat
        sse = (t(resi) \%*\% resi)/enob
        d1 = log(det(sse))
        aic[p + 1] = d1 + (2 * p * ksq)/nT
        bic[p + 1] = d1 + (log(nT) * p * ksq)/nT
        hq[p + 1] = d1 + (2 * log(log(nT)) * p * ksq)/nT
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
    }
    if (output) {
        n1 = length(aic) - 1
        cri = cbind(c(0:n1), aic, bic, hq, c(0, Mstat), c(0, 
            pv))
        colnames(cri) <- c("p", "AIC", "BIC", "HQ", "M(p)", "p-value")
        cat("Summary table: ", "\n")
        print(round(cri, 4))
    }
    VARorder <- list(aic = aic, aicor = aicor, bic = bic, bicor = bicor, 
        hq = hq, hqor = hqor, Mstat = Mstat, Mpv = pv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
