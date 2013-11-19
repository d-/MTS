\name{Vpmiss}
\alias{Vpmiss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial Missing of a VARMA Model
%%  ~~function to do ... ~~
}
\description{Assuming that the data is only partially missing, 
this program estimates those missing values. The model is assumed to be 
known.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Vpmiss(zt, piwgt, sigma, tmiss, mdx, cnst = NULL, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{piwgt}{pi-weights of the model in the form 
piwgt[pi0, pi1, pi2, ....]
%%     ~~Describe \code{piwgt} here~~
}
  \item{sigma}{Residual covariance matrix
%%     ~~Describe \code{sigma} here~~
}
  \item{tmiss}{Time index of the partilly missing data point
%%     ~~Describe \code{tmiss} here~~
}
  \item{mdx}{A k-dimensional indicator with "0" denoting missing 
component and ""1" denoting observed value.
%%     ~~Describe \code{mdx} here~~
}
  \item{cnst}{Constant term of the model
%%     ~~Describe \code{cnst} here~~
}
  \item{output}{values of the partially missing data
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Estimates of the missing values
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

\seealso{Vmiss
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("q-gdp-ukcaus.txt",header=T)
gdp=log(da[,3:5])
m1=VAR(gdp,1)
piwgt=m1$Phi; cnst=m1$Ph0; Sig=m1$Sigma
mdx=c(0,1,1)
m2=Vpmiss(gdp,piwgt,Sig,50,mdx,cnst)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (zt, piwgt, sigma, tmiss, mdx, cnst = NULL, output = T) 
{
    if (!is.matrix(zt)) 
        zt = as.matrix(zt)
    if (!is.matrix(piwgt)) 
        piwgt = as.matrix(piwgt)
    if (!is.matrix(sigma)) 
        sigma = as.matrix(sigma)
    k = dim(zt)[2]
    nT = dim(zt)[1]
    if (length(mdx) < 1) 
        mdx = rep(0, k)
    miss = c(1:k)[mdx == 0]
    nmiss = c(1:k)[mdx == 1]
    nm = length(miss)
    arrange = c(miss, nmiss)
    ARR = diag(rep(1, k))
    ARR = ARR[, arrange]
    ARRi = solve(ARR)
    sig = ARR \%*\% sigma \%*\% t(ARR)
    m1 = eigen(sig)
    va = sqrt(m1$values)
    P = m1$vectors
    di = diag(1/va)
    Sroot = P \%*\% di \%*\% t(P)
    if (length(cnst) < 0) {
        cnst = matrix(0, 1, k)
    }
    else {
        cnst = matrix(cnst, 1, k)
    }
    ncnst = t(ARR \%*\% t(cnst))
    lpi = dim(piwgt)[2]
    lags = floor(lpi/k)
    npiwgt = piwgt
    for (i in 1:lags) {
        icnt = (i - 1) * k
        tmp = piwgt[, (icnt + 1):(icnt + k)]
        tmp1 = ARR \%*\% tmp \%*\% ARRi
        npiwgt[, (icnt + 1):(icnt + k)] = tmp1
    }
    Y = NULL
    X = NULL
    Tpiwgt = t(npiwgt)
    zt1 = zt[, arrange]
    zt1[tmiss, 1:nm] = rep(0, nm)
    yobs = rep(0, k)
    if (k > nm) 
        yobs = matrix(zt1[tmiss, (nm + 1):k], (k - nm), 1)
    zt1[tmiss, ] = rep(0, k)
    iend = min(nT, (tmiss + lags))
    wk = matrix(0, 1, lpi)
    icnt = 0
    jend = min(tmiss - 1, lags)
    if (jend > 0) {
        for (j in 1:jend) {
            wk[1, (icnt + 1):(icnt + k)] = zt1[tmiss - j, ]
            icnt = icnt + k
        }
    }
    yt = wk \%*\% Tpiwgt + ncnst
    xt = diag(rep(1, k))
    Y1 = Sroot \%*\% t(yt)
    X1 = Sroot \%*\% xt
    if (k > nm) {
        Y = Y1 - X1[, (nm + 1):k] \%*\% yobs
    }
    else {
        Y = Y1
    }
    X = X1[, 1:nm]
    if (nm == 1) 
        X = matrix(X1[, 1], k, 1)
    Tmax = min(nT, tmiss + lags)
    iend = min(lags, Tmax - tmiss)
    if (iend > 0) {
        for (i in 1:iend) {
            wk[1, ] = c(zt1[tmiss + i - 1, ], wk[1, 1:(lpi - 
                k)])
            yt = zt1[(tmiss + i), ] - wk \%*\% Tpiwgt - ncnst
            kst = (i - 1) * k
            xt = npiwgt[, (kst + 1):(kst + k)]
            Y1 = Sroot \%*\% t(yt)
            X1 = Sroot \%*\% as.matrix(xt)
            if (k > nm) {
                Y = rbind(Y, Y1 - X1[, (nm + 1):k] \%*\% yobs)
            }
            else {
                Y = rbind(Y, Y1)
            }
            if (nm == 1) {
                X1 = matrix(X1[, 1], k, 1)
                X = rbind(X, X1)
            }
            else {
                X = rbind(X, X1[, 1:nm])
            }
        }
    }
    xpx = t(X) \%*\% X
    xpy = t(X) \%*\% Y
    xpxi = solve(xpx)
    est = xpxi \%*\% xpy
    if (output) {
        cat("Estimate of missing value at time index", tmiss, 
            "\n")
        cat("The missing idicator is: ", mdx, "\n")
        print(round(est, 5))
    }
    return(est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
