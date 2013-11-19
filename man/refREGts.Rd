\name{refREGts}
\alias{refREGts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining a Regression Model with Time Series Errors
%%  ~~function to do ... ~~
}
\description{Refines a fitted REGts by setting simultaneously parameters 
with t-ratios less than the threshold (in modulus) to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refREGts(m1, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An output object from the REGts command or refREGts command
%%     ~~Describe \code{m1} here~~
}
  \item{thres}{Threshold value for individual t-ratio. Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a REGts model
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

\seealso{refVAR, refVARMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-gasoil.txt",header=T)
zt=da[,3:4]; xt=da[,5:6]
m1=REGts(zt,2,xt)
m2=refREGts(1m,thres=1.2)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m1, thres = 1) 
{
    zt = m1$data
    xt = m1$xt
    p = m1$aror
    include.m = m1$include.mean
    coef = m1$coef
    se.coef = m1$se.coef
    k = dim(zt)[2]
    n1 = dim(coef)[1]
    kx = dim(xt)[2]
    if (include.m) 
        kx = kx + 1
    icnt = 0
    par = NULL
    separ = NULL
    fix1 = matrix(0, kx, k)
    for (i in 1:k) {
        tt = coef[1:kx, i]/se.coef[1:kx, i]
        idx = c(1:kx)[abs(tt) > thres]
        ii = length(idx)
        if (ii > 0) {
            par = c(par, coef[idx, i])
            separ = c(separ, se.coef[idx, i])
            fix1[idx, i] = 1
            icnt = icnt + ii
        }
    }
    r2 = n1 - kx
    Phi = coef[(kx + 1):n1, ]
    sePhi = se.coef[(kx + 1):n1, ]
    fix2 = matrix(0, r2, k)
    for (i in 1:k) {
        tt = Phi[, i]/sePhi[, i]
        jdx = c(1:r2)[abs(tt) > thres]
        jj = length(jdx)
        if (jj > 0) {
            par = c(par, Phi[jdx, i])
            separ = c(separ, sePhi[jdx, i])
            icnt = icnt + jj
            fix2[jdx, i] = 1
        }
    }
    fix = rbind(fix1, fix2)
    mm = REGts(zt, p, xt, include.mean = include.m, fixed = fix, 
        par = par, se.par = separ)
    coef = mm$coef
    se.coef = mm$se.coef
    Phi = t(mm$Phi)
    se.Phi = t(mm$se.Phi)
    beta = t(mm$beta)
    se.beta = t(mm$se.beta)
    resi = mm$residuals
    sig = mm$Sigma
    refREGts <- list(data = zt, xt = xt, aror = p, include.mean = include.m, 
        Phi = t(Phi), se.Phi = t(se.Phi), beta = t(beta), se.beta = t(se.beta), 
        residuals = resi, Sigma = sig, coef = coef, se.coef = se.coef)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
