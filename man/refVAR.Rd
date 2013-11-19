\name{refVAR}
\alias{refVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining a VAR Model
%%  ~~function to do ... ~~
}
\description{Refine a fitted VAR model by removing simultaneously insignificant parameters
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refVAR(model, fixed = NULL, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object of the command VAR or the refVAR command
%%     ~~Describe \code{model} here~~
}
  \item{fixed}{A logical matrix for VAR polynomial structure
%%     ~~Describe \code{fixed} here~~
}
  \item{thres}{Threshold used to set parameter to zero. Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{Refine a VAR fitting by setting all estimates with t-ratio 
less than the threshold (in modulus) to zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of VAR estimation
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 2)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{VAR
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("q-gdp-ukcaus.txt",header=T)
gdp=log(da[,3:5])
zt=diffM(gdp)
m1=VAR(zt,3)
m2=refVAR(m1,thres=1.0)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, fixed = NULL, thres = 1) 
{
    x = as.matrix(model$data)
    nT = dim(x)[1]
    k = dim(x)[2]
    p = model$order
    if (p < 1) 
        p = 1
    cnst = model$cnst
    fix = fixed
    if (length(fixed) == 0) {
        coef = as.matrix(model$coef)
        secoef = as.matrix(model$secoef)
        nr = dim(coef)[1]
        nc = dim(coef)[2]
        for (i in 1:nr) {
            for (j in 1:nc) {
                if (secoef[i, j] < 10^(-8)) 
                  secoef[i, j] = 1
            }
        }
        fix = matrix(1, nr, k)
        xmtx = NULL
        ist = p + 1
        y = x[ist:nT, ]
        ne = nT - p
        if (cnst) 
            xmtx = matrix(1, ne, 1)
        for (j in 1:p) {
            xmtx = cbind(xmtx, x[(ist - j):(nT - j), ])
        }
        xmtx = as.matrix(xmtx)
        for (j in 1:k) {
            tt = abs(coef[, j]/secoef[, j])
            idx = c(1:nr)[tt == min(tt)]
            idx1 = idx[1]
            if (tt[idx1] < thres) 
                fix[idx, j] = 0
        }
        for (j in 1:k) {
            npar = sum(fix[, j])
            while (npar > 0) {
                jdx = c(1:nr)[fix[, j] == 1]
                xp = as.matrix(xmtx[, jdx])
                nxp = dim(xp)[2]
                m1 = lm(y[, j] ~ -1 + xp)
                m2 = summary(m1)
                est = m1$coefficients
                se1 = sqrt(diag(m2$cov.unscaled)) * m2$sigma
                tt = abs(est/se1)
                idx = c(1:nxp)[tt == min(tt)]
                idx1 = idx[1]
                if (tt[idx1] < thres) {
                  fix[jdx[idx], j] = 0
                  npar = sum(fix[, j])
                }
                else {
                  npar = 0
                }
            }
        }
    }
    mm = VAR(x, p, output = T, include.mean = cnst, fixed = fix)
    refVAR <- list(data = mm$data, order = p, cnst = cnst, coef = mm$coef, 
        aic = mm$aic, bic = mm$bic, hq = mm$hq, residuals = mm$residuals, 
        secoef = mm$secoef, Sigma = mm$Sigma, Phi = mm$Phi, Ph0 = mm$Ph0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
