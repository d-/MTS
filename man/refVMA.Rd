\name{refVMA}
\alias{refVMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining VMA Models
%%  ~~function to do ... ~~
}
\description{Refines a fitted VMA model by setting insignificant parameters to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refVMA(model, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object from the command VMA or the refVMA command
%%     ~~Describe \code{model} here~~
}
  \item{thres}{A threshold for individual t-ratio of parameter estimate. 
Default is 2.
%%     ~~Describe \code{thres} here~~
}
}
\details{The program simultaneously sets all estimates with t-ratios 
less than the threshold (in modulus) to zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of VMA estimation
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

\seealso{VMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
theta=matrix(c(0.4,0.5,0.0,0.6),2,2)
sig=matrix(c(3,1,1,1),2,2)
m1=VARMAsim(400,malags=c(1),theta=theta,sigma=sig)
zt=m1$series
m2=VMA(zt,q=1,include.mean=F)
m3=refVMA(m2,thres=1.8)
names(m3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, thres = 2) 
{
    x = model$data
    q = model$MAorder
    cnst = model$cnst
    coef = as.matrix(model$coef)
    secoef = as.matrix(model$secoef)
    nr = dim(coef)[1]
    nc = dim(coef)[2]
    for (i in 1:nc) {
        idx = is.na(secoef[, i])
        jdx = c(1:nr)[idx == T]
        secoef[jdx, i] = 0.01
    }
    fix = matrix(0, nr, nc)
    for (j in 1:nc) {
        tt = coef[, j]/secoef[, j]
        idx = c(1:nr)[abs(tt) >= thres]
        fix[idx, j] = 1
    }
    if (cnst) {
        tt = coef[1, ]/secoef[1, ]
        idx = c(1:nc)[abs(tt) > 1]
        if (length(idx) > 0) 
            fix[1, idx] = 1
    }
    mm = VMA(x, q = q, include.mean = cnst, fixed = fix)
    refVMA <- list(data = x, MAorder = q, cnst = cnst, coef = mm$coef, 
        secoef = mm$secoef, residuals = mm$residuals, Sigma = mm$Sigma, 
        aic = mm$aic, bic = mm$bic, mu = mm$mu, Theta = mm$Theta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
