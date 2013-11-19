\name{refsVARMA}
\alias{refsVARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining a Seasonal VARMA Model
%%  ~~function to do ... ~~
}
\description{Refines a fitted seasonal VARMA model by setting insignificant 
estimates to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refsVARMA(model, thres = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object of the sVARMA command or the refsVARMA 
command
%%     ~~Describe \code{model} here~~
}
  \item{thres}{Threshold for individual t-ratio. Default is 0.8. 

%%     ~~Describe \code{thres} here~~
}
}
\details{The command removes simultaneously all parameters with t-ratio 
less than the threshold in modulus.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a sesonal VARMA fitting
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

\seealso{sVARMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m4=refsVARMA(m3,thres=1.2)
names(m4)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, thres = 0.8) 
{
    x = model$data
    order = model$order
    sorder = model$sorder
    s = model$period
    cnst = model$cnst
    swi = model$switch
    coef = as.matrix(model$coef)
    secoef = as.matrix(model$secoef)
    nr = dim(coef)[1]
    nc = dim(coef)[2]
    for (j in 1:nc) {
        idx = is.na(secoef[, j])
        jdx = c(1:nr)[idx == T]
        secoef[jdx, j] = 0.01
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
    mm = sVARMA(x, order, sorder, s, include.mean = cnst, fixed = fix, 
        switch = swi)
    refsVARMA <- list(data = x, coef = mm$coef, secoef = mm$secoef, 
        order = mm$order, sorder = mm$sorder, period = mm$period, 
        cnst = cnst, residuals = mm$residuals, regPhi = mm$regPhi, 
        seaPhi = mm$seaPhi, regTheta = mm$regTheta, seaTheta = mm$seaTheta, 
        Ph0 = mm$Pho, Sigma = mm$Sigma, aic = mm$aic, bic = mm$bic, 
        switch = mm$switch)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
