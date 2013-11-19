\name{refECMvar1}
\alias{refECMvar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining ECM for a VAR process
%%  ~~function to do ... ~~
}
\description{Performs constrained least squares estimation of 
a ECM VAR model with known co-integrated processes
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refECMvar1(m1, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An object of the ECMvar1 command or the refECMvar1 command
%%     ~~Describe \code{m1} here~~
}
  \item{thres}{Threshold for individual t-ratio. Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Constrained estimation results of an ECM VAR model
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 5)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ECMvar1, refECMvar
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m2=refECMvar1(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m1, thres = 1) 
{
    x = m1$data
    wt = m1$wt
    p = m1$arorder
    include.con = m1$include.const
    coef = m1$coef
    secoef = m1$secoef
    idm = dim(coef)[1]
    k = dim(coef)[2]
    fix = matrix(0, idm, k)
    for (i in 1:k) {
        tra = coef[, i]/secoef[, i]
        idx = c(1:idm)[abs(tra) > thres]
        fix[idx, i] = 1
    }
    mm = ECMvar1(x, p, wt, include.const = include.con, fixed = fix)
    beta = mm$coef
    sdbeta = mm$secoef
    aic = mm$aic
    bic = mm$bic
    resi = mm$residuals
    sse = mm$Sigma
    refECMvar1 <- list(data = x, wt = wt, arorder = p, include.const = include.con, 
        coef = beta, aic = aic, bic = bic, residuals = resi, 
        secoef = sdbeta, Sigma = sse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line