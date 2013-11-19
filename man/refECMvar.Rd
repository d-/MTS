\name{refECMvar}
\alias{refECMvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining Error-Correction Model for VAR series
%%  ~~function to do ... ~~
}
\description{Refining an estimated ECM VAR model by setting insignificant 
estimates to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refECMvar(m1, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An object of the ECMvar command or the refECMvar command
%%     ~~Describe \code{m1} here~~
}
  \item{thres}{Threshold for individual t-ratio. The default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{Set simultaneously all estimates with t-ratio less than 
the threshold to zero (in modulus).
%%  ~~ If necessary, more details than the description above ~~
}
\value{Constrained estimation results of a ECM VAR model
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m2=refECMvar(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m1, thres = 1) 
{
    x = m1$data
    m = m1$ncoint
    p = m1$arorder
    include.con = m1$include.const
    alpha = m1$alpha
    se.alpha = m1$se.alpha
    beta = m1$beta
    se.beta = m1$se.beta
    Phip = m1$Phip
    se.Phip = m1$se.Phip
    mm = ECMvar(x, p, beta, include.const = include.con, fixed = NULL)
    alpha = mm$alpha
    se.alpha = mm$se.alpha
    beta = mm$beta
    se.beta = mm$se.beta
    Phip = mm$Phip
    se.Phip = mm$se.Phip
    aic = mm$aic
    bic = mm$bic
    resi = mm$residuals
    sse = mm$Sigma
    refECMvar <- list(data = x, ncoint = m, arorder = p, include.const = include.const, 
        alpha = alpha, se.alpha = se.alpha, beta = beta, se.beta = se.beta, 
        aic = aic, bic = bic, residuals = resi, Phip = phip, 
        se.Phip = se.phip, Sigma = sse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
