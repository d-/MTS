\name{mtCopulaVol}
\alias{mtCopulaVol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Volatility Matrices of Multivariate t-Copula Model
%%  ~~function to do ... ~~
}
\description{Compute the volatility matrices of a fitted mtCopula model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mtCopulaVol(par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{paramer vector
%%     ~~Describe \code{par} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 7)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{An internal function, not intended as a general command
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (par) 
{
    if (mtCopulaincludeth0) {
        th1 = par[4:length(par)]
    }
    else {
        th1 <- mtCopulath0
    }
    par0 = 1 - par[2] - par[3]
    tmp0 = par0 * th1
    Theta <- mtCopulaTheta
    ncor <- mtCopulancor
    grp <- mtCopulagrp
    m <- mtCopulam
    ist = m + 1
    nT = dim(Theta)[1]
    Xmtx = the2Xmtx(th1, grp)
    Rt = t(Xmtx) \%*\% Xmtx
    rho.t = NULL
    TH = NULL
    for (i in 1:m) {
        rho.t = rbind(rho.t, c(Rt))
        TH = rbind(TH, th1)
    }
    for (t in ist:nT) {
        thet = Theta[t, ]
        tht = tmp0 + TH[t - 1, ] * par[2] + par[3] * thet
        TH = rbind(TH, tht)
        Xmtx = the2Xmtx(tht, grp)
        Rt = t(Xmtx) \%*\% Xmtx
        rho.t = rbind(rho.t, c(Rt))
    }
    mtCopulaVol <- list(rho.t = rho.t, angles = TH)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
