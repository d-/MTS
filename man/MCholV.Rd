\name{MCholV}
\alias{MCholV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Cholesky Volatility Model
%%  ~~function to do ... ~~
}
\description{Use Cholesky decomposition to obtain multivariate volatility models
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MCholV(rtn, size = 36, lambda = 0.96, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns.
%%     ~~Describe \code{rtn} here~~
}
  \item{size}{The initial sample size used to start recursive least squares 
estimation
%%     ~~Describe \code{size} here~~
}
  \item{lambda}{The exponential smoothing parameter. Default is 0.96.
%%     ~~Describe \code{lambda} here~~
}
  \item{p}{VAR order for the mean equation. Default is 0.
%%     ~~Describe \code{p} here~~
}
}
\details{Use recursive least squares to perform the time-varying 
Cholesky decomposition. The least squares estimates are then smoothed 
via the exponentially weighted moving-average method with decaying rate 0.96.
University GARCH(1,1) model is used for the innovations of each 
linear regression. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{Sigma.t: Volatility matrices, betat: Time-varying beta's.
bt = innovation series, and Vol: volatility series of innovations.
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{fGarch
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-ibmspko-6111.txt",header=T)
rtn=log(da[,2:4]+1)
require(fGarch)
m3=MCholV(rtn)
names(m3)
Sample means: 0.007728
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
