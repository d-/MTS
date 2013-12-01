\name{EWMAvol}
\alias{EWMAvol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponentially Weighted Moving-Average for Volatility
%%  ~~function to do ... ~~
}
\description{Use exponentially weighted moving-average method to 
compute the volatility matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EWMAvol(rtn, lambda = 0.96)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns, assuming the 
mean is zero
%%     ~~Describe \code{rtn} here~~
}
  \item{lambda}{Smoothing parameter. The deafult is 0.96. If lambda is negative, then multivariate Gaussian likelihood is used to estimate the smoothing parameter.
%%     ~~Describe \code{lambda} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Volatility matrices
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-dec125910-6111.txt",header=T)
rtn=log(da[,2:4]+1)
m1=VAR(rtn,1)
at=m1$residuals
m2=EWMAvol(at,lambda=-0.1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
