\name{mtCopula}
\alias{mtCopula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mulivariate t-Copula for Volatility
%%  ~~function to do ... ~~
}
\description{Fits a t-copula to k-dimensional standardized return series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mtCopula(rt, g1, g2, grp = NULL, th0 = NULL, m = 0, include.th0 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{A T-by-k data matrix of k standardized time series (after 
univariate volatility modeling)
%%     ~~Describe \code{rt} here~~
}
  \item{g1}{lamda1 parameter, nononegative and less than 1
%%     ~~Describe \code{g1} here~~
}
  \item{g2}{lambda2 parameter, nonnegative and satisfying lambda1+lambda2 < 1.
%%     ~~Describe \code{g2} here~~
}
  \item{grp}{a vector to indicate the number of assets divided into groups.
Default means each individual asset forms a group.
%%     ~~Describe \code{grp} here~~
}
  \item{th0}{initial estimate of theta0
%%     ~~Describe \code{th0} here~~
}
  \item{m}{number of lags used to estimate the local theta-angles
%%     ~~Describe \code{m} here~~
}
  \item{include.th0}{A logical switch to include theta0 in estimation. 
Default is to inlcude.
%%     ~~Describe \code{include.th0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{estimates: parameter estimates. Hessian: Hessian matrix.
rho.t: correlation matrices. theta.t: angel mtrices
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
da=read.table("m-ibmspko-6111.txt",header=T)
rtn=log(da[,-1]+1)
m1=dccPre(rtn,cond.dist="std")
Vol=m1$marVol; eta=m1$sresi
m2=mtCopula(eta,0.8,0.04)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
