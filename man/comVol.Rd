\name{comVol}
\alias{comVol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Common Volatility
%%  ~~function to do ... ~~
}
\description{Compute the principal volatility components after a VAR(p) 
fitting
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
comVol(rtn, m = 10, p = 1, stand = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns
%%     ~~Describe \code{rtn} here~~
}
  \item{m}{The number of lags used to compute generalized cross-Kurtosis 
matrix
%%     ~~Describe \code{m} here~~
}
  \item{p}{VAR order for the mean equation
%%     ~~Describe \code{p} here~~
}
  \item{stand}{A logical switch to standardize the returns
%%     ~~Describe \code{stand} here~~
}
}
\details{Perform a VAR(p) fits, if any. Then, use the residual series to 
perform principal volatility component analysis. The ARCH test statistics 
are also computed for the sample principal components
%%  ~~ If necessary, more details than the description above ~~
}
\value{residuals: The residuals of VAR fit, which are are data used in 
principal volatility component analysis. Values: eigenvalues. Vectors: 
eigenvectors. M: the transformation matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 7)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay and Y.B. Hu
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
da=read.table("w-7fx.txt",header=T)
fx=log(da[,4:10])
rtn=diffM(fx)
m1=comVol(rtn,p=5)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
