\name{Eccm}
\alias{Eccm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extended Cross-Correlation Matrices
%%  ~~function to do ... ~~
}
\description{Compute the extended cross-correlation matrices and the 
associated two-way table of p-values of multivariate Ljung-Box 
statistics for specifying the order of a VARMA process.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Eccm(zt, maxp = 5, maxq = 6, include.mean = FALSE, rev = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{Data matrix (T-by-k) of a vector time series, where T is the sample size and 
k is the dimension. 
%%     ~~Describe \code{zt} here~~
}
  \item{maxp}{Maximum AR order entertained. Default is 5.
%%     ~~Describe \code{maxp} here~~
}
  \item{maxq}{Maximum MA order entertained. Default is 6.
%%     ~~Describe \code{maxq} here~~
}
  \item{include.mean}{A logical switch controling the mean vector in estimation. 
Default assumes zero mean.
%%     ~~Describe \code{include.mean} here~~
}
  \item{rev}{A logical switch to control the cross-correlation matrices used 
to compute the multivariate Ljung-Box statistics. Traditional way is to 
compute test statistics from lag-1 to lag-m. If rev = TRUE, then the 
test statistics are compute from lag-(m-1) to lag-m, from lag-(m-2) to lage m, etc.
%%     ~~Describe \code{rev} here~~
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
zt=matrix(rnorm(900),300,3)
m1=Eccm(zt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
