\name{VMAorder}
\alias{VMAorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Order specification of a VMA process
%%  ~~function to do ... ~~
}
\description{Performs multivariate Ljung-Box tests to specify the 
order of a VMA process
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMAorder(x, lag = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data matrix of the observed k-dimensional time series. 
Each column represents a time series.
%%     ~~Describe \code{x} here~~
}
  \item{lag}{The maximum VMA order entertained. Default is 20.
%%     ~~Describe \code{lag} here~~
}
}
\details{For a given lag, the command computes the Ljung-Box 
statistic for testing rho_j = ... = rho_lag = 0, where 
j = 1, 2, ..., lag. For a VMA(q) process, the 
Ljung-Box statistics should be significant for the first q lags, 
and insignificant thereafter.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014)
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
\examples{zt=matrix(rnorm(600),200,3)
VMAorder(zt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
