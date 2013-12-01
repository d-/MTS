\name{ccm}
\alias{ccm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-Correlation Matrices
%%  ~~function to do ... ~~
}
\description{Computes sample cross-correlation matrices of 
a multivariate time series, including simplified ccm matrix 
and p-value plot of Ljung-Box statistics.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ccm(x, lags = 12, level = FALSE, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix of vector time series, each column represents a series.
%%     ~~Describe \code{x} here~~
}
  \item{lags}{
The number of lags of CCM to be computed. Default is 12.
%%     ~~Describe \code{lags} here~~
}
  \item{level}{
A logical switch. When level=T, numerical values of CCM is printed. 
Default is no print.
%%     ~~Describe \code{level} here~~
}
  \item{output}{
A logical switch. If ouput=F, no output is given. Default is with output.
%%     ~~Describe \code{output} here~~
}
}
\details{
The p-value of Ljung-Box statistics does not include any adjustment 
in degrees of freedom.
%%  ~~ If necessary, more details than the description above ~~
}
\value{For simplified CCM
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
xt=matrix(rnorm(1500),500,3)
ccm(xt)
ccm(xt,lag=20)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
