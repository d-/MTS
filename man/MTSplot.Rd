\name{MTSplot}
\alias{MTSplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Time Series Plot
%%  ~~function to do ... ~~
}
\description{Provides time plots of a vector time series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MTSplot(data, caltime = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Vector time series data matrix
%%     ~~Describe \code{data} here~~
}
  \item{caltime}{Calendar time. Default is NULL, that is, using 
time index
%%     ~~Describe \code{caltime} here~~
}
}
\details{Provides time plots of a vector time series. The output frame 
depends on the dimension of the time series
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
xt=matrix(rnorm(1500),500,3)
MTSplot(xt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
