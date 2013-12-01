\name{iniSCM}
\alias{iniSCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial Estimates of SCM Approach
%%  ~~function to do ... ~~
}
\description{Obtain initial parameter estimates of a model specified by the 
SCM approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
iniSCM(da, at, locTmtx, locAR, locMA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{at}{A T-by-k initial estimates of the innovations
%%     ~~Describe \code{at} here~~
}
  \item{locTmtx}{Locator of 1 for the transformation matrix
%%     ~~Describe \code{locTmtx} here~~
}
  \item{locAR}{Locator of AR polynomial structure
%%     ~~Describe \code{locAR} here~~
}
  \item{locMA}{Locator of MA polynomial structure
%%     ~~Describe \code{locMA} here~~
}
}
\details{Initial function used by SCMfit, not intended as a general 
command
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 4)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{SCMfit
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
