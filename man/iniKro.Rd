\name{iniKro}
\alias{iniKro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial Estimates for Hronecker Index Approach
%%  ~~function to do ... ~~
}
\description{Obtain initial parameter estimates of a VARMA model specified 
by the Kronecker index approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
iniKro(da, at, ARid, MAid, include.mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{at}{A T-by-k preliminary estimates of the innovations
%%     ~~Describe \code{at} here~~
}
  \item{ARid}{AR polynomial structure
%%     ~~Describe \code{ARid} here~~
}
  \item{MAid}{MA polynomial structure
%%     ~~Describe \code{MAid} here~~
}
  \item{include.mean}{A logical switch to include the mean vector
%%     ~~Describe \code{include.mean} here~~
}
}
\details{Initial function, not intended as a general command
%%  ~~ If necessary, more details than the description above ~~
}
\value{Initial parameter estimates
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

\seealso{Kronfit
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
