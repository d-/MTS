\name{Kronid}
\alias{Kronid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kronecker Index Identification
%%  ~~function to do ... ~~
}
\description{Find the Kronecker indices of a k-dimensional time series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Kronid(x, plag = 5, crit = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data matrix (T-by-k) of a k-dimensional time series
%%     ~~Describe \code{x} here~~
}
  \item{plag}{The number of lags used to represent the past vector. Default 
is 5.
%%     ~~Describe \code{plag} here~~
}
  \item{crit}{Type-I error used in testing for zero canonical correlations. 
Deafult is 0.05.
%%     ~~Describe \code{crit} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{index: Kronecker indices; tests: Chi-square test statistics
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{da=read.table("flourc.txt")
Kronid(da,plag=3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
