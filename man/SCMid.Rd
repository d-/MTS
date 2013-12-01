\name{SCMid}
\alias{SCMid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scalar Component Identification
%%  ~~function to do ... ~~
}
\description{Find the overall order of a VARMA process via the 
scalar component model approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMid(zt, maxp = 5, maxq = 5, h = 0, crit = 0.05, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{The T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{maxp}{Maximum AR order entertained. Default is 5.
%%     ~~Describe \code{maxp} here~~
}
  \item{maxq}{Maximum MA order entertained. Default is 5.
%%     ~~Describe \code{maxq} here~~
}
  \item{h}{The additional past lags used in canonical correlation analysis.
Default is 0.
%%     ~~Describe \code{h} here~~
}
  \item{crit}{Type-I error of the test
%%     ~~Describe \code{crit} here~~
}
  \item{output}{A logical switch to control the output.
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Nmtx: The table of the numbers of zero canonical correlations, and 
DDmtx: 
the table of diagonal differences.
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
\examples{zt=read.table("flourc.txt")
m1=SCMid(zt)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
