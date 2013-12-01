\name{SCMid2}
\alias{SCMid2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scalar Component Model Identification
%%  ~~function to do ... ~~
}
\description{Provides detailed analysis of scalar component models 
for a specified VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMid2(zt, maxp = 2, maxq = 2, h = 0, crit = 0.05, sseq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{The T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{maxp}{Maximum AR order specified. Default is 2.
%%     ~~Describe \code{maxp} here~~
}
  \item{maxq}{Maximum MA order specified. Default is 2.
%%     ~~Describe \code{maxq} here~~
}
  \item{h}{The additional past lags used in canonical correlation analysis. 
Default is zero.
%%     ~~Describe \code{h} here~~
}
  \item{crit}{Type-I error used in testing
%%     ~~Describe \code{crit} here~~
}
  \item{sseq}{The search sequence for SCM components
%%     ~~Describe \code{sseq} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Tmatrix: The transformation matrix T. SCMorder: The orders of 
SCM.
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

\seealso{SCMid
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{zt=read.table("flourc.txt")
m2=SCMid2(zt)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
