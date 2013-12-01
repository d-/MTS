\name{diffM}
\alias{diffM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Difference of multivariate time series
%%  ~~function to do ... ~~
}
\description{Performs the difference operation of a vector time series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diffM(zt, d = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A vector time series (T by k, with sample size T and dimension k)
%%     ~~Describe \code{zt} here~~
}
  \item{d}{Order of differencing. Deafult is d=1. 
%%     ~~Describe \code{d} here~~
}
}
\details{When d = 1, the command is equivalent to 
apply(xt,2,diff)
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
\author{Ruey S Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{da=read.table("q-gdp-ukcaus.txt",header=T)
zt=log(da[,3:5])
xt=diffM(zt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
