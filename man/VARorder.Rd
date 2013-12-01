\name{VARorder}
\alias{VARorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR Order Specification
%%  ~~function to do ... ~~
}
\description{Computes information criteria and the sequential Chi-square 
statistics for a vector autoregressive process
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARorder(x, maxp = 13, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data matrix of dimension T-by-k with T being the sample size and 
k the number of time series
%%     ~~Describe \code{x} here~~
}
  \item{maxp}{The maximum VAR order entertained. Default is 13.
%%     ~~Describe \code{maxp} here~~
}
  \item{output}{A logical switch to control the output. Default is to provide output
%%     ~~Describe \code{output} here~~
}
}
\details{For a given maxp, the command computes Akaike, Bayesin and 
Hannan-Quinn informtation criteria for various VAR models using the 
data from t=maxp+1 to T. It also computes the Tiao-Box sequential Chi-square 
statistics and their p-values.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2013)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{VARorderI
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{da=read.table("q-gdp-ukcaus.txt",header=T)
zt=diffM(log(da[,3:5]))
VARorder(zt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
