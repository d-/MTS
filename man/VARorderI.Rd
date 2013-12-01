\name{VARorderI}
\alias{VARorderI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR order specification I
%%  ~~function to do ... ~~
}
\description{This program is similar to VARorder, but it uses 
observations from t=p+1 to T to compute the information criteria for 
a given VAR(p) model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARorderI(x, maxp = 13, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of vector time series
%%     ~~Describe \code{x} here~~
}
  \item{maxp}{The maximum VAR order entertained
%%     ~~Describe \code{maxp} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{For a given VAR(p) model, the program uses observations 
from t=p+1 to T to compute the information criteria.
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
\author{Ruey S Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{VARorder
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
