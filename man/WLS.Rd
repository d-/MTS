\name{WLS}
\alias{WLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Least Squares 
%%  ~~function to do ... ~~
}
\description{Perform weighted least squares estimation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WLS(y, x, size, xpxi = NULL, xpy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Data vector of dependent variable
%%     ~~Describe \code{y} here~~
}
  \item{x}{Data matrix of explanatory variables
%%     ~~Describe \code{x} here~~
}
  \item{size}{Window size used to estimate the weights
%%     ~~Describe \code{size} here~~
}
  \item{xpxi}{Inverse of the X'X matrix
%%     ~~Describe \code{xpxi} here~~
}
  \item{xpy}{X'y matrix
%%     ~~Describe \code{xpy} here~~
}
}
\details{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
