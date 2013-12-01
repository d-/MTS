\name{RLS}
\alias{RLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursive Least squares
%%  ~~function to do ... ~~
}
\description{Compute recursive least squares estimation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RLS(y, x, ist = 30, xpxi = NULL, xpy0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data of dependent variable
%%     ~~Describe \code{y} here~~
}
  \item{x}{data matrix of regressors
%%     ~~Describe \code{x} here~~
}
  \item{ist}{initial number of data points used to start the estimation
%%     ~~Describe \code{ist} here~~
}
  \item{xpxi}{Inverse of the X'X matrix
%%     ~~Describe \code{xpxi} here~~
}
  \item{xpy0}{Initial value of X'y.
%%     ~~Describe \code{xpy0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Recursive least squares estimates
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
\note{This function is used internally, but can be used as a command too.
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
