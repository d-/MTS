\name{mq}
\alias{mq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Ljung-Box Q Statistics
%%  ~~function to do ... ~~
}
\description{Computes the multivariate Ljung-Box statistics for 
cross-correlation matrices
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mq(x, lag = 24, adj = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data matrix of a vector time series or residual series of a 
fitted multivariate model.
%%     ~~Describe \code{x} here~~
}
  \item{lag}{The number of cross-correlation matrices used. Default is 24.
%%     ~~Describe \code{lag} here~~
}
  \item{adj}{Adjustment for the degrees of freedom for the Ljung-Box 
statistics. This is used for residual series. Default is zero.
%%     ~~Describe \code{adj} here~~
}
}
\details{Computes the multivariate Ljung-Box statistics and their 
p-values. For model checking, the subcommand adj can be used to 
adjust the degrees of freedom of the Chi-square statistics.
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x=matrix(rnorm(1500),500,3)
mq(x)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
