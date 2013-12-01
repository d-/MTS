\name{iniEST}
\alias{iniEST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial Estimation
%%  ~~function to do ... ~~
}
\description{Performs least squares estimation to obtain initial parameter 
estimates.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
iniEST(y, x, p, q, include.mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Dependent data matrix
%%     ~~Describe \code{y} here~~
}
  \item{x}{Data matrix of explanatory variables
%%     ~~Describe \code{x} here~~
}
  \item{p}{AR order
%%     ~~Describe \code{p} here~~
}
  \item{q}{MA order
%%     ~~Describe \code{q} here~~
}
  \item{include.mean}{A logical switch to include mean vector
%%     ~~Describe \code{include.mean} here~~
}
}
\details{An internal function, not intended as a general command
%%  ~~ If necessary, more details than the description above ~~
}
\value{Initial estimates
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 3)
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
