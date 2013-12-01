\name{dccFit}
\alias{dccFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dynamic Cross-Correlation Fitting
%%  ~~function to do ... ~~
}
\description{Fits a DCC model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dccFit(rt, type = "TseTsui", theta = c(0.9, 0.02), ub = c(0.92, 0.079999), lb = c(0.4, 1e-04), cond.dist = "std", df = 7, m = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{The T-by-k data matrix of k-dimensional asset returns
%%     ~~Describe \code{rt} here~~
}
  \item{type}{A logical switch to specify the type of DCC model. 
Type="TseTsui" 
for Tse and Tsui's DCC model.
%%     ~~Describe \code{type} here~~
}
  \item{theta}{The initial parameter values for theta1 and theta2
%%     ~~Describe \code{theta} here~~
}
  \item{ub}{Upper bound of parameters
%%     ~~Describe \code{ub} here~~
}
  \item{lb}{Lower bound of parameters
%%     ~~Describe \code{lb} here~~
}
  \item{cond.dist}{Conditional innovation distribution
%%     ~~Describe \code{cond.dist} here~~
}
  \item{df}{degrees of freedom
%%     ~~Describe \code{df} here~~
}
  \item{m}{For Tse and Tsui method, m denotes the number of returns 
used in local correlation matric estimation
%%     ~~Describe \code{m} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Time varying cross-correlation matrices
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 7)
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
# continued from dccPre command demonstration
rtn1=m1$sresi
Vol=m1$marVol
m2=dccFit(rtn1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
