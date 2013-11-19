\name{tfm1}
\alias{tfm1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transfer Function Model
%%  ~~function to do ... ~~
}
\description{Estimation of a general transfer function model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tfm1(y, x, orderN, orderX)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Data vector of dependent variable
%%     ~~Describe \code{y} here~~
}
  \item{x}{Data vector of input (or independent) variable
%%     ~~Describe \code{x} here~~
}
  \item{orderN}{Order (p,d,q) of the disturbance component
%%     ~~Describe \code{orderN} here~~
}
  \item{orderX}{Order (r,s,b) of the transfer function model, where 
r and s are the degrees of denominator and numerator polynomials and 
b is the delay
%%     ~~Describe \code{orderX} here~~
}
}
\details{Perform estimation of a general transfer function model
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{estimate }{Coefficient estimates}
\item{sigma2 }{Residual variance sigma-square}
\item{residuals }{Residual series}
\item{varcoef }{Variance of the estimates}
\item{Nt }{The disturbance series}
}
\references{Box and Jenkins (1976, Chapter 12)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{tfm
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