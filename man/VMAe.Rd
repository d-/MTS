\name{VMAe}
\alias{VMAe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of VMA Model with Exact likelihood
%%  ~~function to do ... ~~
}
\description{Estimation of a VMA(q) model using the exact likelihood method. 
Multivariate Gaussian likelihood function is used.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMAe(da, q = 1, include.mean = T, coef0 = NULL, secoef0 = NULL, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix (T-by-k) for a k-dimensional VMA process
%%     ~~Describe \code{da} here~~
}
  \item{q}{The order of a VMA model
%%     ~~Describe \code{q} here~~
}
  \item{include.mean}{A logicak switch to include the mean vector in estimation. Default is to include the mean vector.
%%     ~~Describe \code{include.mean} here~~
}
  \item{coef0}{
%%     ~~Describe \code{coef0} here~~
}
  \item{secoef0}{
%%     ~~Describe \code{secoef0} here~~
}
  \item{fixed}{
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{
%%     ~~Describe \code{details} here~~
}
  \item{thres}{
%%     ~~Describe \code{thres} here~~
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
\author{
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
