\name{VMA}
\alias{VMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vector Moving Averge Model
%%  ~~function to do ... ~~
}
\description{Performs VMA estimation using the conditional multivariate 
Gaussian likelihood function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMA(da, q = 1, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix of a k-dimensional VMA process with each column 
containing one time series
%%     ~~Describe \code{da} here~~
}
  \item{q}{The order of VMA model
%%     ~~Describe \code{q} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. 
The default is to include the mean vector in estimation.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used to fix parameter to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch to select parameters to be included in estimation
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control the amount of output
%%     ~~Describe \code{details} here~~
}
  \item{thres}{Threshold for t-ratio used to fix parameeter to zero. Deault is 2.
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
\references{Tsay (2014, Chapter 3).
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
\examples{da=read.table("m-dec15678-6111.txt",header=T)
rtn=cbind(da$dec5,da$dec8)
zt=log(rtn+1)*100
m1=VMA(zt,q=1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
