\name{SCMfit}
\alias{SCMfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting a VARMA Model via SCM Approch
%%  ~~function to do ... ~~
}
\description{Perform estimation of a VARMA model specified via 
the SCM approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMfit(da, scms, Tdx, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 1, ref = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{The T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{scms}{A k-by-2 matrix of the orders of SCMs
%%     ~~Describe \code{scms} here~~
}
  \item{Tdx}{A k-dimensional vector for locating "1" of each row in the 
transformation matrix.
%%     ~~Describe \code{Tdx} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. 
Default is to include mean vector.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to set parameters to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch for preliminary estimation. 
Default is false. 
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control details of output
%%     ~~Describe \code{details} here~~
}
  \item{thres}{Threshold for individual t-ratio when setting parameters 
to zero. Default is 1.
%%     ~~Describe \code{thres} here~~
}
  \item{ref}{A switch to use SCMmod in model specification.
%%     ~~Describe \code{ref} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of estimation, including model specification and residuals
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 4)
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
\examples{zt=read.table("flourc.txt")
ord=matrix(c(1,1,2,0,0,0),3,2)
Tdx=c(2,3,1)
m1=SCMfit(zt,ord,Tdx)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
