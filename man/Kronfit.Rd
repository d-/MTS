\name{Kronfit}
\alias{Kronfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of VARMA Model via Kronecker Index
%%  ~~function to do ... ~~
}
\description{Perform estimation of a VARMA model specified by Kronecker 
indices
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Kronfit(da, kidx, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix (T-by-k) of a k-dimensional time series
%%     ~~Describe \code{da} here~~
}
  \item{kidx}{The set of Kronecker indices
%%     ~~Describe \code{kidx} here~~
}
  \item{include.mean}{A logical switch for including the mean vector in 
estimation. Default is to include the mean vector.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used to set zero parameter constraints. 
This is used mainly in the command refKronfit.
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A lofical switch for a preliminary estimation.
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control output.
%%     ~~Describe \code{details} here~~
}
  \item{thres}{A threshold for t-ratio in setting parameter to zero. 
Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of estimation, including details of model specification and 
residuals
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

\seealso{refKronfit
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{zt=read.table("flourc.txt")
kdx=c(1,1,1)
m1=Kronfit(da,kdx)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
