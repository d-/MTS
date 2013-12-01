\name{VMAs}
\alias{VMAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VMA Model with Selected Lags
%%  ~~function to do ... ~~
}
\description{Performs conditional maximum likelihood estimation 
of a VMA model with selected lags in the model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMAs(da, malags, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k matrix of a k-dimensional time series with T observations
%%     ~~Describe \code{da} here~~
}
  \item{malags}{A vector consisting of non-zero MA lags
%%     ~~Describe \code{malags} here~~
}
  \item{include.mean}{A logical switch to include the mean vector
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to fix coefficients to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch concerning initial estimation
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control output level
%%     ~~Describe \code{details} here~~
}
  \item{thres}{A threshold value for setting coefficient estimates to zero
%%     ~~Describe \code{thres} here~~
}
}
\details{A modified version of VMA model by allowing the users to select 
non-zero MA lags 
%%  ~~ If necessary, more details than the description above ~~
}
\value{VMA estimation
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

\seealso{VMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-ibmsp-6111.txt",header=T)
rtn=log(da[,2:3]+1)
m1=VMAs(rtn,malags=c(1,3))
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
