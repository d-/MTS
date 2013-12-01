\name{sVARMA}
\alias{sVARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Seasonal VARMA Model
%%  ~~function to do ... ~~
}
\description{Performs conditional maximum likelihood estimation of 
a seasonal VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sVARMA(da, order, sorder, s, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 0.8, switch = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k data matrix of a k-dimensional seasonal time series
%%     ~~Describe \code{da} here~~
}
  \item{order}{Regular order (p,d,q) of the model
%%     ~~Describe \code{order} here~~
}
  \item{sorder}{Seasonal order (P,D,Q) of the model
%%     ~~Describe \code{sorder} here~~
}
  \item{s}{Seasonality. s=4 for quarterly data and s=12 for monthly series
%%     ~~Describe \code{s} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. Deafult is 
to include the mean
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to set zero parameter constraints
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A lofical switch for preliminary estimation. 
Default is none.
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch for output
%%     ~~Describe \code{details} here~~
}
  \item{thres}{Threshold value for individual t-ratio used to set zero parameter constraints
%%     ~~Describe \code{thres} here~~
}
  \item{switch}{A logical switch to exchange the ordering of the regular 
and seasonal VMA factors. Default is theta(B)*Theta(B).
%%     ~~Describe \code{switch} here~~
}
}
\details{Estimation of a sesonal VARMA model
%%  ~~ If necessary, more details than the description above ~~
}
\value{Estimation results of a seasonal VARMA model, including model specification
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 6)
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
da=read.table("m-hsoldhst6312.txt",header=T)
zt=log(da[,3:4])
m3=sVARMA(zt,order=c(0,1,2),sorder=c(0,1,1),s=12)
names(m3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
