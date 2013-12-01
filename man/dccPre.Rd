\name{dccPre}
\alias{dccPre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Preliminary Fitting of Dynamic Cross Correlations
%%  ~~function to do ... ~~
}
\description{This program fits marginally GARCH models to each component 
of a vector return so that the standardized return series can be obtained
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dccPre(rtn, include.mean = T, p = 0, cond.dist = "norm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensinal asset returns
%%     ~~Describe \code{rtn} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. 
Deafult is to include the mean.
%%     ~~Describe \code{include.mean} here~~
}
  \item{p}{VAR order for the mean equation
%%     ~~Describe \code{p} here~~
}
  \item{cond.dist}{The conditional distribution of the innovations.
Default is Gaussian. 
%%     ~~Describe \code{cond.dist} here~~
}
}
\details{The program uses fGarch package to estimate univariate GARCH 
model for each return series.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Volatility series of each marginal return series
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
\note{fGarch package is used
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-ibmspko-6111.txt",header=T)
rtn=log(da[,2:4]+1)
m1=dccPre(rtn,include.mean=T,p=0)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
