\name{MBcluster}
\alias{MBcluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model-Based Clustering Analysis
%%  ~~function to do ... ~~
}
\description{Perform clustering analysis via the MCMC methods
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MBcluster(data, p, K, mcmc, prior = NULL, differ = FALSE, start.values = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A T-by-k data matrix ofk-dimensional  time series
%%     ~~Describe \code{data} here~~
}
  \item{p}{Univariate AR order
%%     ~~Describe \code{p} here~~
}
  \item{K}{The number of clusters
%%     ~~Describe \code{K} here~~
}
  \item{mcmc}{A vector consisting of the number of burn-in and the total 
number of itertations of MCMC methods.
%%     ~~Describe \code{mcmc} here~~
}
  \item{prior}{prior distribution 
%%     ~~Describe \code{prior} here~~
}
  \item{differ}{Differencing the data
%%     ~~Describe \code{differ} here~~
}
  \item{start.values}{Initial values
%%     ~~Describe \code{start.values} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Parameter estimates for each clusters and members of each cluster
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 6)
%% ~put references to the literature/web site here ~
}
\author{Yongning Wang
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
da=read.table("m-unempstatesAdj.txt",header=T)
zt=diffM(da)
iter=c(2000,5000)
m1=MBcluster(zt,4,4,mcmc=iter)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
