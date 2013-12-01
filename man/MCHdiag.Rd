\name{MCHdiag}
\alias{MCHdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Conditional Heteroscadastic Model Cheking
%%  ~~function to do ... ~~
}
\description{Apply four portmanteau test statistics to check the 
validity of a fitted multivariate volatility model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MCHdiag(at, Sigma.t, m = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{at}{A T-by-k matrix of residuals for a k-dimensional asset returns
%%     ~~Describe \code{at} here~~
}
  \item{Sigma.t}{The fitted volatility matrices. The dimension is 
T-by-k^2 matrix
%%     ~~Describe \code{Sigma.t} here~~
}
  \item{m}{The number of lags used in the tests. Default is 10.
%%     ~~Describe \code{m} here~~
}
}
\details{The four test statistics are given in Tsay (2014, Chapter 7)
%%  ~~ If necessary, more details than the description above ~~
}
\value{Test statistics and their p-values
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
da=read.table("m-ibmsp-6111.txt",header=T)
rtn=log(da[,2:3]+1)
m1a=BEKK11(rtn)
names(m1a)
Sigma.t=m1a$Sigma.t
at=cbind(rtn[,1]-0.00776,rtn[,2]-0.00565)
MCHdiag(at,Sigma.t)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
