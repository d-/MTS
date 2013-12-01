\name{REGts}
\alias{REGts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regression Model with Time Series Errors
%%  ~~function to do ... ~~
}
\description{Perform maximum likelihood estimation of a multivariate 
linear regression model with time-series errors
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
REGts(zt, p, xt, include.mean = T, fixed = NULL, par = NULL, se.par = NULL, details = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{p}{The VAR order
%%     ~~Describe \code{p} here~~
}
  \item{xt}{A T-by-v data matrix of independent variables, where 
v denotes the number of independent variables (excluding constant 1).
%%     ~~Describe \code{xt} here~~
}
  \item{include.mean}{A logical switch to include the constant term
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used to set parameters to zero
%%     ~~Describe \code{fixed} here~~
}
  \item{par}{Initial parameter estimates of the beta coefficients
%%     ~~Describe \code{par} here~~
}
  \item{se.par}{Standard errors of the parameters in par.
%%     ~~Describe \code{se.par} here~~
}
  \item{details}{A logical switch to control the output
%%     ~~Describe \code{details} here~~
}
}
\details{Perform maximum likelihood estimation of a multivariate linear 
regression model with time series errors. Use multivariate linear regression 
to obtain initial estimates of regression coefficients if not provided
%%  ~~ If necessary, more details than the description above ~~
}
\value{Estimation results of a regression model with time series errors.
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
da=read.table("m-gasoil.txt",header=T)
zt=da[,3:4]; xt=da[,5:6]
m1=REGts(zt,2,xt)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
