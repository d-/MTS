\name{VARXorder}
\alias{VARXorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VARX Order Specification
%%  ~~function to do ... ~~
}
\description{Specifies the orders of a VARX model, including AR order 
and the number of lags of exogenous variables
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARXorder(x, exog, maxp = 13, maxm = 3, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{x} here~~
}
  \item{exog}{A T-by-m data matrix of exogenous variables
%%     ~~Describe \code{exog} here~~
}
  \item{maxp}{The maximum VAR order detertained
%%     ~~Describe \code{maxp} here~~
}
  \item{maxm}{The maximum lags of exogenous variables entertained
%%     ~~Describe \code{maxm} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{Computes the information criteria of a VARX process
%%  ~~ If necessary, more details than the description above ~~
}
\value{Information criteria
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
m1=VARXorder(zt,xt,maxp=11,maxm=2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
