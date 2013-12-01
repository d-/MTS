\name{VARX}
\alias{VARX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR Model with Exogenous Variables
%%  ~~function to do ... ~~
}
\description{Estimation of a VARX model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARX(zt, p, xt = NULL, m = 0, include.mean = T, fixed = NULL, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{p}{The VAR order
%%     ~~Describe \code{p} here~~
}
  \item{xt}{A T-by-m data matrix of m exogenous variables
%%     ~~Describe \code{xt} here~~
}
  \item{m}{The number of exogenous variables
%%     ~~Describe \code{m} here~~
}
  \item{include.mean}{A logical switch to include the constant vector. 
Default is to include the constant.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix for setting parameters to zero.
%%     ~~Describe \code{fixed} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{Performs least squares estimation of a VARX(p,s) model
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a VARX model, including model and data
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
m1=VARX(zt,2,xt,1)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
