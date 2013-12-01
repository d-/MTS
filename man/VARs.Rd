\name{VARs}
\alias{VARs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR Model of Selected Lags
%%  ~~function to do ... ~~
}
\description{This is a modified version of VAR command by allowing the 
users to specify which AR lags to be included in the model. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARs(x, lags, include.mean = T, output = T, fixed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of k-dimensional time series with T observations
%%     ~~Describe \code{x} here~~
}
  \item{lags}{A vector of non-zero AR lags. For instance, lags=c(1,3) 
denotes a VAR(3) model with Phi2 = 0.
%%     ~~Describe \code{lags} here~~
}
  \item{include.mean}{A logical switch to include the mean vector
%%     ~~Describe \code{include.mean} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
  \item{fixed}{A logical matrix to fix parameters to zero.
%%     ~~Describe \code{fixed} here~~
}
}
\details{Performs VAR estimation by allowing certain lag 
coefficient matrices being zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{VAR estimation
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ruey S. Tsay (Chapter 2)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{VAR command
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("q-gdp-ukcaus.txt",header=T)
zt=log(da[,3:5])
m1=VARs(zt,lags=c(1,2,4))
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
