\name{VARpred}
\alias{VARpred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VAR Prediction
%%  ~~function to do ... ~~
}
\description{Computes the forecasts of a VAR model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARpred(model, h = 1, orig = 0, Out.level = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object of a VAR or refVAR command
%%     ~~Describe \code{model} here~~
}
  \item{h}{Forecast horizon, a positive integer
%%     ~~Describe \code{h} here~~
}
  \item{orig}{Forecast origin. Default is zero meaning the forecast origin 
is the last data point
%%     ~~Describe \code{orig} here~~
}
  \item{Out.level}{A logical switch to control output
%%     ~~Describe \code{Out.level} here~~
}
}
\details{Computes point forecasts and the associated variances of forecast errors
%%  ~~ If necessary, more details than the description above ~~
}
\value{Points forecasts, standard errors of forecasts, and mean squared 
forecast errors if available
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 2)
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
da=read.table("q-gdp-ukcaus.txt",header=T)
gdp=log(da[,3:5])
zt=diffM(gdp)
m1=VAR(zt,p=2)
VARpred(m1,4)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
