\name{VARMApred}
\alias{VARMApred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VARMA Prediction
%%  ~~function to do ... ~~
}
\description{Compute forecasts and their associate forecast error covariances 
of a VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMApred(model, h = 1, orig = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A fitted VARMA model
%%     ~~Describe \code{model} here~~
}
  \item{h}{Number of steps of frecasts, i.e., forecast horizon.
%%     ~~Describe \code{h} here~~
}
  \item{orig}{Forecast origin. Default is end of the sample.
%%     ~~Describe \code{orig} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Forecasts and their standard errors
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m2=VARMA(zt,p=3,q=1) # where zt is the growth rates of PCE and DSPI.
m2p=VARMApred(m2,h=6)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
