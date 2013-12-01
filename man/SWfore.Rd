\name{SWfore}
\alias{SWfore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stock-Watson Diffusion Index FOrecasts
%%  ~~function to do ... ~~
}
\description{Uses the diffusion index approach of Stock and Watson to 
compute out-of-sample forecasts
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SWfore(y, x, orig, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The scalar variable of interest 
%%     ~~Describe \code{y} here~~
}
  \item{x}{The data matrix (T-by-k) of the observed explanatory variables
%%     ~~Describe \code{x} here~~
}
  \item{orig}{Forecast origin
%%     ~~Describe \code{orig} here~~
}
  \item{m}{The number of diffusion index used
%%     ~~Describe \code{m} here~~
}
}
\details{Performs PCA on X at the forecast origin. Then, fit a linear 
regression model to obtain the coefficients of prediction equation. 
Use the prediction equation to produce forecasts and compute forecast 
errors, if any.
%%  ~~ If necessary, more details than the description above ~~
}
\value{coef: beta coefficient. yhat: predictions. MSE: mean squared error if 
available. loadings: loading matrix. DFindex: diffusion indices
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
da=read.table("m-unempstatesAdj.txt",header=T)
drate=diffM(da)
yt=drate[5:415,1]
x=cbind(drate[4:414,],drate[3:413,],drate[2:412,],drate[1:411,])
m1=SWfore(y,x,350,10)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
