\name{VARXpred}
\alias{VARXpred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VARX Model Prediction
%%  ~~function to do ... ~~
}
\description{Computes point forecasts of a VARX model. The values of 
exogenous variables must be given.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARXpred(m1, newxt = NULL, hstep = 1, orig = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An output object of VARX or refVARX command
%%     ~~Describe \code{m1} here~~
}
  \item{newxt}{The data matrix of exogenous variables needed in forecasts.
%%     ~~Describe \code{newxt} here~~
}
  \item{hstep}{Forecast horizon
%%     ~~Describe \code{hstep} here~~
}
  \item{orig}{Forecast origin. Default is 0, meaning the last data point.
%%     ~~Describe \code{orig} here~~
}
}
\details{Uses the provided exogenous variables and the model to 
compute forecasts
%%  ~~ If necessary, more details than the description above ~~
}
\value{Point forecasts
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
T=dim(da)[1]
ori=T-5
zt1=zt[1:ori,]; xt1=xt[1:ori,]
m1=VARX(zt1,2,xt1,1)
nxt=xt[(ori+1):T,]
m2=VARXpred(m1,newxt=nxt,hstep=5)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
