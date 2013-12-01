\name{refREGts}
\alias{refREGts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining a Regression Model with Time Series Errors
%%  ~~function to do ... ~~
}
\description{Refines a fitted REGts by setting simultaneously parameters 
with t-ratios less than the threshold (in modulus) to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refREGts(m1, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An output object from the REGts command or refREGts command
%%     ~~Describe \code{m1} here~~
}
  \item{thres}{Threshold value for individual t-ratio. Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a REGts model
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

\seealso{refVAR, refVARMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-gasoil.txt",header=T)
zt=da[,3:4]; xt=da[,5:6]
m1=REGts(zt,2,xt)
m2=refREGts(1m,thres=1.2)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
