\name{refVMA}
\alias{refVMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining VMA Models
%%  ~~function to do ... ~~
}
\description{Refines a fitted VMA model by setting insignificant parameters to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refVMA(model, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object from the command VMA or the refVMA command
%%     ~~Describe \code{model} here~~
}
  \item{thres}{A threshold for individual t-ratio of parameter estimate. 
Default is 2.
%%     ~~Describe \code{thres} here~~
}
}
\details{The program simultaneously sets all estimates with t-ratios 
less than the threshold (in modulus) to zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of VMA estimation
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

\seealso{VMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
theta=matrix(c(0.4,0.5,0.0,0.6),2,2)
sig=matrix(c(3,1,1,1),2,2)
m1=VARMAsim(400,malags=c(1),theta=theta,sigma=sig)
zt=m1$series
m2=VMA(zt,q=1,include.mean=F)
m3=refVMA(m2,thres=1.8)
names(m3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
