\name{BEKK11}
\alias{BEKK11}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{BEKK Model 
%%  ~~function to do ... ~~
}
\description{Estimation of BEKK(1,1) Model. Only k = 2 or 3 is 
available
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BEKK11(rt, include.mean = T, cond.dist = "normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{A T-by-k data matrix of k-dimensional asset returns
%%     ~~Describe \code{rt} here~~
}
  \item{include.mean}{A logical switch to include a constant vector in the 
mean equation
%%     ~~Describe \code{include.mean} here~~
}
  \item{cond.dist}{Conditional innovation distribution. 
Deafult is Gaussian.
%%     ~~Describe \code{cond.dist} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 7)
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
> da=read.table("m-ibmsp-6111.txt",header=T)
> rtn=log(da[,2:3]+1)
> m1=BEKK11(rtn)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
