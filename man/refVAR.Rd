\name{refVAR}
\alias{refVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining a VAR Model
%%  ~~function to do ... ~~
}
\description{Refine a fitted VAR model by removing simultaneously insignificant parameters
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refVAR(model, fixed = NULL, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object of the command VAR or the refVAR command
%%     ~~Describe \code{model} here~~
}
  \item{fixed}{A logical matrix for VAR polynomial structure
%%     ~~Describe \code{fixed} here~~
}
  \item{thres}{Threshold used to set parameter to zero. Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{Refine a VAR fitting by setting all estimates with t-ratio 
less than the threshold (in modulus) to zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of VAR estimation
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

\seealso{VAR
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("q-gdp-ukcaus.txt",header=T)
gdp=log(da[,3:5])
zt=diffM(gdp)
m1=VAR(zt,3)
m2=refVAR(m1,thres=1.0)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
