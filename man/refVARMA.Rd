\name{refVARMA}
\alias{refVARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining VARMA Estimation
%%  ~~function to do ... ~~
}
\description{Refines a fitted VARMA model by setting insignificant 
estimates to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refVARMA(model, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object from the command VARMA or the command 
refVARMA
%%     ~~Describe \code{model} here~~
}
  \item{thres}{A threshold value for individual t-ratio of the estimates.
%%     ~~Describe \code{thres} here~~
}
}
\details{The program simultaneously sets estimates with t-ratios less 
than the threshold (in modulus) to zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of VARMA estimates
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

\seealso{VARMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("q-gdp-ukcaus.txt",header=T)
gdp=log(da[,3:5])
zt=diffM(gdp)
m1=VARMA(zt,p=1,q=1)
m2=refVARMA(m1,thres=1.0)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
