\name{VARMA}
\alias{VARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vector Autoregressive Moving-Average Models
%%  ~~function to do ... ~~
}
\description{Performs conditional maximum likelihood estimation of 
a VARMA model. Multivariate Gaussian likelihood function is used.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMA(da, p = 0, q = 0, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data matrix (T-by-k) of a k-dimensional time series with sample 
size T.
%%     ~~Describe \code{da} here~~
}
  \item{p}{AR order
%%     ~~Describe \code{p} here~~
}
  \item{q}{MA order
%%     ~~Describe \code{q} here~~
}
  \item{include.mean}{A logical switch to control estimation of the mean vector.
Dafault is to include the mean in estimation.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix to control zero coefficients in estimation. 
It is mainly used by the command refVARMA.
%%     ~~Describe \code{fixed} here~~
}
  \item{prelim}{A logical switch to control preliminary estimation. 
Deafult is none.
%%     ~~Describe \code{prelim} here~~
}
  \item{details}{A logical switch to control the amount of output.
%%     ~~Describe \code{details} here~~
}
  \item{thres}{A threshold used to set zero parameter constraints based on 
individual t-ratio.
Default is 2.
%%     ~~Describe \code{thres} here~~
}
}
\details{Coefficient estimates, model specification, and residuals are returned.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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

\seealso{refVARMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
d1=read.table("m-pce.txt",header=T); d2=read.table("m-dspi.txt",header=T)
z1=diff(log(d1$pce)); z2=diff(log(d2$dspi))
zt=cbind(z1,z2)
colnames(zt) <- c("pce","dspi")
m2=VARMA(zt,p=3,q=1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
