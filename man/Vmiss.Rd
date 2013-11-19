\name{Vmiss}
\alias{Vmiss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Missing Values of A Vector Time Series
%%  ~~function to do ... ~~
}
\description{Assuming that the model is known, this program estimates 
the value of a missing data point. The whole data point is missing.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Vmiss(zt, piwgt, sigma, tmiss, cnst = NULL, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{piwgt}{The pi-weights of a VARMA model defined as
piwgt=[pi0, pi1, pi2, ....]
%%     ~~Describe \code{piwgt} here~~
}
  \item{sigma}{Positive definite covariance matrix of the innovations
%%     ~~Describe \code{sigma} here~~
}
  \item{tmiss}{Time index of the missing data point
%%     ~~Describe \code{tmiss} here~~
}
  \item{cnst}{Constant term of the model
%%     ~~Describe \code{cnst} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{Use least squares method to estimate a missing data point. 
The missing is random.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Estimates of the missing values
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

\seealso{Vpmiss
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	da=read.table("q-gdp-ukcaus.txt",header=T)
	gdp=log(da[,3:5])
	m1=VAR(gdp,3)
	piwgt=m1$Phi; Sig=m1$Sigma; cnst=m1$Ph0
	m2=Vmiss(gdp,piwgt,Sig,50,cnst)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line