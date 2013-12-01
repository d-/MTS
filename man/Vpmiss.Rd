\name{Vpmiss}
\alias{Vpmiss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial Missing of a VARMA Model
%%  ~~function to do ... ~~
}
\description{Assuming that the data is only partially missing, 
this program estimates those missing values. The model is assumed to be 
known.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Vpmiss(zt, piwgt, sigma, tmiss, mdx, cnst = NULL, output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
%%     ~~Describe \code{zt} here~~
}
  \item{piwgt}{pi-weights of the model in the form 
piwgt[pi0, pi1, pi2, ....]
%%     ~~Describe \code{piwgt} here~~
}
  \item{sigma}{Residual covariance matrix
%%     ~~Describe \code{sigma} here~~
}
  \item{tmiss}{Time index of the partilly missing data point
%%     ~~Describe \code{tmiss} here~~
}
  \item{mdx}{A k-dimensional indicator with "0" denoting missing 
component and ""1" denoting observed value.
%%     ~~Describe \code{mdx} here~~
}
  \item{cnst}{Constant term of the model
%%     ~~Describe \code{cnst} here~~
}
  \item{output}{values of the partially missing data
%%     ~~Describe \code{output} here~~
}
}
\details{
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

\seealso{Vmiss
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("q-gdp-ukcaus.txt",header=T)
gdp=log(da[,3:5])
m1=VAR(gdp,1)
piwgt=m1$Phi; cnst=m1$Ph0; Sig=m1$Sigma
mdx=c(0,1,1)
m2=Vpmiss(gdp,piwgt,Sig,50,mdx,cnst)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
