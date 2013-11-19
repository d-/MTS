\name{hfactor}
\alias{hfactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constrained Factor Model
%%  ~~function to do ... ~~
}
\description{Performs constraint factor model analysis with a given 
constrained matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hfactor(X, H, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A T-by-k data matrix of the  observed k-dimensional time series
%%     ~~Describe \code{X} here~~
}
  \item{H}{The constrained matrix with each column representing a constraint
%%     ~~Describe \code{H} here~~
}
  \item{r}{The number of common factor
%%     ~~Describe \code{r} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of the traditional PCA and constrained factor models are given
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 6). Tsai and Tsay (2010, JASA)
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
da=read.table("m-tenstocks.txt",header=T)
rtn=log(da[,2:11]+1) # compute log returns
h1=c(1,1,1,1,rep(0,6)) # specify the constraints
h2=c(0,0,0,0,1,1,1,0,0,0)
h3=c(rep(0,7),1,1,1)
H=cbind(h1,h2,h3)
m1=hfactor(rtn,H,3)
[1] "Data are individually"
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
