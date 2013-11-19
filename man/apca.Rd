\name{apca}
\alias{apca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Asymptotic Principal Component Analysis
%%  ~~function to do ... ~~
}
\description{Perform asymptotic PCA for a data set
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
apca(da, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{A T-by-k data set matrix, where T is the sample size and k is the 
dimension
%%     ~~Describe \code{da} here~~
}
  \item{m}{The number of common factors
%%     ~~Describe \code{m} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{square root of eigenvalues, loading and common factors
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
rtn=read.table("m-sp100y2011.txt",header=T)
sp100=apca(rtn,3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (da, m) 
{
    if (!is.matrix(da)) 
        da = as.matrix(da)
    if (m < 1) 
        m = 1
    nT = dim(da)[1]
    k = dim(da)[2]
    if (k <= nT) {
        da = t(da)
        nT = dim(da)[1]
        k = dim(da)[2]
    }
    m1 = princomp(t(da), cor = F, rotation = "none")
    print(summary(m1))
    factors = matrix(m1$loadings, nT, nT)
    factors = factors[, 1:m]
    loadings = m1$scores[, 1:m]
    sdev = m1$sdev
    apca <- list(sdev = sdev, factors = factors, loadings = loadings)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line