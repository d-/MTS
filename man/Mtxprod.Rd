\name{Mtxprod}
\alias{Mtxprod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix Product
%%  ~~function to do ... ~~
}
\description{Compute the product of two matrices
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Mtxprod(Mtx, sMtx, p, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mtx}{
%%     ~~Describe \code{Mtx} here~~
}
  \item{sMtx}{
%%     ~~Describe \code{sMtx} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{P}{
%%     ~~Describe \code{P} here~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Mtx, sMtx, p, P) 
{
    if (!is.matrix(Mtx)) 
        Mtx = as.matrix(Mtx)
    if (!is.matrix(sMtx)) 
        sMtx = as.matrix(sMtx)
    k = dim(Mtx)[1]
    kp = dim(Mtx)[2]
    kP = dim(sMtx)[2]
    pMtx = Mtx
    for (i in 1:P) {
        ii = (i - 1) * k
        m2 = sMtx[, (ii + 1):(ii + k)]
        pMtx = cbind(pMtx, m2)
        for (j in 1:p) {
            jdx = (j - 1) * k
            m1 = Mtx[, (jdx + 1):(jdx + k)]
            pMtx = cbind(pMtx, -m1 \%*\% m2)
        }
    }
    pMtx
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
