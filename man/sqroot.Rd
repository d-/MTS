\name{sqroot}
\alias{sqroot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Square Root of a Matrix
%%  ~~function to do ... ~~
}
\description{Computes the positive-definite square root of a matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sqroot(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A positive-definite matrix
%%     ~~Describe \code{m} here~~
}
}
\details{Uses spectral decomposition to compute the square root of 
a positive-definite matrix
%%  ~~ If necessary, more details than the description above ~~
}
\value{Square root matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
mtx=matrix(c(4,1,1,1),2,2)
m1=sqroot(mts)
m1
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m) 
{
    if (!is.matrix(m)) 
        m = as.matrix(m)
    m = (t(m) + m)/2
    e1 = eigen(m)
    dd = diag(sqrt(e1$values))
    P = e1$vectors
    sqRmtx = P \%*\% dd \%*\% t(P)
    sqroot <- list(sqRt = sqRmtx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
