\name{msqrt}
\alias{msqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix Square Root
%%  ~~function to do ... ~~
}
\description{Compute the symmetric square root of a positive definite 
matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
msqrt(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A positive definite matrix
%%     ~~Describe \code{M} here~~
}
}
\details{Use spectral decomposition to compute the square root of 
a positive definite matrix
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
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{This command is used in some of the MTS functions.
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
function (M) 
{
    if (!is.matrix(M)) 
        M = as.matrix(M)
    n1 = nrow(M)
    if (n1 == 1) {
        Mh = sqrt(M)
        Mhinv = 1/Mh
    }
    if (n1 > 1) {
        M = (M + t(M))/2
        m1 = eigen(M)
        V = m1$vectors
        eiv = sqrt(m1$values)
        L = diag(eiv)
        Linv = diag(1/eiv)
        Mh = V \%*\% L \%*\% t(V)
        Mhinv = V \%*\% Linv \%*\% t(V)
    }
    list(mtxsqrt = Mh, invsqrt = Mhinv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
