\name{VechM}
\alias{VechM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix constructed from output of the Vech Command
%%  ~~function to do ... ~~
}
\description{Restores the symmetric matrix from the Vech command
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VechM(vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{A vector
%%     ~~Describe \code{vec} here~~
}
}
\details{This command re-construct a symmetric matrix from output 
of the Vech command
%%  ~~ If necessary, more details than the description above ~~
}
\value{A symmetric matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Appendix A)
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
v1=c(2,1,3)
m1=VechM(v1)
m1
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vec) 
{
    m = length(vec)
    k = (-1 + sqrt(1 + 8 * m))/2
    mtx = diag(rep(1, k))
    ist = 0
    for (j in 1:k) {
        mtx[j:k, j] = vec[(ist + 1):(ist + k - j + 1)]
        ist = ist + (k - j + 1)
    }
    for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
            mtx[i, j] = mtx[j, i]
        }
    }
    mtx
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
