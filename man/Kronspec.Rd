\name{Kronspec}
\alias{Kronspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kronecler Index Specification
%%  ~~function to do ... ~~
}
\description{For a given set of Kronecker indices, the program 
specifies a VARMA model. It gives details of parameter specification.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Kronspec(kdx, output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kdx}{A set of Kronecker indices
%%     ~~Describe \code{kdx} here~~
}
  \item{output}{A logical switch to control output.  Deafult is with 
output.
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{PhiID: specification of the AR matrix polynomial. 
ThetaID: specification of the MA matrix polynomial.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 4)
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
\examples{kdx=c(2,1,1)
m1=Kronspec(kdx)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (kdx, output = TRUE) 
{
    k = length(kdx)
    if (output) {
        cat("Kronecker indices: ", kdx, "\n")
        cat("Dimension: ", k, "\n")
    }
    KK = sort(kdx, index.return = T)
    idx = KK$ix
    p = KK$x[k]
    q = p
    mx = (p + 1) * k
    Theta = matrix(2, k, k * (q + 1))
    for (i in 1:k) {
        Theta[i, i] = 1
        if (kdx[i] < q) {
            jj = (kdx[i] + 1) * k
            Theta[i, (jj + 1):mx] = 0
        }
    }
    if (k > 1) {
        for (i in 1:(k - 1)) {
            Theta[i, (i + 1):k] = 0
        }
    }
    Phi = Theta
    if (k > 1) {
        for (i in 2:k) {
            for (j in 1:(i - 1)) {
                if (kdx[j] <= kdx[i]) 
                  Phi[i, j] = 0
            }
        }
    }
    Theta[1:k, 1:k] = Phi[1:k, 1:k]
    for (i in 1:k) {
        for (j in 1:k) {
            if (kdx[i] > kdx[j]) {
                for (ii in 1:(kdx[i] - kdx[j])) Phi[i, ii * k + 
                  j] = 0
            }
        }
    }
    if (output) {
        cat("Notation: ", "\n")
        cat(" 0: fixed to 0", "\n")
        cat(" 1: fixed to 1", "\n")
        cat(" 2: estimation", "\n")
        cat("AR coefficient matrices: ", "\n")
        print(Phi)
        cat("MA coefficient matrices: ", "\n")
        print(Theta)
    }
    Kronspec <- list(PhiID = Phi, ThetaID = Theta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
