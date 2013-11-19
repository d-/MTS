\name{SCMmod}
\alias{SCMmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SCM Model Specification
%%  ~~function to do ... ~~
}
\description{For a given set of SCMs and locator of transformation matrix, 
the program specifies a VARMA model via SCM approach for estimation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMmod(order, Ivor, output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{order}{A k-by-2 matrix of the orders of SCM
%%     ~~Describe \code{order} here~~
}
  \item{Ivor}{A k-dimensioal vector to indicate the location of "1" 
for each component z(i,t) in the transformation matrix.
%%     ~~Describe \code{Ivor} here~~
}
  \item{output}{A logical switch to control output.
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Tmtx: specification of the transformation matrix T. 
ARpar: Specification of the VAR parameters. MApar: specification of the 
VMA parameters.
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
\examples{ord=matrix(c(0,1,1,0,0,1),3,2)
Ivor=c(3,1,2)
m1=SCMmod(ord,Ivor,TRUE)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (order, Ivor, output) 
{
    if (!is.matrix(order)) 
        order = as.matrix(order)
    k = dim(order)[1]
    p = max(order[, 1])
    q = max(order[, 2])
    kp = k * p
    kq = k * q
    if (output) 
        cat("VARMA order (p,q) = (", p, ",", q, ")", "\n")
    phi = NULL
    theta = NULL
    if (p > 0) 
        phi = matrix(0, k, kp)
    if (q > 0) 
        theta = matrix(0, k, kq)
    for (i in 1:k) {
        pi = order[i, 1]
        qi = order[i, 2]
        if (pi > 0) {
            for (j in 1:pi) {
                jst = k * (j - 1)
                phi[i, (jst + 1):(jst + k)] = 2
            }
        }
        if (qi > 0) {
            for (j in 1:qi) {
                jst = k * (j - 1)
                theta[i, (jst + 1):(jst + k)] = 2
            }
        }
    }
    for (i in 1:k) {
        pi = order[i, 1]
        qi = order[i, 2]
        for (j in 1:i) {
            pj = order[j, 1]
            qj = order[j, 2]
            mm = min(pi - pj, qi - qj)
            if (mm > 0) {
                for (jj in 1:mm) {
                  jdx = (jj - 1) * k + j
                  theta[i, jdx] = 0
                }
            }
        }
    }
    if (output) {
        cat("VAR matrices: ", "\n")
        print(phi)
        cat("VMA matrices: ", "\n")
        print(theta)
    }
    Tmtx = matrix(2, k, k)
    for (i in 1:k) {
        k1 = Ivor[i]
        Tmtx[i, k1] = 1
    }
    for (i in 2:k) {
        pi = order[i, 1]
        qi = order[i, 2]
        for (j in 1:(i - 1)) {
            k1 = Ivor[j]
            pj = order[j, 1]
            qj = order[j, 2]
            mink = min(pi - pj, qi - qj)
            if (mink > -1) 
                Tmtx[i, k1] = 0
        }
    }
    if (output) {
        cat("Transformation matrix: ", "\n")
        print(Tmtx)
    }
    SCMmod <- list(Tmtx = Tmtx, ARpar = phi, MApar = theta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
