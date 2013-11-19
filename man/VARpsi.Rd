\name{VARpsi}
\alias{VARpsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Psi-Weights of a VAR Model
%%  ~~function to do ... ~~
}
\description{Computes the psi-weight matrices of a VAR model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARpsi(Phi, lag = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficients in the form 
Phi=[Phi1, Phi2, ..., Phip]
%%     ~~Describe \code{Phi} here~~
}
  \item{lag}{Number of psi-weight lags 
%%     ~~Describe \code{lag} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Psi-weights of a VAR model
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 2)
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
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
m1=VARpsi(p1,4)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Phi, lag = 5) 
{
    k = nrow(Phi)
    m = ncol(Phi)
    p = floor(m/k)
    Si = diag(rep(1, k))
    if (p < 1) 
        p = 1
    if (lag < 1) 
        lag = 1
    for (i in 1:lag) {
        if (i <= p) {
            idx = (i - 1) * k
            tmp = Phi[, (idx + 1):(idx + k)]
        }
        else {
            tmp = matrix(0, k, k)
        }
        jj = i - 1
        jp = min(jj, p)
        if (jp > 0) {
            for (j in 1:jp) {
                jdx = (j - 1) * k
                idx = (i - j) * k
                w1 = Phi[, (jdx + 1):(jdx + k)]
                w2 = Si[, (idx + 1):(idx + k)]
                tmp = tmp + w1 \%*\% w2
            }
        }
        Si = cbind(Si, tmp)
    }
    VARpsi <- list(psi = Si)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
