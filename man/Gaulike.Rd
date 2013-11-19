\name{Gaulike}
\alias{Gaulike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gaussian Likelihood Function
%%  ~~function to do ... ~~
}
\description{Compute the log-likelihood function of a Gaussian model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Gaulike(par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{parameter vector
%%     ~~Describe \code{par} here~~
}
}
\details{An internal function, not to be used as a general command.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014)
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (par) 
{
    p = OrderN[1]
    q = OrderN[3]
    r = OrderX[1]
    s = OrderX[2]
    b = OrderX[3]
    c0 = par[1]
    ome = par[2:(2 + s)]
    if (r > 0) 
        del = par[(2 + s + 1):(2 + s + r)]
    if (p > 0) 
        phi = par[(3 + r + s):(r + s + 2 + p)]
    if (q > 0) 
        theta = par[(3 + r + s + p):(2 + r + s + p + q)]
    N = length(dY)
    ist = r + 1
    N1t = dY - c0
    Nt = N1t[ist:N]
    if (r > 0) {
        for (j in 1:r) {
            Nt = Nt - del[j] * N1t[(ist - j):(N - j)]
        }
    }
    ist = b + s + 1
    N = length(Nt)
    N1t = Nt[ist:N] - ome[1] * dX[(ist - b + r):(N - b + r)]
    if (s > 0) {
        for (j in 1:s) {
            N1t = N1t - ome[j + 1] * dX[(ist - j - b + r):(N - 
                j - b + r)]
        }
    }
    N1 = length(N1t)
    resi = N1t[(p + 1):N1]
    if (p > 0) {
        for (j in 1:p) {
            resi = resi - phi[j] * N1t[(p + 1 - j):(N1 - j)]
        }
    }
    if (q > 0) 
        resi = filter(resi, theta, method = "r", init = rep(0, 
            q))
    Gaulike = resi
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
