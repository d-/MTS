\name{BEKK11vol}
\alias{BEKK11vol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Volatility of BEKK(1,1) Model
%%  ~~function to do ... ~~
}
\description{Computes the volatility matrices of a fitted BEKK(1,1) model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BEKK11vol(rtn, est, include.mean = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns
%%     ~~Describe \code{rtn} here~~
}
  \item{est}{Parameters of a BEKK model, including the mean, if any.
%%     ~~Describe \code{est} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. Default is to include the mean
%%     ~~Describe \code{include.mean} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Volatility matrices
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rtn, est, include.mean = T) 
{
    if (!is.matrix(rtn)) 
        rtn = as.matrix(rtn)
    nT = dim(rtn)[1]
    k = dim(rtn)[2]
    if (k > 3) 
        k = 3
    mu1 = 0
    mu2 = 0
    mu3 = 0
    icnt = 0
    if (k == 2) {
        if (include.mean) {
            mu1 = est[1]
            mu2 = est[2]
            icnt = k
        }
        A0 = matrix(c(est[icnt + 1], est[icnt + 2], 0, est[icnt + 
            3]), 2, 2)
        A0A0t = A0 \%*\% A0
        A1 = matrix(c(est[(icnt + 4):(icnt + 7)]), 2, 2)
        B1 = matrix(c(est[(icnt + 8):(icnt + 11)]), 2, 2)
        resi = cbind(rtn[, 1] - mu1, rtn[, 2] - mu2)
        Cov1 = cov(rtn)
        Sig = matrix(c(Cov1), 1, 4)
        res = resi \%*\% t(A1)
        ArchP = A1at(res)
    }
    if (k == 3) {
        if (include.mean) {
            mu1 = est[1]
            mu2 = est[2]
            mu3 = est[3]
            icnt = k
        }
        A0 = matrix(c(est[icnt + 1], est[icnt + 2], est[icnt + 
            3], 0, est[icnt + 4], est[icnt + 5], 0, 0, est[icnt + 
            6]), k, k)
        A0A0t = A0 \%*\% t(A0)
        icnt = icnt + 6
        A1 = matrix(c(est[(icnt + 1):(icnt + 9)]), k, k)
        icnt = icnt + 9
        B1 = matrix(c(est[(icnt + 1):(icnt + 9)]), k, k)
        resi = cbind(rtn[, 1] - mu1, rtn[, 2] - mu2, rtn[, 3] - 
            mu3)
        Cov1 = cov(rtn)
        Sig = matrix(c(Cov1), 1, k^2)
        res = resi \%*\% t(A1)
        ArchP = A1at(res)
    }
    for (t in 2:nT) {
        Sigt = A0A0t + matrix(ArchP[t - 1, ], k, k) + B1 \%*\% 
            matrix(Sig[t - 1, ], k, k) \%*\% t(B1)
        Sigt = (Sigt + t(Sigt))/2
        Sig = rbind(Sig, c(Sigt))
    }
    BEKK11vol <- list(volmtx = Sig)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
