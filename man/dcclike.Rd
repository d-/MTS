\name{dcclike}
\alias{dcclike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood of a DCC Model
%%  ~~function to do ... ~~
}
\description{evaluate the log likelihood function of a DCC model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dcclike(par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{parameter vector
%%     ~~Describe \code{par} here~~
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
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{This is an internal function, not designed to be a general command
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
    rt <- DCCrtn
    m <- mused
    cond.dist <- cond.D
    loccor <- Loccor
    type <- DCCtype
    nT = dim(rt)[1]
    k = dim(rt)[2]
    R0 = cor(rt)
    th0 = Vech(R0)
    theta1 = par[1]
    theta2 = par[2]
    theta0 = 1 - theta1 - theta2
    wk0 = theta0 * th0
    loccor = as.matrix(loccor)
    if (type == "TseTsui") {
        ist = m + 1
        one1 = matrix(1, m, 1)
        TH = one1 \%*\% th0
        icnt = 0
        TH1 = NULL
        for (i in 1:(k - 1)) {
            icnt = icnt + 1
            TH1 = cbind(TH1, rep(1, (nT - m)))
            for (j in (i + 1):k) {
                icnt = icnt + 1
                bi = wk0[icnt] + theta2 * loccor[ist:nT, icnt]
                bb = filter(bi, theta1, "r", init = th0[icnt])
                TH1 = cbind(TH1, bb)
            }
        }
        TH1 = cbind(TH1, rep(1, (nT - m)))
        TH = rbind(TH, TH1)
    }
    else {
        ist = 2
        k1 = k * (k + 1)/2
        TH = matrix(th0, 1, k1)
        TH1 = NULL
        for (i in 1:k1) {
            bi = wk0[i] + theta2 * loccor[ist:nT, i]
            bb = filter(bi, theta1, "r", init = th0[i])
            TH1 = cbind(TH1, bb)
        }
        TH = rbind(TH, TH1)
        for (t in ist:nT) {
            Qt = VechM(TH[t, ])
            dd = sqrt(diag(Qt))
            for (i in 1:k) {
                for (j in 1:i) {
                  Qt[i, j] = Qt[i, j]/(dd[i] * dd[j])
                }
            }
            TH[t, ] = Vech(Qt)
        }
    }
    llike = 0
    if (cond.dist == "norm") {
        for (t in ist:nT) {
            sigma1 = VechM(TH[t, ])
            l1 = dmvnorm(rt[t, ], mean = rep(0, k), sigma = sigma1, 
                log = TRUE)
            llike = llike - l1
        }
    }
    else {
        v = par[3]
        k1 = log(gamma((v + k)/2))
        k2 = log(gamma(v/2))
        k3 = k * log((v - 2) * pi)/2
        for (t in ist:nT) {
            Rt = VechM(TH[t, ])
            Rt = as.matrix(Rt)
            yt = matrix(rt[t, ], k, 1)
            b1 = chol(Rt)
            b1 = t(b1)
            b1inv = Lminv(b1)
            d1 = prod(diag(b1))^2
            yt = matrix(rt[t, ], k, 1)
            P1 = b1inv \%*\% yt
            tmp1 = t(P1) \%*\% P1
            llike = llike - k1 + k2 + 0.5 * log(d1) + k3 + (v + 
                k)/2 * log(1 + tmp1/(v - 2))
        }
    }
    llike
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
