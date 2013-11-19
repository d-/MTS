\name{dccRho}
\alias{dccRho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Correlation Matrices of DCC Model
%%  ~~function to do ... ~~
}
\description{Given the parameter values, the program compute the conditional correlation matrices 
of a DCC model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dccRho(par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{Parameter vector consisting of theta_1 and theta_2
%%     ~~Describe \code{par} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Conditional cross-correlation matrices
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
\note{An internal program, not intended as a general R command
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
        for (t in ist:nT) {
            v1 = wk0 + theta1 * TH[t - 1, ] + theta2 * loccor[t, 
                ]
            TH = rbind(TH, v1)
        }
    }
    else {
        ist = 2
        TH = matrix(th0, 1, k * (k + 1)/2)
        for (t in ist:nT) {
            v1 = wk0 + theta1 * TH[t - 1, ] + theta2 * loccor[t, 
                ]
            TH = rbind(TH, v1)
        }
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
    rho.t = NULL
    for (t in 1:nT) {
        Rt = VechM(TH[t, ])
        rho.t = rbind(rho.t, c(Rt))
    }
    rho.t
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
