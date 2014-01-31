\name{archTest}
\alias{archTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ARCH test
%%  ~~function to do ... ~~
}
\description{
Perform tests to check the conditional heteroscedasticity in a time series. The Ljung-Box statistics of squared series and a rank-based Ljung-Box 
test are used.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
archTest(rt, lag = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{A scalar time series. If rt is a matrix, only the first column is used.
%%     ~~Describe \code{rt} here~~
}
  \item{lag}{
The number of lags of ACF used in the Ljung-Box statistics. The default is 10.
%%     ~~Describe \code{lag} here~~
}
}
\details{
The Ljung-Box statistics based on the squared series are computed first. 
The rank series of the squared time series is than used to test the 
conditional heteroscedasticity.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Test results and their p-values
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ Ruey Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{MarchTest
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
rt=rnorm(200)
archTest(rt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rt, lag = 10) 
{
    at = rt
    if (is.matrix(at)) 
        at = at[, 1]
    m1 = acf(at^2, lag.max = lag, plot = F)
    acf = m1$acf[2:(lag + 1)]
    nT = length(at)
    c1 = c(1:lag)
    deno = rep(nT, lag) - c1
    Q = sum(acf^2/deno) * nT * (nT + 2)
    pv1 = 1 - pchisq(Q, lag)
    cat("Q(m) of squared series(LM test): ", "\n")
    cat("Test statistic: ", Q, " p-value: ", pv1, "\n")
    rk = rank(at^2)
    m2 = acf(rk, lag.max = lag, plot = F)
    acf = m2$acf[2:(lag + 1)]
    mu = -(rep(nT, lag) - c(1:lag))/(nT * (nT - 1))
    v1 = rep(5 * nT^4, lag) - (5 * c(1:lag) + 9) * nT^3 + 9 * 
        (c(1:lag) - 2) * nT^2 + 2 * c(1:lag) * (5 * c(1:lag) + 
        8) * nT + 16 * c(1:lag)^2
    v1 = v1/(5 * (nT - 1)^2 * nT^2 * (nT + 1))
    QR = sum((acf - mu)^2/v1)
    pv2 = 1 - pchisq(QR, lag)
    cat("Rank-based Test: ", "\n")
    cat("Test statistic: ", QR, " p-value: ", pv2, "\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
