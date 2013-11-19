\name{VMAorder}
\alias{VMAorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Order specification of a VMA process
%%  ~~function to do ... ~~
}
\description{Performs multivariate Ljung-Box tests to specify the 
order of a VMA process
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VMAorder(x, lag = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data matrix of the observed k-dimensional time series. 
Each column represents a time series.
%%     ~~Describe \code{x} here~~
}
  \item{lag}{The maximum VMA order entertained. Default is 20.
%%     ~~Describe \code{lag} here~~
}
}
\details{For a given lag, the command computes the Ljung-Box 
statistic for testing rho_j = ... = rho_lag = 0, where 
j = 1, 2, ..., lag. For a VMA(q) process, the 
Ljung-Box statistics should be significant for the first q lags, 
and insignificant thereafter.
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
\examples{zt=matrix(rnorm(600),200,3)
VMAorder(zt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, lag = 20) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nr = dim(x)[1]
    nc = dim(x)[2]
    g0 = var(x)
    ginv = solve(g0)
    qm = NULL
    for (i in 1:lag) {
        x1 = x[(i + 1):nr, ]
        x2 = x[1:(nr - i), ]
        g = cov(x1, x2)
        g = g * (nr - i - 1)/(nr - 1)
        h = t(g) \%*\% ginv \%*\% g \%*\% ginv
        qmi = nr * nr * sum(diag(h))/(nr - i)
        qm = c(qmi, qm)
    }
    tst = rev(cumsum(qm))
    ksq = nc * nc
    df = ksq * lag
    QM = NULL
    for (i in 1:lag) {
        pv = 1 - pchisq(tst[i], df)
        QM = rbind(QM, c(i, tst[i], pv))
        df = df - ksq
    }
    pvs = QM[, 3]
    dimnames(QM) = list(names(pvs), c("  j  ", "  Q(j,m) ", " p-value"))
    cat("Q(j,m) Statistics: ", "\n")
    printCoefmat(QM, digits = 3)
    par(mfcol = c(1, 1))
    plot(pvs, ylim = c(0, 1), xlab = "j", ylab = "prob", main = "p-values: Q(j,m) Statistics")
    abline(h = c(0))
    lines(rep(0.05, lag), lty = 2, col = "blue")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
