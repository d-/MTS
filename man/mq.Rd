\name{mq}
\alias{mq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Ljung-Box Q Statistics
%%  ~~function to do ... ~~
}
\description{Computes the multivariate Ljung-Box statistics for 
cross-correlation matrices
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mq(x, lag = 24, adj = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data matrix of a vector time series or residual series of a 
fitted multivariate model.
%%     ~~Describe \code{x} here~~
}
  \item{lag}{The number of cross-correlation matrices used. Default is 24.
%%     ~~Describe \code{lag} here~~
}
  \item{adj}{Adjustment for the degrees of freedom for the Ljung-Box 
statistics. This is used for residual series. Default is zero.
%%     ~~Describe \code{adj} here~~
}
}
\details{Computes the multivariate Ljung-Box statistics and their 
p-values. For model checking, the subcommand adj can be used to 
adjust the degrees of freedom of the Chi-square statistics.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2013)
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
x=matrix(rnorm(1500),500,3)
mq(x)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, lag = 24, adj = 0) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nr = nrow(x)
    nc = ncol(x)
    g0 = var(x)
    ginv = solve(g0)
    qm = 0
    QM = NULL
    df = 0
    for (i in 1:lag) {
        x1 = x[(i + 1):nr, ]
        x2 = x[1:(nr - i), ]
        g = cov(x1, x2)
        g = g * (nr - i - 1)/(nr - 1)
        h = t(g) \%*\% ginv \%*\% g \%*\% ginv
        qm = qm + nr * nr * sum(diag(h))/(nr - i)
        df = df + nc * nc
        dff = df - adj
        mindeg = nc^2 - 1
        pv = 1
        if (dff > mindeg) 
            pv = 1 - pchisq(qm, dff)
        QM = rbind(QM, c(i, qm, dff, pv))
    }
    pvs = QM[, 4]
    dimnames(QM) = list(names(pvs), c("  m  ", "    Q(m) ", "   df  ", 
        " p-value"))
    cat("Ljung-Box Statistics: ", "\n")
    printCoefmat(QM, digits = 3)
    par(mfcol = c(1, 1))
    plot(pvs, ylim = c(0, 1), xlab = "m", ylab = "prob", main = "p-values of Ljung-Box statistics")
    abline(h = c(0))
    lines(rep(0.05, lag), lty = 2, col = "blue")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
