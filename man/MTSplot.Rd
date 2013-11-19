\name{MTSplot}
\alias{MTSplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Time Series Plot
%%  ~~function to do ... ~~
}
\description{Provides time plots of a vector time series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MTSplot(data, caltime = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Vector time series data matrix
%%     ~~Describe \code{data} here~~
}
  \item{caltime}{Calendar time. Default is NULL, that is, using 
time index
%%     ~~Describe \code{caltime} here~~
}
}
\details{Provides time plots of a vector time series. The output frame 
depends on the dimension of the time series
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
xt=matrix(rnorm(1500),500,3)
MTSplot(xt)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, caltime = NULL) 
{
    if (!is.matrix(data)) 
        data = as.matrix(data)
    if (is.ts(data)) {
        plot(data)
    }
    else {
        nT = dim(data)[1]
        tdx = c(1:nT)
        if (length(caltime) > 1) 
            tdx = caltime
        k = dim(data)[2]
        par(mfcol = c(k, 1))
        k1 = k
        if (k == 4) {
            par(mfcol = c(2, 2))
            k1 = 4
        }
        if (k > 4) {
            par(mfcol = c(3, 2))
            k1 = 6
        }
        jcnt = 0
        for (j in 1:k) {
            plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
                type = "l")
            jcnt = jcnt + 1
            if ((jcnt == k1) && (k > 3)) {
                jcnt = 0
                cat("Hit return for more plots: ", "\n")
                readline()
            }
        }
    }
    par(mfcol = c(1, 1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
