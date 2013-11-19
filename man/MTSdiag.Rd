\name{MTSdiag}
\alias{MTSdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Time Series Diagnostic Check
%%  ~~function to do ... ~~
}
\description{Performs model checking for a fitted multivariate time series model, including residual cross-correlation matrices, 
multivariate Ljung-Box tests for residuals, and residual plots
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MTSdiag(model, gof = 24, adj = 0, level = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A fitted multivariate time series model
%%     ~~Describe \code{model} here~~
}
  \item{gof}{The number of lags of residual cross-correlation matrices
%%     ~~Describe \code{gof} here~~
}
  \item{adj}{
The adjustment for degrees of freedom of Ljung-Box statistics. 
Typically, the number of fitted coefficients of the model. Default is zero.
%%     ~~Describe \code{adj} here~~
}
  \item{level}{Logical switch for printing residual cross-correlation matrices
%%     ~~Describe \code{level} here~~
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
\author{Ruey S Tsay
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
function (model, gof = 24, adj = 0, level = F) 
{
    resi = model$residuals
    colnames(resi) <- colnames(model$data)
    ccm(resi, lag = gof, level = level)
    cat("Hit Enter to compute MQ-statistics:", "\n")
    readline()
    mq(resi, lag = gof, adj = adj)
    cat("Hit Enter to obtain residual plots:", "\n")
    readline()
    MTSplot(resi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
