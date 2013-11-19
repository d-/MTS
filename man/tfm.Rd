\name{tfm}
\alias{tfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transfer Function Model
%%  ~~function to do ... ~~
}
\description{Estimates a transform function model. This program does not allow 
rational transfer function model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tfm(y, x, b = 0, s = 1, p = 0, q = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Data vector of dependent (output) variable
%%     ~~Describe \code{y} here~~
}
  \item{x}{Data vector of independent variable
%%     ~~Describe \code{x} here~~
}
  \item{b}{deadtime or delay
%%     ~~Describe \code{b} here~~
}
  \item{s}{The order of the transfer function polynomial
%%     ~~Describe \code{s} here~~
}
  \item{p}{AR order of the disturbance
%%     ~~Describe \code{p} here~~
}
  \item{q}{MA order of the disturbance
%%     ~~Describe \code{q} here~~
}
}
\details{The model entertained is y_t = c_0+ v(B)x_t + n_t. 
v(B) = 1- v1*B - ... - vs*B^s, and n_t is an ARMA(p,q) process.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{coef }{Coefficient estimates of the transfer function}
\item{se.coef }{Standard errors of the transfer function coefficients}
\item{coef.arma }{Coefficient estimates of ARMA models}
\item{se.arma }{Standard errors of ARMA coefficients}
\item{nt }{The disturbance series}
\item{residuals }{The residual series}
%% ...
}
\references{Box and Jenkins (1970)
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
function (y, x, b = 0, s = 1, p = 0, q = 0) 
{
    nT = length(y)
    T1 = length(x)
    nT = min(nT, T1)
    mx = b + s
    ist = mx + 1
    y1 = y[ist:nT]
    X = x[(s + 1):(nT - b)]
    if (s > 0) {
        for (i in 1:s) {
            X = cbind(X, x[(s + 1 - i):(nT - b - i)])
        }
    }
    nx = ncol(X)
    m1 = arima(y1, order = c(p, 0, q), xreg = X)
    se = sqrt(diag(m1$var.coef))
    coef.arma = NULL
    se.arma = NULL
    pq = p + q
    if (pq > 0) {
        coef.arma = m1$coef[1:pq]
        se.arma = se[1:pq]
        p1 = cbind(coef.arma, se.arma)
        cat("ARMA coefficients & s.e.:", "\n")
        print(t(p1), digits = 3)
    }
    v = m1$coef[(pq + 1):(pq + 1 + nx)]
    se.v = se[(pq + 1):(pq + 1 + nx)]
    pr = cbind(v, se.v)
    cat("Transfer function coefficients & s.e.:", "\n")
    print(t(pr), digits = 3)
    res = m1$residuals
    beta = matrix(v[2:(nx + 1)], nx, 1)
    nt = y1 - v[1] - X \%*\% beta
    tfm <- list(coef = v, se.coef = se.v, coef.arma = coef.arma, 
        se.arma = se.arma, nt = nt, residuals = res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
