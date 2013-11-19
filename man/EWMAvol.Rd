\name{EWMAvol}
\alias{EWMAvol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponentially Weighted Moving-Average for Volatility
%%  ~~function to do ... ~~
}
\description{Use exponentially weighted moving-average method to 
compute the volatility matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EWMAvol(rtn, lambda = 0.96)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns, assuming the 
mean is zero
%%     ~~Describe \code{rtn} here~~
}
  \item{lambda}{Smoothing parameter. The deafult is 0.96. If lambda is negative, then multivariate Gaussian likelihood is used to estimate the smoothing parameter.
%%     ~~Describe \code{lambda} here~~
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
\references{Tsay (2014, Chapter 7)
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
da=read.table("m-dec125910-6111.txt",header=T)
rtn=log(da[,2:4]+1)
m1=VAR(rtn,1)
at=m1$residuals
m2=EWMAvol(at,lambda=-0.1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rtn, lambda = 0.96) 
{
    if (!is.matrix(rtn)) 
        rtn = as.matrix(rtn)
    T = dim(rtn)[1]
    k = dim(rtn)[2]
    x = scale(rtn, center = TRUE, scale = FALSE)
    XX <<- x
    if (lambda <= 0) {
        par = c(lambda = 0.96)
        S = 10^{
            -5
        }
        lowerb = c(lambda = S)
        upperb = c(lambda = 1 - S)
        fit = nlminb(start = par, objective = MGAUS, lower = lowerb, 
            upper = upperb)
        epsilon = 1e-04 * fit$par
        npar = length(par)
        Hessian = matrix(0, ncol = npar, nrow = npar)
        for (i in 1:npar) {
            for (j in 1:npar) {
                x1 = x2 = x3 = x4 = fit$par
                x1[i] = x1[i] + epsilon[i]
                x1[j] = x1[j] + epsilon[j]
                x2[i] = x2[i] + epsilon[i]
                x2[j] = x2[j] - epsilon[j]
                x3[i] = x3[i] - epsilon[i]
                x3[j] = x3[j] + epsilon[j]
                x4[i] = x4[i] - epsilon[i]
                x4[j] = x4[j] - epsilon[j]
                Hessian[i, j] = (MGAUS(x1) - MGAUS(x2) - MGAUS(x3) + 
                  MGAUS(x4))/(4 * epsilon[i] * epsilon[j])
            }
        }
        se.coef = sqrt(diag(solve(Hessian)))
        tval = fit$par/se.coef
        matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
        dimnames(matcoef) = list(names(tval), c(" Estimate", 
            " Std. Error", " t value", "Pr(>|t|)"))
        cat("\nCoefficient(s):\n")
        printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
        lambda = fit$par
    }
    if (lambda > 0) {
        h1 = 1 - lambda
        Sigt = cov(XX)
        V1 = c(Sigt)
        for (t in 2:T) {
            xx = as.numeric(x[t - 1, ])
            for (i in 1:k) {
                Sigt[i, ] = h1 * xx * xx[i] + lambda * Sigt[i, 
                  ]
            }
            V1 = rbind(V1, c(Sigt))
        }
    }
    EWMAvol <- list(Sigma.t = V1, return = rtn, lambda = lambda)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
