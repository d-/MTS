\name{VAR}
\alias{VAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vector Autoregressive Model
%%  ~~function to do ... ~~
}
\description{Perform least squares estimation of a VAR model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VAR(x, p = 1, output = T, include.mean = T, fixed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k matrix of k-dimensional time series
%%     ~~Describe \code{x} here~~
}
  \item{p}{Order of VAR model. Default is 1.
%%     ~~Describe \code{p} here~~
}
  \item{output}{A logical switch to control output. Default is with output.
%%     ~~Describe \code{output} here~~
}
  \item{include.mean}{A logical switch. It is true if mean vector is estimated.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used in constrained estimation. 
It is used mainly in model simplifcation, e.g., removing insignificant 
estimates.
%%     ~~Describe \code{fixed} here~~
}
}
\details{To remove insignificant estimates, one specifies a threshold 
for individual t-ratio. The fixed matrix is then defined automatically to 
identify those parameters for removal. 
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

\seealso{refVAR command
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, p = 1, output = T, include.mean = T, fixed = NULL) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    Tn = dim(x)[1]
    k = dim(x)[2]
    if (p < 1) 
        p = 1
    idm = k * p
    ne = Tn - p
    ist = p + 1
    y = x[ist:Tn, ]
    if (include.mean) {
        idm = idm + 1
        xmtx = cbind(rep(1, ne), x[p:(Tn - 1), ])
    }
    else {
        xmtx = x[p:(Tn - 1), ]
    }
    if (p > 1) {
        for (i in 2:p) {
            xmtx = cbind(xmtx, x[(ist - i):(Tn - i), ])
        }
    }
    ndim = ncol(xmtx)
    if (length(fixed) == 0) {
        paridx = matrix(1, ndim, k)
    }
    else {
        paridx = fixed
    }
    res = NULL
    beta = matrix(0, ndim, k)
    sdbeta = matrix(0, ndim, k)
    npar = 0
    for (i in 1:k) {
        idx = c(1:ndim)[paridx[, i] == 1]
        resi = y[, i]
        if (length(idx) > 0) {
            xm = as.matrix(xmtx[, idx])
            npar = npar + dim(xm)[2]
            xpx = t(xm) \%*\% xm
            xpxinv = solve(xpx)
            xpy = t(xm) \%*\% as.matrix(y[, i], ne, 1)
            betai = xpxinv \%*\% xpy
            beta[idx, i] = betai
            resi = y[, i] - xm \%*\% betai
            nee = dim(xm)[2]
            sse = sum(resi * resi)/(Tn - p - nee)
            dd = diag(xpxinv)
            sdbeta[idx, i] = sqrt(dd * sse)
        }
        res = cbind(res, resi)
    }
    sse = t(res) \%*\% res/(Tn - p)
    aic = 0
    bic = 0
    hq = 0
    Phi = NULL
    Ph0 = NULL
    jst = 0
    if (include.mean) {
        Ph0 = beta[1, ]
        se = sdbeta[1, ]
        if (output) {
            cat("Constant term:", "\n")
            cat("Estimates: ", Ph0, "\n")
            cat("Std.Error: ", se, "\n")
        }
        jst = 1
    }
    if (include.mean) {
        for (i in 1:k) {
            if (abs(Ph0[i]) > 1e-08) 
                npar = npar - 1
        }
    }
    if (output) 
        cat("AR coefficient matrix", "\n")
    for (i in 1:p) {
        phi = t(beta[(jst + 1):(jst + k), ])
        se = t(sdbeta[(jst + 1):(jst + k), ])
        if (output) {
            cat("AR(", i, ")-matrix", "\n")
            print(phi, digits = 3)
            cat("standard error", "\n")
            print(se, digits = 3)
        }
        jst = jst + k
        Phi = cbind(Phi, phi)
    }
    if (output) {
        cat(" ", "\n")
        cat("Residuals cov-mtx:", "\n")
        print(sse)
        cat(" ", "\n")
    }
    dd = det(sse)
    d1 = log(dd)
    aic = d1 + (2 * npar)/Tn
    bic = d1 + log(Tn) * npar/Tn
    hq = d1 + 2 * log(log(Tn)) * npar/Tn
    if (output) {
        cat("det(SSE) = ", dd, "\n")
        cat("AIC = ", aic, "\n")
        cat("BIC = ", bic, "\n")
        cat("HQ  = ", hq, "\n")
    }
    VAR <- list(data = x, cnst = include.mean, order = p, coef = beta, 
        aic = aic, bic = bic, hq = hq, residuals = res, secoef = sdbeta, 
        Sigma = sse, Phi = Phi, Ph0 = Ph0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
