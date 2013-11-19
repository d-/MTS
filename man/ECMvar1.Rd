\name{ECMvar1}
\alias{ECMvar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Error-Correction VAR Model 1
%%  ~~function to do ... ~~
}
\description{Perform least-squares estimation of an ECM VAR(p) model with 
known co-integrating processes
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ECMvar1(x, p, wt, include.const = FALSE, fixed = NULL, output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional co-integrated VAR process
%%     ~~Describe \code{x} here~~
}
  \item{p}{VAR order
%%     ~~Describe \code{p} here~~
}
  \item{wt}{A T-by-m data matrix of m-dimensional co-integrated process
%%     ~~Describe \code{wt} here~~
}
  \item{include.const}{A logical switch to include a constant term. 
Default is no constant.
%%     ~~Describe \code{include.const} here~~
}
  \item{fixed}{A logical matrix to set zero parameter constraints
%%     ~~Describe \code{fixed} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of least-squares estimation of an ECM mdoel for VAR process.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 5)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ECMvar
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-bnd.txt")
bdn=da[,4:5]
wt=bnd[,1]-0.886*bnd[,2]
m1=ECMvar1(bnd,3,wt)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, p, wt, include.const = FALSE, fixed = NULL, output = TRUE) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    if (p < 1) 
        p = 1
    nT = dim(x)[1]
    k = dim(x)[2]
    dx = x[2:nT, ] - x[1:(nT - 1), ]
    dx = rbind(rep(0, k), dx)
    wt = as.matrix(wt)
    m = dim(wt)[2]
    wtadj = wt
    idm = k * (p - 1) + m
    if (include.const) {
        idm = idm + 1
    }
    else {
        wtadj = wt - matrix(1, nT, 1) \%*\% matrix(apply(wt, 2, 
            mean), 1, m)
    }
    ist = max(1, p)
    ne = nT - ist + 1
    y = dx[ist:nT, ]
    xmtx = wtadj[(ist - 1):(nT - 1), ]
    if (include.const) 
        xmtx = cbind(xmtx, rep(1, (nT - ist + 1)))
    if (p > 1) {
        for (i in 2:p) {
            ii = i - 1
            xmtx = cbind(xmtx, dx[(ist - ii):(nT - ii), ])
        }
    }
    y = as.matrix(y)
    xmtx = as.matrix(xmtx)
    sdbeta = matrix(1, idm, k)
    if (length(fixed) < 1) {
        xpx = t(xmtx) \%*\% xmtx
        xpxinv = solve(xpx)
        xpy = t(xmtx) \%*\% y
        beta = xpxinv \%*\% xpy
        yhat = xmtx \%*\% beta
        resi = y - yhat
        sse = (t(resi) \%*\% resi)/ne
        dd = diag(xpxinv)
        sdbeta = matrix(1, idm, k)
        for (i in 1:k) {
            sdbeta[, i] = sqrt(sse[i, i] * dd)
        }
        npar = idm * k
    }
    else {
        resi = NULL
        beta = matrix(0, idm, k)
        sdbeta = matrix(1, idm, k)
        npar = 0
        for (i in 1:k) {
            idx = c(1:idm)[fixed[, i] == 1]
            npi = length(idx)
            cat("Equation: ", i, " npar = ", npi, "\n")
            npar = npar + npi
            if (npi > 0) {
                xm = xmtx[, idx]
                xpx = t(xm) \%*\% xm
                xpxinv = solve(xpx)
                dd = diag(xpxinv)
                xpy = t(xm) \%*\% y[, i]
                betai = xpxinv \%*\% xpy
                beta[idx, i] = betai
                res = y[, i] - xm \%*\% betai
                resi = cbind(resi, res)
                se2 = sum(res^2)/ne
                sdbeta[idx, i] = sqrt(dd * se2)
            }
        }
    }
    sse = (t(resi) \%*\% resi)/ne
    aic = 0
    bic = 0
    if (output) {
        alpha = beta[1:m, ]
        icnt = m
        if (include.const) {
            icnt = m + 1
            c = beta[icnt, ]
        }
        se = sdbeta[1:m, ]
        cat("alpha: ", "\n")
        print(t(alpha), digits = 3)
        cat("standard error", "\n")
        print(t(se), digits = 3)
        if (include.const) {
            cat("constant term:", "\n")
            print(c, digits = 3)
            se = sdbeta[icnt, ]
            cat("standard error", "\n")
            print(se, digits = 3)
        }
        if (p > 1) {
            cat("AR coefficient matrix", "\n")
            jst = icnt
            for (i in 1:(p - 1)) {
                cat("AR(", i, ")-matrix", "\n")
                phi = t(beta[(jst + 1):(jst + k), ])
                se = t(sdbeta[(jst + 1):(jst + k), ])
                print(phi, digits = 3)
                cat("standard error", "\n")
                print(se, digits = 3)
                jst = jst + k
            }
        }
        cat("-----", "\n")
        cat("Residuals cov-mtx:", "\n")
        print(sse)
        cat("      ", "\n")
        dd = det(sse)
        cat("det(sse) = ", dd, "\n")
        d1 = log(dd)
        aic = d1 + (2 * npar)/nT
        bic = d1 + log(nT) * npar/nT
        cat("AIC = ", aic, "\n")
        cat("BIC = ", bic, "\n")
    }
    ECMvar1 <- list(data = x, wt = wt, arorder = p, include.const = include.const, 
        coef = beta, aic = aic, bic = bic, residuals = resi, 
        secoef = sdbeta, Sigma = sse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
