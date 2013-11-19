\name{BEKK11}
\alias{BEKK11}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{BEKK Model 
%%  ~~function to do ... ~~
}
\description{Estimation of BEKK(1,1) Model. Only k = 2 or 3 is 
available
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BEKK11(rt, include.mean = T, cond.dist = "normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{A T-by-k data matrix of k-dimensional asset returns
%%     ~~Describe \code{rt} here~~
}
  \item{include.mean}{A logical switch to include a constant vector in the 
mean equation
%%     ~~Describe \code{include.mean} here~~
}
  \item{cond.dist}{Conditional innovation distribution. 
Deafult is Gaussian.
%%     ~~Describe \code{cond.dist} here~~
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
> da=read.table("m-ibmsp-6111.txt",header=T)
> rtn=log(da[,2:3]+1)
> m1=BEKK11(rtn)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rt, include.mean = T, cond.dist = "normal") 
{
    if (!is.matrix(rt)) 
        rt = as.matrix(rt)
    T = dim(rt)[1]
    k = dim(rt)[2]
    if (k > 3) {
        cat("Program Note: Dimension is limited to 3", "\n")
        k = 3
    }
    RTN <<- rt[, 1:k]
    include.M <<- include.mean
    mu = apply(RTN, 2, mean)
    Cov1 <<- cov(RTN)
    if (k == 2) {
        A1 = matrix(c(0.1, 0.02, 0.02, 0.1), k, k)
        B1 = matrix(c(0.8, 0.1, 0.1, 0.8), k, k)
        m1 = chol(Cov1)
        if (include.mean) {
            par = c(mu1 = mu[1], mu2 = mu[2], A011 = m1[1, 1], 
                A021 = m1[1, 2], A022 = m1[2, 2], A11 = A1[1, 
                  1], A21 = A1[2, 1], A12 = A1[1, 2], A22 = A1[2, 
                  2], B11 = B1[1, 1], B21 = B1[2, 1], B12 = B1[1, 
                  2], B22 = B1[2, 2])
            S = 1e-06
            S1 = -0.5
            c1 = c(mu1 = -10 * abs(mu[1]), mu2 = -10 * abs(mu[2]), 
                A011 = m1[1, 1] * 0.2, A021 = m1[1, 2] * 0.2, 
                A022 = m1[2, 2] * 0.2, A11 = S, A21 = S1, A12 = S1, 
                A22 = S, B11 = S, B21 = S1, B12 = S1, B22 = S)
            c2 = c(mu1 = 10 * abs(mu[1]), mu2 = 10 * abs(mu[2]), 
                A011 = m1[1, 1] * 1.1, A021 = m1[1, 2] * 1.1, 
                A022 = m1[2, 2] * 1.1, A11 = 1 - S, A21 = -S1, 
                A12 = -S1, A22 = 1 - S, B11 = 1 - S, B21 = -S1, 
                B12 = -S1, B22 = 1 - S)
        }
        else {
            par = c(A011 = m1[1, 1], A021 = m1[1, 2], A022 = m1[2, 
                2], A11 = A1[1, 1], A21 = A1[2, 1], A12 = A1[1, 
                2], A22 = A1[2, 2], B11 = B1[1, 1], B21 = B1[2, 
                1], B12 = B1[1, 2], B22 = B1[2, 2])
            S = 1e-06
            S1 = -0.5
            c1 = c(A011 = m1[1, 1] * 0.2, A021 = m1[1, 2] * 0.2, 
                A022 = m1[2, 2] * 0.2, A11 = S, A21 = S1, A12 = S1, 
                A22 = S, B11 = S, B21 = S1, B12 = S1, B22 = S)
            c2 = c(A011 = m1[1, 1] * 1.1, A021 = m1[1, 2] * 1.1, 
                A022 = m1[2, 2] * 1.1, A11 = 1 - S, A21 = -S1, 
                A12 = -S1, A22 = 1 - S, B11 = 1 - S, B21 = -S1, 
                B12 = -S1, B22 = 1 - S)
        }
    }
    if (k == 3) {
        A1 = matrix(c(0.1, 0.02, 0.02, 0.02, 0.1, 0.02, 0.02, 
            0.02, 0.1), k, k)
        B1 = matrix(c(0.8, 0.02, 0.02, 0.02, 0.8, 0.02, 0.02, 
            0.02, 0.8), k, k)
        m1 = chol(Cov1)
        if (include.mean) {
            par = c(mu1 = mu[1], mu2 = mu[2], mu3 = mu[3], A011 = m1[1, 
                1], A021 = m1[1, 2], A031 = m1[1, 3], A022 = m1[2, 
                2], A032 = m1[2, 3], A033 = m1[3, 3], A11 = A1[1, 
                1], A21 = A1[2, 1], A31 = A1[3, 1], A12 = A1[1, 
                2], A22 = A1[2, 2], A32 = A1[3, 2], A13 = A1[1, 
                3], A23 = A1[2, 3], A33 = A1[3, 3], B11 = B1[1, 
                1], B21 = B1[2, 1], B31 = B1[3, 1], B12 = B1[1, 
                2], B22 = B1[2, 2], B32 = B1[3, 2], B13 = B1[1, 
                3], B23 = B1[2, 3], B33 = B1[3, 3])
            S = 1e-06
            S1 = -0.5
            c1 = c(mu1 = -10 * abs(mu[1]), mu2 = -10 * abs(mu[2]), 
                mu3 = -10 * abs(mu[3]), A011 = m1[1, 1] * 0.2, 
                A021 = m1[1, 2] * 0.2, A031 = m1[1, 3] * 0.2, 
                A022 = m1[2, 2] * 0.2, A032 = m1[2, 3] * 0.2, 
                A033 = m1[3, 3] * 0.2, A11 = S, A21 = S1, A31 = S1, 
                A12 = S1, A22 = S, A32 = S1, A13 = S1, A23 = S1, 
                A33 = S, B11 = S, B21 = S1, B31 = S1, B12 = S1, 
                B22 = S, B32 = S1, B13 = S1, B23 = S1, B33 = S)
            c2 = c(mu1 = 10 * abs(mu[1]), mu2 = 10 * abs(mu[2]), 
                mu3 = 10 * abs(mu[3]), A011 = m1[1, 1] * 1.1, 
                A021 = m1[1, 2] * 1.1, A031 = m1[1, 3] * 1.1, 
                A022 = m1[2, 2] * 1.1, A032 = m1[2, 3] * 1.1, 
                A033 = m1[3, 3] * 1.1, A11 = 1 - S, A21 = -S1, 
                A31 = -S1, A12 = -S1, A22 = 1 - S, A32 = -S1, 
                A13 = -S1, A23 = -S1, A33 = 1 - S, B11 = 1 - 
                  S, B21 = -S1, B31 = -S1, B12 = -S1, B22 = 1 - 
                  S, B32 = -S1, B13 = -S1, B23 = -S1, B33 = 1 - 
                  S)
        }
        else {
            par = c(A011 = m1[1, 1], A021 = m1[1, 2], A031 = m1[1, 
                3], A022 = m1[2, 2], A032 = m1[2, 3], A033 = m1[3, 
                3], A11 = A1[1, 1], A21 = A1[2, 1], A31 = A1[3, 
                1], A12 = A1[1, 2], A22 = A1[2, 2], A32 = A1[3, 
                2], A13 = A1[1, 3], A23 = A1[2, 3], A33 = A1[3, 
                3], B11 = B1[1, 1], B21 = B1[2, 1], B31 = B1[3, 
                1], B12 = B1[1, 2], B22 = B1[2, 2], B32 = B1[3, 
                2], B13 = B1[1, 3], B23 = B1[2, 3], B33 = B1[3, 
                3])
            S = 1e-06
            S1 = -0.5
            c1 = c(A011 = m1[1, 1] * 0.2, A021 = m1[1, 2] * 0.2, 
                A031 = m1[1, 3] * 0.2, A022 = m1[2, 2] * 0.2, 
                A032 = m1[2, 3] * 0.2, A033 = m1[3, 3] * 0.2, 
                A11 = S, A21 = S1, A31 = S1, A12 = S1, A22 = S, 
                A32 = S1, A13 = S1, A23 = S1, A33 = S, B11 = S, 
                B21 = S1, B31 = S1, B12 = S1, B22 = S, B32 = S1, 
                B13 = S1, B23 = S1, B33 = S)
            c2 = c(A011 = m1[1, 1] * 1.1, A021 = m1[1, 2] * 1.1, 
                A031 = m1[1, 3] * 1.1, A022 = m1[2, 2] * 1.1, 
                A032 = m1[2, 3] * 1.1, A033 = m1[3, 3] * 1.1, 
                A11 = 1 - S, A21 = -S1, A31 = -S1, A12 = -S1, 
                A22 = 1 - S, A32 = -S1, A13 = -S1, A23 = -S1, 
                A33 = 1 - S, B11 = 1 - S, B21 = -S1, B31 = -S1, 
                B12 = -S1, B22 = 1 - S, B32 = -S1, B13 = -S1, 
                B23 = -S1, B33 = 1 - S)
        }
    }
    cat("Initial estimates: ", par, "\n")
    cat("Lower limits: ", c1, "\n")
    cat("Upper limits: ", c2, "\n")
    fit = nlminb(start = par, objective = mlikeG, lower = c1, 
        upper = c2)
    epsilon = 3e-04 * fit$par
    if (k == 3) 
        epsilon = 5e-04 * fit$par
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
            Hessian[i, j] = (mlikeG(x1) - mlikeG(x2) - mlikeG(x3) + 
                mlikeG(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    est = fit$par
    se.coef = sqrt(diag(solve(Hessian)))
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", 
        " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
    m2 = BEKK11vol(RTN, est, include.mean = include.M)
    Sigma.t = m2$volmtx
    BEKK11 <- list(estimates = est, HessianMtx = Hessian, Sigma.t = Sigma.t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
