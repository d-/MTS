\name{ECMvar}
\alias{ECMvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Error-Correction VAR Models
%%  ~~function to do ... ~~
}
\description{Performs estimation of an Error-Correction VAR(p) model using 
Quasi Maximum Likelihood Method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ECMvar(x, p, ibeta, include.const = FALSE, fixed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional co-integrated VAR process
%%     ~~Describe \code{x} here~~
}
  \item{p}{VAR order
%%     ~~Describe \code{p} here~~
}
  \item{ibeta}{Initial estimate of the co-integrating matrix.
%%     ~~Describe \code{ibeta} here~~
}
  \item{include.const}{A logical switch to include a constant term in the model. The default is no constant
%%     ~~Describe \code{include.const} here~~
}
  \item{fixed}{A logical matrix to set zero parameter constraints.
%%     ~~Describe \code{fixed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of estimation of an ECM VAR model, including 
alpha and beta matrices of co-integrating information.
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-bnd.txt")
bnd=da[,4:5]
beta=matrix(c(1,-0.886),2,1)
m1=ECMvar(bnd,3,ibeta=beta)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, p, ibeta, include.const = FALSE, fixed = NULL) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    if (!is.matrix(ibeta)) 
        ibeta = as.matrix(ibeta)
    if (p < 1) 
        p = 1
    m = dim(ibeta)[2]
    cat("Order p: ", p, " Co-integrating rank: ", m, "\n")
    nT = dim(x)[1]
    k = dim(x)[2]
    dx = x[2:nT, ] - x[1:(nT - 1), ]
    dx = rbind(rep(0, k), dx)
    wt = x \%*\% ibeta
    m1 = ECMvar1(x, p, wt, include.const = include.const, output = FALSE)
    est = m1$coef
    se.est = m1$secoef
    alpha = t(est[1:m, ])
    se.alpha = t(se.est[1:m, ])
    icnt = m
    idm = dim(est)[1]
    if (idm > icnt) {
        phip = est[(icnt + 1):idm, ]
        se.phip = se.est[(icnt + 1):idm, ]
    }
    par = c(alpha)
    separ = c(se.alpha)
    par = c(par, c(ibeta[(m + 1):k, ]))
    separ = c(separ, rep(1/sqrt(nT), (k - m) * m))
    par = c(par, c(phip))
    separ = c(separ, c(se.phip))
    npar = length(par)
    ECMx <<- x
    Ncoint <<- m
    ECMp <<- p
    ECMfix <<- fixed
    ECMconst <<- include.const
    m2 = ECMxmtx(x, p, m, include.const)
    ECMy <<- m2$y
    ECMxm <<- m2$xm
    cat("Number of parameters: ", length(par), "\n")
    cat("initial estimates: ", par, "\n")
    lowerBounds = par
    upperBounds = par
    mult = 1.5
    for (j in 1:npar) {
        lowerBounds[j] = par[j] - mult * separ[j]
        upperBounds[j] = par[j] + mult * separ[j]
    }
    cat("Par. Lower-bounds: ", lowerBounds, "\n")
    cat("Par. Upper-bounds: ", upperBounds, "\n")
    details = FALSE
    if (details) {
        fit = nlminb(start = par, objective = LECMvar, lower = lowerBounds, 
            upper = upperBounds, control = list(trace = 3))
    }
    else {
        fit = nlminb(start = par, objective = LECMvar, lower = lowerBounds, 
            upper = upperBounds)
    }
    epsilon = 1e-04 * fit$par
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
            Hessian[i, j] = (LECMvar(x1) - LECMvar(x2) - LECMvar(x3) + 
                LECMvar(x4))/(4 * epsilon[i] * epsilon[j])
        }
    }
    est = fit$par
    cat("Final   Estimates: ", est, "\n")
    se.coef = sqrt(diag(solve(Hessian)))
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", 
        " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
    km = k * m
    kmm = (k - m) * m
    alpha = matrix(est[1:km], k, m)
    se.alpha = matrix(se.coef[1:km], k, m)
    beta1 = matrix(est[(km + 1):(km + kmm)], (k - m), m)
    se.beta1 = matrix(se.coef[(km + 1):(km + kmm)], (k - m), 
        m)
    Im = diag(rep(1, m))
    beta = rbind(Im, beta1)
    se.beta = rbind(Im, se.beta1)
    icnt = km + kmm
    idm = k * (p - 1)
    if (include.const) 
        idm = idm + 1
    if (icnt < npar) {
        phip = matrix(est[(icnt + 1):npar], idm, k)
        se.phip = matrix(se.coef[(icnt + 1):npar], idm, k)
    }
    cat("alpha: ", "\n")
    print(alpha, digits = 3)
    cat("standard error", "\n")
    print(se.alpha, digits = 3)
    cat("beta: ", "\n")
    print(beta, digits = 3)
    cat("standard error", "\n")
    print(se.beta, digits = 3)
    icnt = 0
    if (include.const) {
        cat("constant term:", "\n")
        print(phip[1, ], digits = 3)
        se = se.phip[1, ]
        cat("standard error", "\n")
        print(se, digits = 3)
        icnt = 1
    }
    if (p > 1) {
        cat("AR coefficient matrix", "\n")
        jst = icnt
        for (i in 1:(p - 1)) {
            cat("AR(", i, ")-matrix", "\n")
            phi = t(phip[(jst + 1):(jst + k), ])
            se = t(se.phip[(jst + 1):(jst + k), ])
            print(phi, digits = 3)
            cat("standard error", "\n")
            print(se, digits = 3)
            jst = jst + k
        }
    }
    Pi = alpha \%*\% t(beta)
    Ome = rbind(t(Pi), phip)
    resi = ECMy - as.matrix(ECMxm) \%*\% Ome
    sse = t(resi) \%*\% resi/nT
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
    ECMvar <- list(data = x, ncoint = m, arorder = p, include.const = include.const, 
        alpha = alpha, se.alpha = se.alpha, beta = beta, se.beta = se.beta, 
        aic = aic, bic = bic, residuals = resi, Phip = phip, 
        se.Phip = se.phip, Sigma = sse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
