\name{VARMApred}
\alias{VARMApred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VARMA Prediction
%%  ~~function to do ... ~~
}
\description{Compute forecasts and their associate forecast error covariances 
of a VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMApred(model, h = 1, orig = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A fitted VARMA model
%%     ~~Describe \code{model} here~~
}
  \item{h}{Number of steps of frecasts, i.e., forecast horizon.
%%     ~~Describe \code{h} here~~
}
  \item{orig}{Forecast origin. Default is end of the sample.
%%     ~~Describe \code{orig} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Forecasts and their standard errors
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 3)
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
m2=VARMA(zt,p=3,q=1) # where zt is the growth rates of PCE and DSPI.
m2p=VARMApred(m2,h=6)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, h = 1, orig = 0) 
{
    x = as.matrix(model$data)
    resi = as.matrix(model$residuals)
    sig = model$Sigma
    Phi = model$Phi
    Theta = model$Theta
    Ph0 = model$Ph0
    p = model$ARorder
    q = model$MAorder
    if (p < 0) 
        p = 0
    if (q < 0) 
        q = 0
    if (h < 1) 
        h = 1
    nT = dim(x)[1]
    k = dim(x)[2]
    T1 = dim(resi)[1]
    if (nT > T1) {
        r1 = matrix(0, (nT - T1), k)
        resi = rbind(r1, resi)
    }
    if (length(Ph0) < 1) 
        Ph0 = rep(0, k)
    if (orig < 1) 
        orig = nT
    if (orig > T) 
        orig = nT
    px = x[1:orig, ]
    presi = resi[1:orig, ]
    psi = diag(rep(1, k))
    wk = c(psi)
    lag = max(1, h)
    for (i in 1:lag) {
        if (i <= p) {
            idx = (i - 1) * k
            tmp = Phi[, (idx + 1):(idx + k)]
        }
        else {
            tmp = matrix(0, k, k)
        }
        if (i <= q) {
            mdx = (i - 1) * k
            tmp = tmp - Theta[, (mdx + 1):(mdx + k)]
        }
        jj = i - 1
        jp = min(jj, p)
        if (jp > 0) {
            for (j in 1:jp) {
                jdx = (j - 1) * k
                idx = (i - j) * k
                w1 = Phi[, (jdx + 1):(jdx + k)]
                w2 = psi[, (idx + 1):(idx + k)]
                tmp = tmp + w1 \%*\% w2
            }
        }
        psi = cbind(psi, tmp)
        wk = cbind(wk, c(tmp))
    }
    sefcst = NULL
    for (j in 1:h) {
        fcst = Ph0
        Sig = sig
        t = orig + j
        if (p > 0) {
            for (ii in 1:p) {
                idx = (ii - 1) * k
                ph = Phi[, (idx + 1):(idx + k)]
                fcst = fcst + matrix(px[(t - ii), ], 1, k) \%*\% 
                  t(ph)
            }
        }
        if (q > 0) {
            for (jj in 1:q) {
                idx = (jj - 1) * k
                if ((t - jj) <= orig) {
                  th = Theta[, (idx + 1):(idx + k)]
                  fcst = fcst - matrix(resi[(t - jj), ], 1, k) \%*\% 
                    t(th)
                }
            }
        }
        px = rbind(px, fcst)
        if (j > 1) {
            Sig = sig
            for (jj in 2:j) {
                jdx = (jj - 1) * k
                wk = psi[, (jdx + 1):(jdx + k)]
                Sig = Sig + wk \%*\% sig \%*\% t(wk)
            }
        }
        sefcst = rbind(sefcst, sqrt(diag(Sig)))
    }
    cat("Predictions at origin ", orig, "\n")
    print(px[(orig + 1):(orig + h), ], digits = 4)
    cat("Standard errors of predictions", "\n")
    if (h == 1) {
        print(sefcst, digits = 4)
    }
    else {
        print(sefcst[1:h, ], digits = 4)
    }
    if (orig < nT) {
        cat("Observations, predictions, and errors: ", "\n")
        tmp = NULL
        jend = min(nT, orig + h)
        for (t in (orig + 1):jend) {
            case = c(t, x[t, ], px[t, ], x[t, ] - px[t, ])
            tmp = rbind(tmp, case)
        }
        colnames(tmp) <- c("time", rep("obs", k), rep("fcst", 
            k), rep("err", k))
        idx = c(1)
        for (j in 1:k) {
            idx = c(idx, c(0, 1, 2) * k + j + 1)
        }
        tmp = tmp[, idx]
        print(tmp, digits = 4)
    }
    VARMApred <- list(pred = px[(orig + 1):(orig + h), ], se.err = sefcst, 
        orig = orig)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
