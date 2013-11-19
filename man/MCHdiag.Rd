\name{MCHdiag}
\alias{MCHdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Conditional Heteroscadastic Model Cheking
%%  ~~function to do ... ~~
}
\description{Apply four portmanteau test statistics to check the 
validity of a fitted multivariate volatility model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MCHdiag(at, Sigma.t, m = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{at}{A T-by-k matrix of residuals for a k-dimensional asset returns
%%     ~~Describe \code{at} here~~
}
  \item{Sigma.t}{The fitted volatility matrices. The dimension is 
T-by-k^2 matrix
%%     ~~Describe \code{Sigma.t} here~~
}
  \item{m}{The number of lags used in the tests. Default is 10.
%%     ~~Describe \code{m} here~~
}
}
\details{The four test statistics are given in Tsay (2014, Chapter 7)
%%  ~~ If necessary, more details than the description above ~~
}
\value{Test statistics and their p-values
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
da=read.table("m-ibmsp-6111.txt",header=T)
rtn=log(da[,2:3]+1)
m1a=BEKK11(rtn)
names(m1a)
Sigma.t=m1a$Sigma.t
at=cbind(rtn[,1]-0.00776,rtn[,2]-0.00565)
MCHdiag(at,Sigma.t)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (at, Sigma.t, m = 10) 
{
    if (!is.matrix(at)) 
        at = as.matrix(at)
    if (!is.matrix(Sigma.t)) 
        Sigma.t = as.matrix(Sigma.t)
    nT = dim(at)[1]
    k = dim(at)[2]
    nT1 = dim(Sigma.t)[1]
    k1 = dim(Sigma.t)[2]
    if ((nT != nT1) || (k1 != k^2)) {
        cat("Inconsistency in dimensions", "\n")
    }
    else {
        et = NULL
        etat = NULL
        for (i in 1:nT) {
            Vt = matrix(Sigma.t[i, ], k, k)
            Vtinv = solve(Vt)
            x = matrix(at[i, ], 1, k)
            tmp = x \%*\% Vtinv \%*\% t(x) - k
            et = c(et, tmp)
            m1 = eigen(Vt)
            P = m1$vectors
            lam = m1$values
            d1 = diag(1/sqrt(lam))
            Vthalf = P \%*\% d1 \%*\% t(P)
            wk = x \%*\% Vthalf
            etat = rbind(etat, wk)
        }
        m1 = acf(et, m, plot = FALSE)
        acf = m1$acf[2:(m + 1)]
        tmp = acf^2/c(rep(nT, m) - c(1:m))
        Q1 = sum(tmp) * nT * (nT + 2)
        pv1 = 1 - pchisq(Q1, m)
        lag = m
        mu = -(rep(nT, lag) - c(1:lag))/(nT * (nT - 1))
        v1 = rep(5 * nT^4, lag) - (5 * c(1:lag) + 9) * nT^3 + 
            9 * (c(1:lag) - 2) * nT^2 + 2 * c(1:lag) * (5 * c(1:lag) + 
            8) * nT + 16 * c(1:lag)^2
        v1 = v1/(5 * (nT - 1)^2 * nT^2 * (nT + 1))
        ret = rank(et)
        m2 = acf(ret, m, plot = FALSE)
        acf = m2$acf[2:(m + 1)]
        Qr = sum((acf - mu)^2/v1)
        pv2 = 1 - pchisq(Qr, m)
        x = etat^2
        g0 = var(x)
        ginv = solve(g0)
        qm = 0
        for (i in 1:lag) {
            x1 = x[(i + 1):nT, ]
            x2 = x[1:(nT - i), ]
            g = cov(x1, x2)
            g = g * (nT - i - 1)/(nT - 1)
            h = t(g) \%*\% ginv \%*\% g \%*\% ginv
            qm = qm + nT * nT * sum(diag(h))/(nT - i)
        }
        QKm = qm
        pv3 = 1 - pchisq(QKm, k^2 * m)
        q95 = quantile(et, 0.95)
        idx = c(1:nT)[et <= q95]
        x = etat[idx, ]^2
        eT = length(idx)
        g0 = var(x)
        ginv = solve(g0)
        qm = 0
        for (i in 1:lag) {
            x1 = x[(i + 1):eT, ]
            x2 = x[1:(eT - i), ]
            g = cov(x1, x2)
            g = g * (eT - i - 1)/(eT - 1)
            h = t(g) \%*\% ginv \%*\% g \%*\% ginv
            qm = qm + eT * eT * sum(diag(h))/(eT - i)
        }
        Qrm = qm
        pv4 = 1 - pchisq(Qrm, k^2 * m)
        cat("Test results: ", "\n")
        cat("Q(m) of et:", "\n")
        cat("Test and p-value: ", c(Q1, pv1), "\n")
        cat("Rank-based test:", "\n")
        cat("Test and p-value: ", c(Qr, pv2), "\n")
        cat("Qk(m) of epsilon_t:", "\n")
        cat("Test and p-value: ", c(QKm, pv3), "\n")
        cat("Robust Qk(m): ", "\n")
        cat("Test and p-value: ", c(Qrm, pv4), "\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
