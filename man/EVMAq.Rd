\name{EVMAq}
\alias{EVMAq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exact Likleihood Function for VMA(q) Model
%%  ~~function to do ... ~~
}
\description{Compute the exact log likelihood function of a VMA(q) model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EVMAq(par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{parameter vector
%%     ~~Describe \code{par} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{log likelihood function
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
\note{This is an internal function, not intended to be a general command
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
function (par) 
{
    zt = VMAdata
    k = dim(zt)[2]
    nT = dim(zt)[1]
    q <- MAq
    mu = rep(0, k)
    icnt = 0
    fix <- fix1
    include.mean <- inc.mean
    Ncnt <- VMAecnt
    ist = 0
    if (include.mean) {
        ist = 1
        jdx = c(1:k)[fix[1, ] == 1]
        icnt = length(jdx)
        if (icnt > 0) 
            mu[jdx] = par[1:icnt]
    }
    for (j in 1:k) {
        zt[, j] = zt[, j] - mu[j]
    }
    kq = k * q
    theta = matrix(0, kq, k)
    for (j in 1:k) {
        idx = c(1:kq)[fix[(ist + 1):(ist + kq), j] == 1]
        jcnt = length(idx)
        if (jcnt > 0) {
            theta[idx, j] = par[(icnt + 1):(icnt + jcnt)]
            icnt = icnt + jcnt
        }
    }
    theta = t(theta)
    k1 = dim(theta)[2]
    Theta = theta[, 1:k]
    Zt = zt
    est = NULL
    if (q > 1) {
        z0 = cbind(diag(1, k * (q - 1)), matrix(0, k * (q - 1), 
            k))
        Theta = rbind(theta, z0)
        Zt = cbind(zt, matrix(0, nT, k * (q - 1)))
    }
    m1 = eigen(Theta)
    V1 = m1$values
    M1 = Mod(V1)
    ich = 0
    for (i in 1:k1) {
        if (M1[i] > 1) {
            V1[i] = 1/V1[i]
            ich = 1
        }
    }
    if (ich > 0) {
        P1 = m1$vectors
        P1i = solve(P1)
        Theta = Re(P1 \%*\% diag(V1) \%*\% P1i)
        beta = t(Theta[1:k, ])
        ist = 0
        nr = kq
        if (include.mean) {
            ist = 1
            nr = kq + 1
        }
        for (j in 1:k) {
            idx = c(1:k1)[fix1[(ist + 1):nr, j] == 1]
            if (length(idx) > 0) {
                est = c(est, beta[idx, j])
            }
        }
        par = est
    }
    if (EstStep) {
        if (Ncnt > 0) {
            ParE = c(ParE[1:Ncnt], est)
        }
        else {
            ParE = est
        }
    }
    theta = Theta[1:k, ]
    at = mFilter(zt, theta)
    sig = t(at) \%*\% at/nT
    da = zt
    sig = (sig + t(sig))/2
    m1 = eigen(sig)
    va = m1$values + 10^(-10)
    VA = diag(1/sqrt(va))
    P = m1$vectors
    SigH = P \%*\% VA \%*\% t(P)
    if (q > 1) 
        SigH = kronecker(diag(rep(1, q)), SigH)
    k1 = k * q
    Psi = diag(rep(1, k1))
    Psi = cbind(Theta, Psi)
    tmp = Theta
    X = -SigH
    tmp1 = -SigH \%*\% tmp
    X = rbind(X, tmp1)
    if (nT > 2) {
        for (i in 2:nT) {
            tmp = tmp \%*\% Theta
            Psi = cbind(tmp, Psi)
            tmp1 = -SigH \%*\% tmp
            X = rbind(X, tmp1)
        }
    }
    vZt = c(t(Zt))
    Y = rep(0, k1)
    nPsi = dim(Psi)[2]
    for (it in 1:nT) {
        iend = it * k1
        wk1 = vZt[1:iend]
        wk2 = Psi[, (nPsi - iend + 1):nPsi]
        wk = wk2 \%*\% as.matrix(wk1, iend, 1)
        Y = c(Y, SigH \%*\% wk)
    }
    XpX = t(X) \%*\% X
    XpXinv = solve(XpX)
    XpY = t(X) \%*\% Y
    a0H = XpXinv \%*\% XpY
    resi = Y - X \%*\% a0H
    SSr = sum(resi^2)
    d1 = det(XpX)
    d2 = det(sig)
    llike = 0.5 * (SSr + nT * log(2 * pi * d2) + log(d2))
    llike
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
