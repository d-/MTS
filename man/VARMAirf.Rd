\name{VARMAirf}
\alias{VARMAirf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Impulse Response Funcation of a VARMA Model
%%  ~~function to do ... ~~
}
\description{Compute the impulse response function of a given VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMAirf(Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, orth = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficients. 
Phi=[Phi1, Phi2, ..., Phip].
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{A k-by-kq matrix of VMA coefficients. 
Theta=[Theta1, Theta2, ..., Thetaq]
%%     ~~Describe \code{Theta} here~~
}
  \item{Sigma}{Covariance matrix (k-by-k) of the innovations.
%%     ~~Describe \code{Sigma} here~~
}
  \item{lag}{Number of lags of impulse response functions to be computed
%%     ~~Describe \code{lag} here~~
}
  \item{orth}{A logical switch to use orthogonal innovations.
Deafult is to perform orthogonalization of the innovations.
%%     ~~Describe \code{orth} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Plots and values of impulse response functions and cumulative 
impulse response function.
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

\seealso{VARMApsi command
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
th1=matrix(c(-0.5,0.2,0.0,-0.6),2,2)
Sig=matrix(c(4,1,1,1),2,2)
m1=VARMAirf(Phi=p1,Theta=th1,Sigma=Sig)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, orth = TRUE) 
{
    q = 0
    p = 0
    k = 0
    if (length(Theta) > 0) {
        k = dim(Theta)[1]
        k1 = dim(Theta)[2]
        q = floor(k1/k)
    }
    if (length(Phi) > 0) {
        k = dim(Phi)[1]
        k1 = dim(Phi)[2]
        p = floor(k1/k)
    }
    if (is.null(Sigma)) {
        Sigma = diag(rep(1, k))
    }
    if (orth) {
        m1 = eigen(Sigma)
        v1 = sqrt(m1$values)
        vv = diag(v1)
        Pmtx = m1$vectors
        Sh = Pmtx \%*\% vv \%*\% t(Pmtx)
    }
    if (k < 1) 
        k = 1
    PSI = diag(rep(1, k))
    if (orth) {
        WGT = c(PSI \%*\% Sh)
    }
    else {
        WGT = c(PSI)
    }
    for (il in 1:lag) {
        ilk = il * k
        tmp = matrix(0, k, k)
        if ((q > 0) && (il <= q)) 
            tmp = -Theta[, (ilk - k + 1):ilk]
        if (p > 0) {
            iend = min(il, p)
            for (j in 1:iend) {
                jdx = (il - j)
                kdx = j * k
                tmp = tmp + Phi[, (kdx - k + 1):kdx] \%*\% PSI[, 
                  (jdx * k + 1):(jdx * k + k)]
            }
        }
        PSI = cbind(PSI, tmp)
        if (orth) {
            WGT = cbind(WGT, c(tmp \%*\% Sh))
        }
        else {
            WGT = cbind(WGT, c(tmp))
        }
    }
    wk1 = WGT
    for (i in 1:k^2) {
        wk1[i, ] = cumsum(WGT[i, ])
    }
    tdx = c(1:(lag + 1)) - 1
    par(mfcol = c(k, k))
    if (orth) {
        gmax = max(WGT)
        gmin = min(WGT)
        cx = (gmax - gmin)/10
        gmax = gmax + cx
        gmin = gmin - cx
        for (j in 1:k^2) {
            plot(tdx, WGT[j, ], type = "l", xlab = "lag", ylab = "IRF", 
                ylim = c(gmin, gmax))
            points(tdx, WGT[j, ], pch = "*")
            title(main = "Orth. innovations")
        }
        cat("Press return to continue ", "\n")
        readline()
        gmax = max(wk1)
        gmin = min(wk1)
        cx = (gmax - gmin)/10
        gmax = gmax + cx
        gmin = gmin - cx
        for (j in 1:k^2) {
            plot(tdx, wk1[j, ], type = "l", xlab = "lag", ylab = "Acu-IRF", 
                ylim = c(gmin, gmax))
            points(tdx, wk1[j, ], pch = "*")
            title(main = "Orth. innovations")
        }
    }
    else {
        gmax = max(WGT)
        gmin = min(WGT)
        cx = (gmax - gmin)/10
        gmax = gmax + cx
        gmin = gmin - cx
        for (j in 1:k^2) {
            plot(tdx, WGT[j, ], type = "l", xlab = "lag", ylab = "IRF", 
                ylim = c(gmin, gmax))
            points(tdx, WGT[j, ], pch = "*")
            title(main = "Orig. innovations")
        }
        cat("Press return to continue ", "\n")
        readline()
        gmax = max(wk1)
        gmin = min(wk1)
        cx = (gmax - gmin)/10
        gmax = gmax + cx
        gmin = gmin - cx
        for (j in 1:k^2) {
            plot(tdx, wk1[j, ], type = "l", xlab = "lag", ylab = "Acu-IRF", 
                ylim = c(gmin, gmax))
            points(tdx, wk1[j, ], pch = "*")
            title(main = "Orig. innovations")
        }
    }
    VARMAirf <- list(psi = PSI, irf = WGT)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
