\name{PIwgt}
\alias{PIwgt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pi Weight Matrices
%%  ~~function to do ... ~~
}
\description{Compute the Pi-weight matrices of a VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PIwgt(Phi = NULL, Theta = NULL, lag = 12, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficients in the form 
[Phi1, Phi2, Phi3, ..., Phip]
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{A k-by-kq matrix of VMA coefficients in the form 
[Theta1, Theta2, ..., Thetaq]
%%     ~~Describe \code{Theta} here~~
}
  \item{lag}{The number of Pi-weight matrices to be computed.
%%     ~~Describe \code{lag} here~~
}
  \item{plot}{A logical switch to plot the Pi-weight matrices
%%     ~~Describe \code{plot} here~~
}
}
\details{The Pi-weight matrices for a VARMA model is 
Pi(B) = inverse(Theta(B)) times Phi(B). 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2013)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{PSIwgt
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Phi1=matrix(0,2,2); Phi2=matrix(c(0.2,-0.6,0.3,1.1),2,2)
Theta1=diag(c(-0.5,-0.4))
Phi=cbind(Phi1,Phi2)
m1=PIwgt(Phi=Phi,Theta=Theta1)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Phi = NULL, Theta = NULL, lag = 12, plot = TRUE) 
{
    m1 = PSIwgt(Phi = Theta, Theta = Phi, lag = lag, plot = FALSE)
    PImtx = m1$psi.weight
    k = dim(PImtx)[1]
    nc = dim(PImtx)[2]
    PImtx[, (k + 1):nc] = -PImtx[, (k + 1):nc]
    lag = floor(nc/k) - 1
    WGT = c(diag(k))
    for (i in 1:lag) {
        cat("Lag: ", i, " pi-matrix", "\n")
        ist = (i - 1) * k
        WGT = cbind(WGT, c(PImtx[, (ist + 1):(ist + k)]))
        print(round(PImtx[, (ist + 1):(ist + k)], 5))
    }
    if (plot) {
        tdx = c(1:(lag + 1)) - 1
        par(mfcol = c(k, k))
        gmax = max(WGT)
        gmin = min(WGT)
        cx = (gmax - gmin)/10
        gmax = gmax + cx
        gmin = gmin - cx
        for (j in 1:k^2) {
            plot(tdx, WGT[j, ], type = "l", xlab = "lag", ylab = "Piwgt", 
                ylim = c(gmin, gmax))
            points(tdx, WGT[j, ], pch = "*")
            title(main = "Pi-weights")
        }
    }
    PIwgt <- list(pi.weight = PImtx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
