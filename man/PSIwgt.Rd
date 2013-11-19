\name{PSIwgt}
\alias{PSIwgt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Psi Wights Matrices
%%  ~~function to do ... ~~
}
\description{Computes the psi-weight matrices of a VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PSIwgt(Phi = NULL, Theta = NULL, lag = 12, plot = TRUE, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficient matrix. 
Phi=[Phi1, Phi1, ..., Phip]
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{A k-by-kq matrix of VMA coefficient matrix. 
Theta=[Theta1, Theta2, ..., Thetaq]
%%     ~~Describe \code{Theta} here~~
}
  \item{lag}{The number of psi-weight matrices to be computed. 
Deafult is 12.
%%     ~~Describe \code{lag} here~~
}
  \item{plot}{A lofical switch to control plotting of the psi-weights.
%%     ~~Describe \code{plot} here~~
}
  \item{output}{A logical switch to control the output.
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{psi-weights matrices
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
phi=matix(c(0.2,-0.6,0.3,1.1),2,2)
theta=matrix(c(-0.5,0.2,0.0,-0.6),2,2)
m1=PSIwgt(Phi=phi,Theta=theta)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Phi = NULL, Theta = NULL, lag = 12, plot = TRUE, output = FALSE) 
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
    if (k < 1) 
        k = 1
    PSI = diag(k)
    WGT = c(PSI)
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
        WGT = cbind(WGT, c(tmp))
    }
    if (output) {
        for (i in 1:lag) {
            cat("Lag: ", i, " psi-matrix", "\n")
            ist = i * k
            print(round(PSI[, (ist + 1):(ist + k)], 5))
        }
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
            plot(tdx, WGT[j, ], type = "l", xlab = "lag", ylab = "Psiwgt", 
                ylim = c(gmin, gmax))
            points(tdx, WGT[j, ], pch = "*")
            title(main = "Psi-weights")
        }
    }
    PSIwgt <- list(psi.weight = PSI, irf = WGT)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
