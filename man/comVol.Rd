\name{comVol}
\alias{comVol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Common Volatility
%%  ~~function to do ... ~~
}
\description{Compute the principal volatility components after a VAR(p) 
fitting
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
comVol(rtn, m = 10, p = 1, stand = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns
%%     ~~Describe \code{rtn} here~~
}
  \item{m}{The number of lags used to compute generalized cross-Kurtosis 
matrix
%%     ~~Describe \code{m} here~~
}
  \item{p}{VAR order for the mean equation
%%     ~~Describe \code{p} here~~
}
  \item{stand}{A logical switch to standardize the returns
%%     ~~Describe \code{stand} here~~
}
}
\details{Perform a VAR(p) fits, if any. Then, use the residual series to 
perform principal volatility component analysis. The ARCH test statistics 
are also computed for the sample principal components
%%  ~~ If necessary, more details than the description above ~~
}
\value{residuals: The residuals of VAR fit, which are are data used in 
principal volatility component analysis. Values: eigenvalues. Vectors: 
eigenvectors. M: the transformation matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 7)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay and Y.B. Hu
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
da=read.table("w-7fx.txt",header=T)
fx=log(da[,4:10])
rtn=diffM(fx)
m1=comVol(rtn,p=5)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rtn, m = 10, p = 1, stand = FALSE) 
{
    if (!is.matrix(rtn)) 
        rtn = as.matrix(rtn)
    x = VARfitC(rtn, p)$residuals
    nT = dim(x)[1]
    k = dim(x)[2]
    if (m < 1) 
        m = 1
    V1 = cov(x)
    m1 = eigen(V1)
    D1 = diag(1/sqrt(m1$values))
    P1 = m1$vectors
    Shalf = P1 \%*\% D1 \%*\% t(P1)
    x1 = x \%*\% Shalf
    A = matrix(0, k, k)
    for (h in 1:m) {
        ist = h + 1
        for (i in 1:k) {
            for (j in i:k) {
                Cmtx = matrix(0, k, k)
                y2 = x1[(ist - h):(nT - h), i] * x1[(ist - h):(nT - 
                  h), j]
                for (ii in 1:k) {
                  for (jj in ii:k) {
                    y1 = x1[ist:nT, ii] * x1[ist:nT, jj]
                    Cmtx[ii, jj] = cov(y1, y2) * (nT - h)/nT
                    Cmtx[jj, ii] = Cmtx[ii, jj]
                  }
                }
                Cmtx = Cmtx * ((nT - h)/nT)
                A = A + Cmtx \%*\% Cmtx
            }
        }
    }
    if (stand) {
        dd = diag(A)
        D = diag(1/sqrt(dd))
        A = D \%*\% A \%*\% D
    }
    else {
        A = A/(k * (k + 1)/2)
    }
    m2 = eigen(A)
    Valu = m2$values
    Prop = Valu/sum(Valu)
    cat("eigen-values: ", Valu, "\n")
    cat("proportion:   ", Prop, "\n")
    Vec = m2$vectors
    Mmtx = Shalf \%*\% Vec
    for (j in 1:k) {
        Mmtx[, j] = Mmtx[, j]/sqrt(sum(Mmtx[, j]^2))
    }
    Tst = NULL
    Tx = x \%*\% Mmtx
    for (i in 1:k) {
        TT = NULL
        mtst10 = archTstC(Tx[, i], 10)
        TT = c(TT, mtst10)
        mtst20 = archTstC(Tx[, i], 20)
        TT = c(TT, mtst20)
        mtst30 = archTstC(Tx[, i], 30)
        TT = c(TT, mtst30)
        Tst = rbind(Tst, c(i, TT))
    }
    cat("Checking: ", "\n")
    cat("Results of individual F-test for ARCH effect", "\n")
    cat("Numbers of lags used: 10, 20, 30", "\n")
    cat("Component,(F-ratio P-val) (F-ratio P-val) (F-ratio P-Val)", 
        "\n")
    print(Tst, digits = 3)
    comVol <- list(residuals = x, values = m2$values, vectors = m2$vectors, 
        M = Mmtx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
