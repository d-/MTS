\name{KronEstimMatrix}
\alias{KronEstimMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kronecker Indices Implied Polynomial Matrices
%%  ~~function to do ... ~~
}
\description{This is another program for specifying VAR and VMA 
matrix polynomials for a given set of Kronecker.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KronEstimMatrix(kron)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kron}{A set of Kronekcer indices
%%     ~~Describe \code{kron} here~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{A PhD Student
%%  ~~who you are~~
}
\note{Not used in the MTS package
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{Kronspec
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (kron) 
{
    k = length(kron)
    maxk = kron[1]
    if (k > 1) {
        for (i in 2:k) {
            if (kron[i] > maxk) {
                maxk = kron[i]
            }
        }
        unk = 0
        cols = (k + 1) * (maxk + 1)
        rows = k
        zt = matrix(0, rows, cols)
        at = matrix(0, rows, cols)
        for (i in 1:cols) {
            if (i\%\%(k + 1) == 0) {
                for (q in 1:k) {
                  zt[q, i] = paste("Lag", (i\%/\%(k + 1) - 1))
                  at[q, i] = paste("Lag", (i\%/\%(k + 1) - 1))
                }
            }
        }
        colpermat = (k + 1)
        for (j in 1:k) {
            paste("Z", j, ",t^", j)
            zt[j, j] = "1"
            if (j - 1 >= 1) {
                for (i in 1:(j - 1)) {
                  mink = min(kron[j] + 1, kron[i])
                  if (mink >= 1) {
                    for (u in 1:mink) {
                      q = paste("Alpha", u, ",", i, "^(", j, 
                        ")", "Z", i, ",t", u - 1 - kron[j])
                      zt[j, (-1 * (u - 1 - kron[j])) * colpermat + 
                        i] = "X"
                      unk = unk + 1
                    }
                  }
                }
            }
            for (i in j:k) {
                mink2 = min(kron[j], kron[i])
                if (mink2 >= 1) {
                  for (u in 1:mink2) {
                    q = paste("Alpha", u, ",", i, "^(", j, ")", 
                      "Z", i, ",t", u - 1 - kron[j])
                    zt[j, (-1 * (u - 1 - kron[j])) * colpermat + 
                      i] = "X"
                    unk = unk + 1
                  }
                }
            }
            paste("A", j, ",t^", j)
            at[j, j] = "1"
            if (j - 1 >= 1) {
                for (i in 1:(j - 1)) {
                  if (kron[i] >= (kron[j] + 1)) {
                    q = paste("Alpha", kron[j] + 1, ",", i, "^(", 
                      j, ")", "A", i, ",t")
                    at[j, i] = "X"
                    unk = unk + 1
                  }
                }
            }
            if (kron[j] > 0) {
                for (i in 1:kron[j]) {
                  for (n in 1:k) {
                    q = paste("u", i, ",", n, "^(", j, ")", "A", 
                      "t-", i)
                    at[j, i * colpermat + n] = "X"
                    unk = unk + 1
                  }
                }
            }
        }
        m = 0
        for (i in 1:k) {
            m = m + kron[i]
        }
        sma = 0
        smb = 0
        smt = 0
        for (j in 1:k) {
            if ((j - 1) > 0) {
                for (i in 1:(j - 1)) {
                  sma = sma + min(kron[j] + 1, kron[i])
                }
            }
            if (k >= (j + 1)) {
                for (i in (j + 1):k) {
                  smb = smb + min(kron[j], kron[i])
                }
            }
            smt = smt + sma + smb
            sma = 0
            smb = 0
        }
        Unknowns = m * (1 + k) + smt
        cat("AR Coefficient matrices:", "\n")
        cat(" X = Requires Estimation", "\n")
        cat(" 1 = Constant", "\n")
        cat(" 0 = No Estimation Required", "\n")
        kp1 = k + 1
        for (i in 1:(maxk + 1)) {
            cat("Lag = ", i - 1, "\n")
            ist = (i - 1) * kp1
            ARmtx = zt[, (ist + 1):(ist + k)]
            print(ARmtx)
        }
        cat(" ", "\n")
        cat("MA Coefficient matrices:", "\n")
        for (i in 1:(maxk + 1)) {
            cat("Lag = ", i - 1, "\n")
            ist = (i - 1) * kp1
            MAmtx = at[, (ist + 1):(ist + k)]
            print(MAmtx)
        }
        cat(" ", "\n")
        print("Number of unknown parameters to be estimated:")
        print(Unknowns)
        print("Total X's:")
        print(unk)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
