\name{Kronid}
\alias{Kronid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kronecker Index Identification
%%  ~~function to do ... ~~
}
\description{Find the Kronecker indices of a k-dimensional time series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Kronid(x, plag = 5, crit = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data matrix (T-by-k) of a k-dimensional time series
%%     ~~Describe \code{x} here~~
}
  \item{plag}{The number of lags used to represent the past vector. Default 
is 5.
%%     ~~Describe \code{plag} here~~
}
  \item{crit}{Type-I error used in testing for zero canonical correlations. 
Deafult is 0.05.
%%     ~~Describe \code{crit} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{index: Kronecker indices; tests: Chi-square test statistics
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 4)
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
\examples{da=read.table("flourc.txt")
Kronid(da,plag=3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, plag = 5, crit = 0.05) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nT = dim(x)[1]
    k = dim(x)[2]
    y = as.matrix(x)
    if (plag < 1) {
        plag = floor(log(nT)) + 1
    }
    iend = nT - plag
    past = y[1:iend, ]
    if (plag > 1) {
        for (i in 2:plag) {
            past = cbind(y[i:(iend + i - 1), ], past)
        }
    }
    kdx = rep(0, k)
    found = rep(0, k)
    h = 0
    ist = plag + 1
    futu1 = as.matrix(y[ist:nT, ])
    cat("h = ", h, "\n")
    for (i in 1:k) {
        cat("Component = ", i, "\n")
        s1 = c(i)
        if (i > 1) {
            fnd = found[1:(i - 1)]
            jdx = c(1:(i - 1))[fnd == 0]
            s1 = c(jdx, i)
        }
        futu = as.matrix(futu1[, s1])
        m1 = cancor(past, futu)
        df = dim(futu)[2]
        dp = dim(past)[2]
        deg = dp - df + 1
        seig = m1$cor[df]^2
        cat("square of the smallest can. corr. = ", seig, "\n")
        tst = -(nT - 1 - 0.5 * (dp + df - 1)) * log(1 - seig)
        pv = 1 - pchisq(tst, deg)
        stat = c(tst, deg, pv)
        cat("    test,   df, &  p-value:", "\n")
        print(round(stat, 3))
        if (i > 1) {
            cstar = cbind(cstar, stat)
        }
        else {
            cstar = stat
        }
        if (pv > crit) {
            found[i] = 1
            kdx[i] = h
            cat("A Kronecker index found", "\n")
        }
    }
    cat("=============", "\n")
    while (sum(found) < k) {
        idim = dim(past)[1]
        h = h + 1
        cat("h = ", h, "\n")
        past = past[1:(idim - 1), ]
        futu = futu[1:(idim - 1), ]
        futu1 = y[(ist + h):nT, ]
        for (ii in 1:k) {
            if (found[ii] == 0) {
                cat("Component = ", ii, "\n")
                futu = cbind(futu, futu1[, ii])
                m1 = cancor(past, futu)
                df = dim(futu)[2]
                dp = dim(past)[2]
                deg = dp - df + 1
                seig = m1$cor[df]^2
                cat("Square of the smallest can. corr. = ", seig, 
                  "\n")
                y1 = futu \%*\% (m1$ycoef[, df])
                x1 = past \%*\% (m1$xcoef[, df])
                m2 = acf(y1, lag = h, plot = F)
                acfy = m2$acf[2:(h + 1)]
                m3 = acf(x1, lag = h, plot = F)
                acfx = m3$acf[2:(h + 1)]
                dsq = 1 + 2 * sum(acfx * acfy)
                seig = seig/dsq
                tst = -(nT - 1 - 0.5 * (dp + df - 1)) * log(1 - 
                  seig)
                pv = 1 - pchisq(tst, deg)
                stat = c(tst, deg, pv, dsq)
                cat("    test,     df, p-value & d-hat:", "\n")
                print(round(stat, 3))
                stat = stat[1:3]
                cstar = cbind(cstar, stat)
                if (pv > crit) {
                  found[ii] = 1
                  kdx[ii] = h
                  futu = futu[, 1:(df - 1)]
                  cat("A Kronecker found", "\n")
                }
            }
        }
        cat("============", "\n")
    }
    cat("   ", "\n")
    cat("Kronecker indexes identified:", "\n")
    print(kdx)
    Kronid <- list(index = kdx, tests = cstar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
