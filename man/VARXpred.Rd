\name{VARXpred}
\alias{VARXpred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{VARX Model Prediction
%%  ~~function to do ... ~~
}
\description{Computes point forecasts of a VARX model. The values of 
exogenous variables must be given.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARXpred(m1, newxt = NULL, hstep = 1, orig = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An output object of VARX or refVARX command
%%     ~~Describe \code{m1} here~~
}
  \item{newxt}{The data matrix of exogenous variables needed in forecasts.
%%     ~~Describe \code{newxt} here~~
}
  \item{hstep}{Forecast horizon
%%     ~~Describe \code{hstep} here~~
}
  \item{orig}{Forecast origin. Default is 0, meaning the last data point.
%%     ~~Describe \code{orig} here~~
}
}
\details{Uses the provided exogenous variables and the model to 
compute forecasts
%%  ~~ If necessary, more details than the description above ~~
}
\value{Point forecasts
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 6)
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
da=read.table("m-gasoil.txt",header=T)
zt=da[,3:4]; xt=da[,5:6]
T=dim(da)[1]
ori=T-5
zt1=zt[1:ori,]; xt1=xt[1:ori,]
m1=VARX(zt1,2,xt1,1)
nxt=xt[(ori+1):T,]
m2=VARXpred(m1,newxt=nxt,hstep=5)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m1, newxt = NULL, hstep = 1, orig = 0) 
{
    zt = m1$data
    xt = m1$xt
    p = m1$aror
    m = m1$m
    Ph0 = m1$Ph0
    Phi = m1$Phi
    Sig = m1$Sigma
    beta = m1$beta
    include.m = m1$include.mean
    nT = dim(zt)[1]
    k = dim(zt)[2]
    dx = dim(xt)[2]
    print(dim(Phi))
    print(dim(beta))
    se = NULL
    if (length(Ph0) < 1) 
        Ph0 = matrix(rep(0, k), k, 1)
    if (hstep < 1) 
        hstep = 1
    if (orig < 1) 
        orig = nT
    if (length(newxt) > 0) {
        if (!is.matrix(newxt)) 
            newxt = as.matrix(newxt)
        h1 = dim(newxt)[1]
        hstep = min(h1, hstep)
        nzt = zt[1:orig, ]
        xt = rbind(xt[1:orig, ], newxt)
        for (i in 1:hstep) {
            tmp = Ph0
            ti = orig + i
            for (i in 1:p) {
                idx = (i - 1) * k
                tmp = tmp + Phi[, (idx + 1):(idx + k)] \%*\% matrix(nzt[ti - 
                  i, ], k, 1)
            }
            if (m > -1) {
                for (j in 0:m) {
                  jdx = j * dx
                  tmp = tmp + beta[, (jdx + 1):(jdx + dx)] \%*\% 
                    matrix(xt[ti - j, ], dx, 1)
                }
            }
            nzt = rbind(nzt, c(tmp))
        }
        mm = VARpsi(Phi, lag = hstep)
        Si = Sig
        se = matrix(sqrt(diag(Si)), 1, k)
        if (hstep > 1) {
            for (i in 2:hstep) {
                idx = (i - 1) * k
                wk = as.matrix(mm$psi[, (idx + 1):(idx + k)])
                Si = Si + wk \%*\% Sig \%*\% t(wk)
                se1 = sqrt(diag(Si))
                se = rbind(se, se1)
            }
        }
        cat("Prediction at origin: ", orig, "\n")
        cat("Point forecasts (starting with step 1): ", "\n")
        print(round(nzt[(orig + 1):(orig + hstep), ], 5))
        cat("Corresponding standard errors: ", "\n")
        print(round(se[1:hstep, ], 5))
    }
    else {
        cat("Need new data for input variables!", "\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
