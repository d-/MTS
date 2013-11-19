\name{SWfore}
\alias{SWfore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stock-Watson Diffusion Index FOrecasts
%%  ~~function to do ... ~~
}
\description{Uses the diffusion index approach of Stock and Watson to 
compute out-of-sample forecasts
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SWfore(y, x, orig, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The scalar variable of interest 
%%     ~~Describe \code{y} here~~
}
  \item{x}{The data matrix (T-by-k) of the observed explanatory variables
%%     ~~Describe \code{x} here~~
}
  \item{orig}{Forecast origin
%%     ~~Describe \code{orig} here~~
}
  \item{m}{The number of diffusion index used
%%     ~~Describe \code{m} here~~
}
}
\details{Performs PCA on X at the forecast origin. Then, fit a linear 
regression model to obtain the coefficients of prediction equation. 
Use the prediction equation to produce forecasts and compute forecast 
errors, if any.
%%  ~~ If necessary, more details than the description above ~~
}
\value{coef: beta coefficient. yhat: predictions. MSE: mean squared error if 
available. loadings: loading matrix. DFindex: diffusion indices
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
da=read.table("m-unempstatesAdj.txt",header=T)
drate=diffM(da)
yt=drate[5:415,1]
x=cbind(drate[4:414,],drate[3:413,],drate[2:412,],drate[1:411,])
m1=SWfore(y,x,350,10)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, x, orig, m) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nT = dim(x)[1]
    k = dim(x)[2]
    if (orig > nT) 
        orig = nT
    if (m > k) 
        m = k
    if (m < 1) 
        m = 1
    x1 = x[1:orig, ]
    me = apply(x1, 2, mean)
    se = sqrt(apply(x1, 2, var))
    x1 = x
    for (i in 1:k) {
        x1[, i] = (x1[, i] - me[i])/se[i]
    }
    V1 = cov(x1[1:orig, ])
    m1 = eigen(V1)
    sdev = m1$values
    M = m1$vectors
    M1 = M[, 1:m]
    Dindex = x1 \%*\% M1
    y1 = y[1:orig]
    DF = Dindex[1:orig, ]
    mm = lm(y1 ~ DF)
    coef = matrix(mm$coefficients, (m + 1), 1)
    yhat = NULL
    MSE = NULL
    if (orig < nT) {
        newx = cbind(rep(1, (nT - orig)), Dindex[(orig + 1):nT, 
            ])
        yhat = newx \%*\% coef
        err = y[(orig + 1):nT] - yhat
        MSE = mean(err^2)
        cat("MSE of out-of-sample forecasts: ", MSE, "\n")
    }
    SWfore <- list(coef = coef, yhat = yhat, MSE = MSE, loadings = M1, 
        DFindex = Dindex)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
