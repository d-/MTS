\name{refVARX}
\alias{refVARX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining a VARX Model
%%  ~~function to do ... ~~
}
\description{Refine a fitted VARX model by setting insignificant parameters 
to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refVARX(m1, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An output object of the VARX commond or the refVARX command
%%     ~~Describe \code{m1} here~~
}
  \item{thres}{A threshold for the individual t-ratio. Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{The program sets simultaneously all estimates with t-ratio less than 
threshold (in modulus) to zero and re-estimate the VARX model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a VARX estimation 
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

\seealso{VARX
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-gasoil.txt",header=T)
zt=da[,3:4]; xt=da[,5:6]
m1=VARX(zt,2,xt,1)
m2=refVARX(m1,thres=1.5)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m1, thres = 1) 
{
    zt = m1$data
    xt = m1$xt
    p = m1$aror
    m = m1$m
    include.m = m1$include.mean
    beta = m1$coef
    se.beta = m1$se.coef
    fix = matrix(0, nrow(beta), ncol(beta))
    for (i in 1:ncol(beta)) {
        tt = beta[, i]/se.beta[, i]
        idx = c(1:nrow(beta))[abs(tt) > thres]
        if (length(idx) > 0) {
            fix[idx, i] = 1
        }
    }
    mm = VARX(zt, p, xt, m, include.mean = include.m, fix)
    Ph0 = mm$Ph0
    Phi = mm$Phi
    Beta = mm$beta
    resi = mm$residuals
    sig = mm$Sigma
    coef = mm$beta
    se.coef = mm$se.coef
    refVARX <- list(data = zt, aror = p, xt = xt, m = m, Ph0 = Ph0, 
        Phi = Phi, beta = Beta, residuals = resi, Sigma = sig, 
        coef = beta, se.coef = se.beta, include.mean = include.m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
