\name{BVAR}
\alias{BVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Vector Autoregression
%%  ~~function to do ... ~~
}
\description{Estimate a VAR(p) model using Bayesian approach, including the 
use of Minnesota prior
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BVAR(zt,p=1,C,V0,n0=5,Phi0=NULL,include.mean=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zt}{A matrix of vector time series, each column represents a series.
%%     ~~Describe \code{zt} here~~
}
  \item{p}{
The AR order. Default is p=1.
%%     ~~Describe \code{p} here~~
}
  \item{C}{
The precision matrix of the coefficient matrix. With constant, 
the dimension of C is (kp+1)-by-(kp+1). The covariance marix of the 
prior for the parameter vec(Beta) is Kronecker(Sigma_a,C-inverse).
%%     ~~Describe \code{C} here~~
}
  \item{V0}{
A k-by-k covariance matrix to be used as prior for the Sigma_a matrix
%%     ~~Describe \code{V0} here~~
}
  \item{n0}{
The degrees of freedom used for prior of the Sigma_a matrix, the covariance matrix of the innovations. Default is n0=5.
%%     ~~Describe \code{n0} here~~
}
  \item{Phi0}{
The prior mean for the parameters. Default is set to NULL, implying that the 
prior means are zero.
%%     ~~Describe \code{Phi0} here~~
}
  \item{include.mean}{
A logical switch controls the constant term in the VAR model. Default is to include the constant term.
%%     ~~Describe \code{include.mean} here~~
}
}
\details{
for a given prior, the program provide the posterior estimates of a 
VAR(p) model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The posterior estimates, their standard errors, t-ratio, and 
the residual covariance matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 2).
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
da=read.table("q-gdp-ukcaus.txt",header=T)
z=log(da[,3:5])
zt=diffM(z)*100
C=0.1*diag(rep(1,7))
V0=diag(rep(1,3))
BVAR(zt,p=2,C,V0)

## The function is currently defined as
"BVAR" <- function(z,p=1,C,V0,n0=5,Phi0=NULL,include.mean=T){
## Perform Bayesian estimation of a VAR(p) model
##
## z: time series (T-by-k)
## p: AR order
## phi0: prior mean for coefficient matrix [k-by-(kp+1)]
## C: precision matrix of coefficient matrix. [(kp+1)-by-(kp+1)]
## (V0,n0): prior input for Sigma_a (inverted Wishart parameters)
##
if(!is.matrix(z))z=as.matrix(z)
if(p < 1) p=1
if(!is.matrix(C))C=as.matrix(C)
if(!is.matrix(V0))V0=as.matrix(V0)
if(n0 < 1)n0=1
##
T=dim(z)[1]
k=dim(z)[2]
idim=k*p+1
if(length(Phi0) <= 0)Phi0=matrix(0,idim,k)
X=NULL
ne=T-p
if(include.mean)X=rep(1,ne)
for (i in 1:p){
X=cbind(X,z[(p+1-i):(T-i),])
}
Z=as.matrix(z[(p+1):T,])
X=as.matrix(X)
### 
XpX=t(X)%*%X
XpY=t(X)%*%Z
## Bayesian Estimate
WpW=XpX+C
WpWinv=solve(WpW)
WpY=XpY+C%*%Phi0
Bbhat=WpWinv%*%WpY
bAhat=Z-X%*%Bbhat
bB=Bbhat-Phi0
S=t(bAhat)%*%bAhat +t(bB)%*%C%*%bB
BSig=(V0+S)/(n0+ne-k-1)
SD=kronecker(BSig,WpWinv)
phi=c(Bbhat)
se=sqrt(diag(SD))
Est=cbind(phi,se,phi/se)
colnames(Est) <- c("Est","s.e.","t-ratio")
cat("Bayesian estimate:","\n")
print(Est)
cat("Covariance matrix: ","\n")
print(BSig)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
