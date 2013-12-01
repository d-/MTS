\name{FEVdec}
\alias{FEVdec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forecast Error Variance Decomposition
%%  ~~function to do ... ~~
}
\description{Computes the forecast error variance decomposition of a 
VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FEVdec(Phi, Theta, Sig, lag = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{VAR coefficient matrices in the form Phi=[Phi1, Phi2, ..., Phip], 
a k-by-kp matrix.
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{VMA coefficient matrices in form form Theta=[Theta1, Theta2, 
..., Thetaq], a k-by-kq matrix.
%%     ~~Describe \code{Theta} here~~
}
  \item{Sig}{The residual covariance matrix Sigma, a k-by-k positive definite 
matrix.
%%     ~~Describe \code{Sig} here~~
}
  \item{lag}{The number of lags of forecast errors variance to be computed. 
Default is 4.
%%     ~~Describe \code{lag} here~~
}
}
\details{Use the psi-weight matrices to compute the forecast error covariance 
and use Cholesky decomposition to perform the decomposition
%%  ~~ If necessary, more details than the description above ~~
}
\value{irf: impulse response matrices, orthirf: orthogonal impulse response matrices, Omega: forecast error variance matrices, and OmegaR: forecast error 
variance decomposition
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
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
theta1=matrix(c(-0.5,0,0,-0.6),2,2)
Sig=matrix(c(3,1,1,1),2,2)
m1=FEVdec(p1,theta1,Sig)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
