\name{VARMAcov}
\alias{VARMAcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Autocovariance Matrices of a VARMA Model
%%  ~~function to do ... ~~
}
\description{Uses psi-weights to compute the autocovariance matrices of 
a VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMAcov(Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, trun = 120)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{A k-by-kp matrix consisting of VAR coefficient matrices, 
Phi = [Phi1, Phi2, ..., Phip].
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{A k-by-kq matrix consisting of VMA coefficient magrices,
Theta = [Theta1, Theta2, ..., Thetaq]
%%     ~~Describe \code{Theta} here~~
}
  \item{Sigma}{Covariance matrix of the innovations (k-by-k).
%%     ~~Describe \code{Sigma} here~~
}
  \item{lag}{Number of cross-covariance matrices to be computed. 
Default is 12.
%%     ~~Describe \code{lag} here~~
}
  \item{trun}{The lags of pis-weights used in calculation. Default is 120.
%%     ~~Describe \code{trun} here~~
}
}
\details{Use psi-weight matrices to compute approximate autocovariance 
matrices of a VARMA model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The theoretical auto-covariance matrices of a given VARMA model.
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
Phi=matrix(c(0.2,-0.6,0.3,1.1),2,2)
Sig=matrix(c(4,1,1,1),2,2)
VARMAcov(Phi=Phi,Sigma=Sig)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
