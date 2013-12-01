\name{VARMAirf}
\alias{VARMAirf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Impulse Response Funcation of a VARMA Model
%%  ~~function to do ... ~~
}
\description{Compute the impulse response function of a given VARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARMAirf(Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, orth = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficients. 
Phi=[Phi1, Phi2, ..., Phip].
%%     ~~Describe \code{Phi} here~~
}
  \item{Theta}{A k-by-kq matrix of VMA coefficients. 
Theta=[Theta1, Theta2, ..., Thetaq]
%%     ~~Describe \code{Theta} here~~
}
  \item{Sigma}{Covariance matrix (k-by-k) of the innovations.
%%     ~~Describe \code{Sigma} here~~
}
  \item{lag}{Number of lags of impulse response functions to be computed
%%     ~~Describe \code{lag} here~~
}
  \item{orth}{A logical switch to use orthogonal innovations.
Deafult is to perform orthogonalization of the innovations.
%%     ~~Describe \code{orth} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Plots and values of impulse response functions and cumulative 
impulse response function.
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

\seealso{VARMApsi command
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
th1=matrix(c(-0.5,0.2,0.0,-0.6),2,2)
Sig=matrix(c(4,1,1,1),2,2)
m1=VARMAirf(Phi=p1,Theta=th1,Sigma=Sig)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
