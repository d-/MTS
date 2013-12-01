\name{tfm}
\alias{tfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transfer Function Model
%%  ~~function to do ... ~~
}
\description{Estimates a transform function model. This program does not allow 
rational transfer function model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tfm(y, x, b = 0, s = 1, p = 0, q = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Data vector of dependent (output) variable
%%     ~~Describe \code{y} here~~
}
  \item{x}{Data vector of independent variable
%%     ~~Describe \code{x} here~~
}
  \item{b}{deadtime or delay
%%     ~~Describe \code{b} here~~
}
  \item{s}{The order of the transfer function polynomial
%%     ~~Describe \code{s} here~~
}
  \item{p}{AR order of the disturbance
%%     ~~Describe \code{p} here~~
}
  \item{q}{MA order of the disturbance
%%     ~~Describe \code{q} here~~
}
}
\details{The model entertained is y_t = c_0+ v(B)x_t + n_t. 
v(B) = 1- v1*B - ... - vs*B^s, and n_t is an ARMA(p,q) process.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{coef }{Coefficient estimates of the transfer function}
\item{se.coef }{Standard errors of the transfer function coefficients}
\item{coef.arma }{Coefficient estimates of ARMA models}
\item{se.arma }{Standard errors of ARMA coefficients}
\item{nt }{The disturbance series}
\item{residuals }{The residual series}
%% ...
}
\references{Box and Jenkins (1970)
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
