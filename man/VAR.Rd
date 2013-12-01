\name{VAR}
\alias{VAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vector Autoregressive Model
%%  ~~function to do ... ~~
}
\description{Perform least squares estimation of a VAR model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VAR(x, p = 1, output = T, include.mean = T, fixed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k matrix of k-dimensional time series
%%     ~~Describe \code{x} here~~
}
  \item{p}{Order of VAR model. Default is 1.
%%     ~~Describe \code{p} here~~
}
  \item{output}{A logical switch to control output. Default is with output.
%%     ~~Describe \code{output} here~~
}
  \item{include.mean}{A logical switch. It is true if mean vector is estimated.
%%     ~~Describe \code{include.mean} here~~
}
  \item{fixed}{A logical matrix used in constrained estimation. 
It is used mainly in model simplifcation, e.g., removing insignificant 
estimates.
%%     ~~Describe \code{fixed} here~~
}
}
\details{To remove insignificant estimates, one specifies a threshold 
for individual t-ratio. The fixed matrix is then defined automatically to 
identify those parameters for removal. 
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
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{refVAR command
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
