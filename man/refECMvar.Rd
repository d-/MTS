\name{refECMvar}
\alias{refECMvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining Error-Correction Model for VAR series
%%  ~~function to do ... ~~
}
\description{Refining an estimated ECM VAR model by setting insignificant 
estimates to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refECMvar(m1, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{An object of the ECMvar command or the refECMvar command
%%     ~~Describe \code{m1} here~~
}
  \item{thres}{Threshold for individual t-ratio. The default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{Set simultaneously all estimates with t-ratio less than 
the threshold to zero (in modulus).
%%  ~~ If necessary, more details than the description above ~~
}
\value{Constrained estimation results of a ECM VAR model
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 5)
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
m2=refECMvar(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
