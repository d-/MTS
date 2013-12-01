\name{BEKK11vol}
\alias{BEKK11vol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Volatility of BEKK(1,1) Model
%%  ~~function to do ... ~~
}
\description{Computes the volatility matrices of a fitted BEKK(1,1) model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BEKK11vol(rtn, est, include.mean = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns
%%     ~~Describe \code{rtn} here~~
}
  \item{est}{Parameters of a BEKK model, including the mean, if any.
%%     ~~Describe \code{est} here~~
}
  \item{include.mean}{A logical switch to include the mean vector. Default is to include the mean
%%     ~~Describe \code{include.mean} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Volatility matrices
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
