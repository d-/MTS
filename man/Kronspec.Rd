\name{Kronspec}
\alias{Kronspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kronecler Index Specification
%%  ~~function to do ... ~~
}
\description{For a given set of Kronecker indices, the program 
specifies a VARMA model. It gives details of parameter specification.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Kronspec(kdx, output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kdx}{A set of Kronecker indices
%%     ~~Describe \code{kdx} here~~
}
  \item{output}{A logical switch to control output.  Deafult is with 
output.
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{PhiID: specification of the AR matrix polynomial. 
ThetaID: specification of the MA matrix polynomial.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 4)
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
\examples{kdx=c(2,1,1)
m1=Kronspec(kdx)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
