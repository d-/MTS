\name{Vech}
\alias{Vech}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Half-Stacking Vector of a Symmetrix Matrix
%%  ~~function to do ... ~~
}
\description{Obtain the half-stacking vector of a symmetric matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Vech(mtx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mtx}{A symmetric matrix
%%     ~~Describe \code{mtx} here~~
}
}
\details{Stacking a matrix into a vector using data on and below the diagonal.
%%  ~~ If necessary, more details than the description above ~~
}
\value{a vector
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
m1=matrix(c(1:9),3,3)
m2=(m1+t(m1))/2
v1=Vech(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
