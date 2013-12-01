\name{SCMmod}
\alias{SCMmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SCM Model Specification
%%  ~~function to do ... ~~
}
\description{For a given set of SCMs and locator of transformation matrix, 
the program specifies a VARMA model via SCM approach for estimation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCMmod(order, Ivor, output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{order}{A k-by-2 matrix of the orders of SCM
%%     ~~Describe \code{order} here~~
}
  \item{Ivor}{A k-dimensioal vector to indicate the location of "1" 
for each component z(i,t) in the transformation matrix.
%%     ~~Describe \code{Ivor} here~~
}
  \item{output}{A logical switch to control output.
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Tmtx: specification of the transformation matrix T. 
ARpar: Specification of the VAR parameters. MApar: specification of the 
VMA parameters.
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
\examples{ord=matrix(c(0,1,1,0,0,1),3,2)
Ivor=c(3,1,2)
m1=SCMmod(ord,Ivor,TRUE)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
