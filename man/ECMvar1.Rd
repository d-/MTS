\name{ECMvar1}
\alias{ECMvar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Error-Correction VAR Model 1
%%  ~~function to do ... ~~
}
\description{Perform least-squares estimation of an ECM VAR(p) model with 
known co-integrating processes
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ECMvar1(x, p, wt, include.const = FALSE, fixed = NULL, output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional co-integrated VAR process
%%     ~~Describe \code{x} here~~
}
  \item{p}{VAR order
%%     ~~Describe \code{p} here~~
}
  \item{wt}{A T-by-m data matrix of m-dimensional co-integrated process
%%     ~~Describe \code{wt} here~~
}
  \item{include.const}{A logical switch to include a constant term. 
Default is no constant.
%%     ~~Describe \code{include.const} here~~
}
  \item{fixed}{A logical matrix to set zero parameter constraints
%%     ~~Describe \code{fixed} here~~
}
  \item{output}{A logical switch to control output
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of least-squares estimation of an ECM mdoel for VAR process.
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

\seealso{ECMvar
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
da=read.table("m-bnd.txt")
bdn=da[,4:5]
wt=bnd[,1]-0.886*bnd[,2]
m1=ECMvar1(bnd,3,wt)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
