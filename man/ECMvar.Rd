\name{ECMvar}
\alias{ECMvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Error-Correction VAR Models
%%  ~~function to do ... ~~
}
\description{Performs estimation of an Error-Correction VAR(p) model using 
Quasi Maximum Likelihood Method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ECMvar(x, p, ibeta, include.const = FALSE, fixed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional co-integrated VAR process
%%     ~~Describe \code{x} here~~
}
  \item{p}{VAR order
%%     ~~Describe \code{p} here~~
}
  \item{ibeta}{Initial estimate of the co-integrating matrix.
%%     ~~Describe \code{ibeta} here~~
}
  \item{include.const}{A logical switch to include a constant term in the model. The default is no constant
%%     ~~Describe \code{include.const} here~~
}
  \item{fixed}{A logical matrix to set zero parameter constraints.
%%     ~~Describe \code{fixed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of estimation of an ECM VAR model, including 
alpha and beta matrices of co-integrating information.
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
da=read.table("m-bnd.txt")
bnd=da[,4:5]
beta=matrix(c(1,-0.886),2,1)
m1=ECMvar(bnd,3,ibeta=beta)
names(m1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
