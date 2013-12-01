\name{refSCMfit}
\alias{refSCMfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining Estimation of VARMA Model via SCM Approch
%%  ~~function to do ... ~~
}
\description{Refine estimation of a VARMA model specified via the SCM 
approach by removing insignificant parameters
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refSCMfit(model, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Name of the model from the SCMfit command or the refSCMfit command
%%     ~~Describe \code{model} here~~
}
  \item{thres}{Threshold for the t-ratio of individual coefficient. 
Default is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of estimation, including model specification and residuals
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

\seealso{SCMfit
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{zt=read.table("flourc.txt")
Ord=matrix(c(1,1,2,0,0,0),3,2)
Tdx=c(2,3,1)
m1=SCMfit(zt,Ord,Tdx)
m2=refSCMfit(m1,thres=0.1)
names(m2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
