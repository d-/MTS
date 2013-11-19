\name{refKronfit}
\alias{refKronfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining VARMA Estimation via Kronecker Index Approach
%%  ~~function to do ... ~~
}
\description{This program performs model simplification of a fitted 
VARMA model via the Kronecker index approach
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refKronfit(model, thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{The name of a model from the command Kronfit or refKronfit
%%     ~~Describe \code{model} here~~
}
  \item{thres}{A threshold for t-ratio of individual parameter estimate.
The deafult is 1.
%%     ~~Describe \code{thres} here~~
}
}
\details{For a given threshold, the program set a parameter to zero if 
its t-ratio (in modulus) is less than the threshold
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a fitted VARMA model, including model specification and 
residuals
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

\seealso{Kronfit
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{zt=read.table("flourc.txt")
kdx=c(1,1,1)
m1=Kronfit(zt,kdx)
m2=refKronfit(m1,thres=1.2)
names(m2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
