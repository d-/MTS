\name{refsVARMA}
\alias{refsVARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Refining a Seasonal VARMA Model
%%  ~~function to do ... ~~
}
\description{Refines a fitted seasonal VARMA model by setting insignificant 
estimates to zero
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refsVARMA(model, thres = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An output object of the sVARMA command or the refsVARMA 
command
%%     ~~Describe \code{model} here~~
}
  \item{thres}{Threshold for individual t-ratio. Default is 0.8. 

%%     ~~Describe \code{thres} here~~
}
}
\details{The command removes simultaneously all parameters with t-ratio 
less than the threshold in modulus.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Results of a sesonal VARMA fitting
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Tsay (2014, Chapter 6)
%% ~put references to the literature/web site here ~
}
\author{Ruey S. Tsay
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{sVARMA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m4=refsVARMA(m3,thres=1.2)
names(m4)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
