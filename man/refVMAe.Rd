\name{refVMAe}
\alias{refVMAe}
\title{Refining VMA Estimation via the Exact Likelihood Method
}
\description{Refines a fitted VMA model via the VMAe command by setting 
insignificant parameters to zero
}
\usage{
refVMAe(model, thres = 1)
}
\arguments{
  \item{model}{An output object of the command VMAe or the command 
refVMAe itself
}
  \item{thres}{A threshold for individual t-ratio of parameter estimates. 
Deafult is 1.
}
}
\details{The program sets simultaneously all estimates with t-ratios 
less than the threshold (in modulus) to zero. 
}
\value{The same as those of the command VMAe.
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VMAe, refVMA
}
\examples{
th1=matrix(c(0.9,0.2,0.0,0.8),2,2)
sig=matrix(c(4,1,1,1),2,2)
m1=VARMAsim(400,malags=c(1),theta=th1,sigma=sig)
zt=m1$series
m2=VMAe(zt,q=1,include.mean=F)
m3=refVMAe(m2,thres=1.8)
names(m3)
}




