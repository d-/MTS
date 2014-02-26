\name{refVMA}
\alias{refVMA}
\title{Refining VMA Models
}
\description{Refines a fitted VMA model by setting insignificant parameters to zero
}
\usage{
refVMA(model, thres = 2)
}
\arguments{
  \item{model}{An output object from the command VMA or the refVMA command
}
  \item{thres}{A threshold for individual t-ratio of parameter estimate. 
Default is 2.
}
}
\details{The program simultaneously sets all estimates with t-ratios 
less than the threshold (in modulus) to zero.
}
\value{The same as those of the command VMA.
}
\references{Tsay (2014, Chapter 3)
}
\author{Ruey S. Tsay
}
\seealso{VMA
}
\examples{
theta=matrix(c(0.4,0.5,0.0,0.6),2,2)
sig=matrix(c(3,1,1,1),2,2)
m1=VARMAsim(400,malags=c(1),theta=theta,sigma=sig)
zt=m1$series
m2=VMA(zt,q=1,include.mean=F)
m3=refVMA(m2,thres=1.8)
names(m3)
}




