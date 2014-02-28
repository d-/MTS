\name{dccFit}
\alias{dccFit}
\title{Dynamic Cross-Correlation Model Fitting
}
\description{Fits a DCC model using either multivariate Gaussian or 
multivariate Student-t innovations. Two types of DCC models are available. The first type is proposed by Engle 
and the other is by Tse and Tsui. Both models appear in the Journal of Business and 
Economic Statistics, 2002.
}
\usage{
dccFit(rt, type = "TseTsui", theta = c(0.9, 0.02), ub = c(0.92, 0.079999), lb = c(0.4, 1e-04), cond.dist = "std", df = 7, m = 0)
}
\arguments{
  \item{rt}{The T-by-k data matrix of k-dimensional standardized asset returns. Typically, they are the 
standardized residuals of the command dccPre.
}
  \item{type}{A logical switch to specify the type of DCC model. 
Type="TseTsui" for Tse and Tsui's DCC model. Type = "Engle" for Engle's DCC model. Default is Tse-Tsui model.
}
  \item{theta}{The initial parameter values for theta1 and theta2
}
  \item{ub}{Upper bound of parameters
}
  \item{lb}{Lower bound of parameters
}
  \item{cond.dist}{Conditional innovation distribution with std for multivariate Student-t innovations.
}
  \item{df}{degrees of freedom of the multivariate Student-t innovations.
}
  \item{m}{For Tse and Tsui method only, m denotes the number of returns 
used in local correlation matrix estimation
}
}
\value{
\item{estimates}{Parameter estimates}
\item{Hessian}{Hessian marix of the estimates}
\item{rho.t}{Time-varying correlation matrices. 
Each row contains elements of a cross-correlation matrix.
}
}
\references{Tsay (2014, Chapter 7). Multivariate Time Series Analysis with R and Financial 
Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{dccPre
}
\examples{
data("mts-examples",package="MTS")
rtn=log(ibmspko[,2:4]+1)
m1=dccPre(rtn,include.mean=TRUE,p=0)
rtn1=m1$sresi
m2=dccFit(rtn1)
# continued from dccPre command demonstration
}




