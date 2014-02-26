\name{tfm1}
\alias{tfm1}
\title{Transfer Function Model 1
}
\description{Estimation of a general transfer function model. The model is 
for one input and one output.
}
\usage{
tfm1(y, x, orderN, orderX)
}
\arguments{
  \item{y}{Data vector of dependent variable
}
  \item{x}{Data vector of input (or independent) variable
}
  \item{orderN}{Order (p,d,q) of the disturbance component
}
  \item{orderX}{Order (r,s,b) of the transfer function model, where 
r and s are the degrees of denominator and numerator polynomials and 
b is the delay
}
}
\details{Perform estimation of a general transfer function model
}
\value{
\item{estimate }{Coefficient estimates}
\item{sigma2 }{Residual variance sigma-square}
\item{residuals }{Residual series}
\item{varcoef }{Variance of the estimates}
\item{Nt }{The disturbance series}
}
\references{Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (1994). Time Series Analysis:
Forecasting and Control, 3rd edition, Prentice Hall, Englewood Cliffs, NJ.
}
\author{Ruey S. Tsay
}
\seealso{tfm
}



