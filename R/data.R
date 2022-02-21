

#' Monthly simple returns of the stocks of International Business Machines
#' (IBM) and Coca Cola (KO) and the S&P Composite index (SP)
#'
#' Monthly simple returns of the stocks of International Business Machines
#' (IBM) and Coca Cola (KO) and the S&P Composite index (SP).  The sample
#' period is from January 1961 to December 2011. The original data were from
#' the Center for Research in Security Prices (CRSP) of the University of
#' Chicago. The files has four columns. They are dates, IBM, SP, and KO.
#'
#'
#' @name ibmspko
#' @docType data
#' @format A 2-d list containing 612x4 observations.  The files has four
#' columns. They are dates, IBM, SP, and KO.
#' @source World Almanac and Book of Facts, 1975, page 406.
#' @keywords datasets
NULL





#' MTS Internal Functions
#'
#' MTS Internal Functions
#'
#' These are not to be called by the user.
#'
#' @aliases Lminv VARchi VARecm VARfore VARirf mFilter VMApred refVARs refVMAs
#' revmq
NULL





#' Multivariate Time Series
#'
#' Multivariate Time Series (MTS) is a general package for analyzing
#' multivariate linear time series and estimating multivariate volatility
#' models. It also handles factor models, constrained factor models, asymptotic
#' principal component analysis commonly used in finance and econometrics, and
#' principal volatility component analysis.  (a) For the multivariate linear
#' time series analysis, the package performs model specification, estimation,
#' model checking, and prediction for many widely used models, including vector
#' AR models, vector MA models, vector ARMA models, seasonal vector ARMA
#' models, VAR models with exogenous variables, multivariate regression models
#' with time series errors, augmented VAR models, and Error-correction VAR
#' models for co-integrated time series. For model specification, the package
#' performs structural specification to overcome the difficulties of
#' identifiability of VARMA models. The methods used for structural
#' specification include Kronecker indices and Scalar Component Models.  (b)
#' For multivariate volatility modeling, the MTS package handles several
#' commonly used models, including multivariate exponentially weighted
#' moving-average volatility, Cholesky decomposition volatility models, dynamic
#' conditional correlation (DCC) models, copula-based volatility models, and
#' low-dimensional BEKK models. The package also considers multiple tests for
#' conditional heteroscedasticity, including rank-based statistics.  (c)
#' Finally, the MTS package also performs forecasting using diffusion index,
#' transfer function analysis, Bayesian estimation of VAR models, and
#' multivariate time series analysis with missing values.Users can also use the
#' package to simulate VARMA models, to compute impulse response functions of a
#' fitted VARMA model, and to calculate theoretical cross-covariance matrices
#' of a given VARMA model.
#'
#' \tabular{ll}{ Package: \tab MTS\cr Type: \tab Package\cr License: \tab
#' Artistic License 2.0\cr }
#'
#' @name MTS-package
#' @aliases MTS-package MTS
#' @docType package
#' @author Ruey S. Tsay and David Wood
NULL





#' Quarterly real gross domestic products of United Kingdom, Canada, and the
#' United States
#'
#' Quarterly real gross domestic products of United Kingdom, Canada, and the
#' United States from the first quarter of 1980 to the second quarter of 2011.
#' The UK and CA data were originally from OECD and the US data from the
#' Federal Reserve Bank of St Louis.
#'
#'
#' @name qgdp
#' @docType data
#' @format A 2-d list containing 126x5 observations. The data set consists of 5
#' columns: name, year, month, UK, CA, and US.
#' @source The data were downloaded from the FRED of the Federal Reserve Bank
#' of St Louis. The UK data were in millions of chained 2006 Pounds, the CA
#' data were in millions of chained 2002 Canadian dollars, and the US data were
#' in millions of chained 2005 dollars.
#' @keywords datasets
NULL





#' Monthly simple returns of ten U.S. stocks
#'
#' Monthly simple returns of ten U.S. stocks. The sample period is from January
#' 2001 to December 2011. Tick symbols of the ten stocks are used as column
#' names for the returns.
#'
#'
#' @name tenstocks
#' @docType data
#' @format A 2-d list containing 132x11 observations.
#' @source The original data were from Center for Research in Security Prices
#' (CRSP) of the University of Chicago.  The first column denotes the dates.
#' @keywords datasets
NULL
