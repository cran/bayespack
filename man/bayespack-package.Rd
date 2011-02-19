\name{bayespack-package}
\alias{bayespack-package}
\alias{bayespack}
\docType{package}
\title{
Numerical Integration for Bayesian Inference
}
\description{
             The bayespack package provides an R interface to the Fortran
	     BAYESPACK integration routines written by Alan Genz. Given
	     an unnormalized posterior distribution, the BAYESPACK
	     routines approximate the normalization constant and the
	     mean and covariance of the posterior. Posterior means of
	     additional functions of the parameters can be calculated as
	     well.
}
\details{
\tabular{ll}{
Package: \tab bayespack\cr
Type: \tab Package\cr
Version: \tab 1.0-2\cr
Date: \tab 2011-02-19\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Alan Genz (BAYESPACK Fortran source code) and
Bjoern Bornkamp (R interface and minor adaptions in Fortran source code)

Maintainer: Bjoern Bornkamp <bbnkmp@gmail.com>
}
\references{
  Genz, A. and Bretz F. (2009). Computation of Multivariate Normal and t
  Probabilities, Springer, Chapter 6.2

  Genz, A., and Kass, R. (1997). Subregion Adaptive Integration of
  Functions Having a Dominant Peak, Journal Computational Graphical
  Statistics, 6, 92-111
}
\keyword{ package }
\examples{
## Stanford Heart Transplant Problem (using default values of banint)
mode <- rep(0,3)
bypk <- banint(exlogPost, start = mode, type = 6)
## basic information
bypk
## more details
summary(bypk)

## now calculate the means of a transformation as well
mnFoo <- function(x, ...) exp(x)
bypk2 <- banint(exlogPost, start = mode, type = 6, mnFns = mnFoo)

## Dugongs Example (from WinBUGS manual)
mode <- c(2.650, 0.963, 1.023, -2.346)
## parameter transformations needed
tr1 <- function(x) (exp(x) * 1 + 0.5)/(1 + exp(x))
mnF <- function(x, ...) c(x[1], x[2], tr1(x[3]), exp(x[4]))
bypk <- banint(exlogPost2, start = mode, type = 6, mnFns = mnF)
## compare to WinBUGS manual
summary(bypk)
}
