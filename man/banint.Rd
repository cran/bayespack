\name{banint}
\alias{banint}
\title{
  Interface to BAYESPACK integration routines
}
\description{
  This function provides an interface to the BAYESPACK integration
  routines written in Fortran by Alan Genz. The code performs a Laplace
  approximation of the posterior distribution, and uses this to
  transform the original integrand. The code then implements classical
  Monte Carlo importance sampling, but also lattice rule methods,
  Gauss-Hermite quadrature, radial-spherical and subregion adaptive
  integration routines (see below or ?viewDoc for details).
}
\usage{
banint(logPost, mnFns = function(x, ...) NULL, mode = NULL,
       start = NULL, method = "Gauss-Hermite", trans = "split-t",
       control = list(), optimctrl = list(), optimMethod = "BFGS", ...)
}
\arguments{
  \item{logPost}{
    Function that returns the (possibly non-normalized) log-posterior density
  }
  \item{mnFns}{
    Function that returns additional means that should be calculated in
    a vector. When the ... argument is used in \code{logPost} also the mnFns function
    should include a ... argument even if not used.
  }
  \item{mode}{
    Mode to assume for Laplace approximation
  }
  \item{start}{
    Starting value for optimization (only needed if the mode is not
    specified) to calculate the Laplace approximation
  }
  \item{method}{
    Integration Method, should be one of "Monte Carlo", "Lattice Rule",
    "Subregion-Adaptive","Mixed Spherical-Radial", "Gauss-Hermite",
    "Modified Gauss-Hermite" and "Stochastic Spherical-Radial".
  }
  \item{trans}{
    Type of transformation to use (default: "split-t"), other possible
    options are "none", multivariate normal ("mvn") and multivariate t
    ("mvt"). Note that for "Mixed Spherical-Radial", "Gauss-Hermite",
    "Modified Gauss-Hermite" and "Stochastic Spherical-Radial" use of
    "mvt" is not possible.
  }
  \item{control}{
    List of control elements for the integration code. Should be a list
    with entries \code{relreq} specifying the relative error tolerance for
    all integration results (default: 0.01), \code{maxvls} the
    maximum number of function evaluations for the integration (default:
    10000) and \code{df} an argument specifying the degrees of freedom
    if a multivariate t transformation is selected (default: df=4). The
    \code{print} entry is a logical, specifying whether information
    should be put out directly from the Fortran code (in a very
    rudimentary form).
  }
  \item{optimctrl}{
    List containing control arguments for the optim function. See help
    page of the optim function for details.
  }
  \item{optimMethod}{
    Optimization method to be used, passed down to the method argument
    of the optim function. See help page of the optim function for details.
  }
  \item{\dots}{
    Additional arguments needed for evaluation of \code{logPost} (or \code{mnFns})
}
}
\details{
  The banint function is an interface to the BAYESPACK Fortran routines
  by Alan Genz. There are 7 integration methods implemented:
  "Monte Carlo", "Lattice Rule","Subregion-Adaptive","Mixed Spherical-Radial",
  "Gauss-Hermite", "Modified Gauss-Hermite" and "Stochastic
  Radial-Spherical". The type of transformations
  that can be used are "none", "mvn", "mvt" and "split-t". See
  the references below for more information. 
}
\value{
  The function returns an object of class bayespack, for which print and
  summary methods exist.

  The object is a list, its most interesting entries are
  \item{method, trans}{The used integration method and transformation}
  \item{m, mn}{Length of parameter vector and additional means}
  \item{mode}{The mode, either calculated or handed over to the code}
  \item{max}{The function value of the log-posterior at the mode}
  \item{nrmcon}{The log-normalization constant (with error estimates)}
  \item{neval}{The number of function evaluations (used for
    integration)}
  \item{inform}{Information whether integration was completed within
    error tolerance}
  \item{means}{The calculated means}
  \item{errors}{The error estimates associated to the means}
  \item{cov}{The calculated covariance matrix}
  \item{modcov}{The modal covariance matrix (calculated from the Hessian
    matrix)}
}
\references{
  Genz, A. and Bretz F. (2009). Computation of Multivariate Normal and t
  Probabilities, Springer, Chapter 6.2

  Genz, A., and Kass, R. (1997). Subregion Adaptive Integration of
  Functions Having a Dominant Peak, Journal Computational Graphical
  Statistics, 6, 92-111
}
\author{
  Alan Genz (BAYESPACK Fortran code), Bjoern Bornkamp (R interface)
}

\seealso{
    \code{\link{exlogPost}}
}
\examples{
## Stanford Heart Transplant Problem (using default values)
strt <- rep(0,3)
bypk <- banint(exlogPost, start = strt, type = 6)
## basic information
bypk
## more details
summary(bypk)

## now calculate the means of a transformation as well
mnFoo <- function(x, ...) exp(x)
bypk2 <- banint(exlogPost, start = strt, type = 6, mnFns = mnFoo)
summary(bypk2)

## now use a different integration routine and transformation
bypk3 <- banint(exlogPost, start = strt, type = 6, mnFns = mnFoo,
                method = "Lattice Rule", control = list(numtrn = "mvn"))

## Dugongs Example (from WinBUGS manual)
mode <- c(2.650, 0.963, 1.023, -2.346)
## parameter transformations needed
tr1 <- function(x) (exp(x) * 1 + 0.5)/(1 + exp(x))
mnF <- function(x, ...) c(x[1], x[2], tr1(x[3]), exp(x[4]))
bypk <- banint(exlogPost2, start = mode, type = 6, mnFns = mnF)
## compare to WinBUGS manual
summary(bypk)
}
\keyword{ distribution }
