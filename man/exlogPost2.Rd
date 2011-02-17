\name{exlogPost2}
\alias{exlogPost2}
\title{
  Additional example log-posterior densities
}
\description{
  This function provides an interface to 6 example log-densities.
}
\usage{
exlogPost2(par, type)
}
\arguments{
  \item{par}{
    Parameter vector
  }
  \item{type}{
    Integer (from 1 to 6) specifying the log-density to use. See Details.
  }
}
\details{
  References and additional information for the different log-densitiess\cr
  \enumerate{
    \item Bivariate skew t density, 2-d\cr
    Bornkamp, B. (2010), Approximating Probability Densities by Iterated
    Laplace Approximations, Journal of Computational and Graphical
    Statistics
    
    \item Banana shaped density, 10-d\cr
    Version used here is described in
    Bornkamp, B. (2010), Approximating Probability Densities by Iterated
    Laplace Approximations, Journal of Computational and Graphical
    Statistics. This function depends on the dmst function from
    the sn R-package.
    
    \item Climate pressure differences, 11-d\cr
    Model used here is described in
    Bornkamp, B. (2010), Approximating Probability Densities by Iterated
    Laplace Approximations, Journal of Computational and Graphical
    Statistics\cr
    The data set has been obtained from the NIST repository:
    www.itl.nist.gov/div898/strd/nls/data/enso.shtml

    \item Pump example, 12-d\cr
    Model used here is described in
    Clark and Gelfand (2006) "Hierarchical modelling for environmental
    sciences", Oxford University Press, p. 22-23
    \cr

    \item Schools example, 7-d\cr
    Model used here is described in
    Gelman, Carlin, Stern and Rubin (2003), "Bayesian Data Analysis",
    Chapman and Hall, p. 299, Table 11.2/Table 12.1.\cr
    A log-Transformation is used for the positive parameters
    \cr

    \item Dugong example, 4-d\cr
    Model used here is described in
    the WinBUGS manual (Volume 2)
    A log-Transformation is used for the positive parameter and
    a generalized logit transform for the parameter that is in [0.5,1]
    \cr
  }
}
\value{
  The function returns the value of the log density for the
  corresponding problem.
}
\author{
  Bjoern Bornkamp
}
\seealso{
  \code{\link{banint}}
}
\examples{
## 6 log densities and their corresp. modes
mode1 <- c(-0.15,  0.17)
exlogPost2(mode1, type = 1)

mode2 <- c(0, 3, rep(0, 8))
exlogPost2(mode2, type = 2)

mode3 <- c(10.5,11.93,44.12,26.82,2.87,1.28,-1.67,0.44,0.11,1.49,0.76)
exlogPost2(mode3, type = 3)

mode4 <- c(-2.792,-2.209,-2.393,-2.145,-0.499,-0.494,-0.151,-0.151,
           0.416,0.68,-0.167,0.092)
exlogPost2(mode4, type = 4)

mode5 <- c(61.29, 65.87, 67.73, 61.15, 64.01, 0.78, 1.197)
exlogPost2(mode5, type = 5)

mode6 <- c(2.650, 0.963, 1.023, -2.346)
exlogPost2(mode6, type = 6)
}
\keyword{ distribution }
