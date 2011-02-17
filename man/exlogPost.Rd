\name{exlogPost}
\alias{exlogPost}
\title{
  Example log-posterior functions
}
\description{
  This function provides an interface to 21 example log-posteriors from
  the postpack.f file from the Fortran BAYESPACK package.
}
\usage{
exlogPost(par, type)
}
\arguments{
  \item{par}{
    Parameter vector
  }
  \item{type}{
    Integer (from 1 to 21) specifying the log-posterior to use. See Details.
  }
}
\details{
  References and additional information for the different log-posteriors\cr
  \enumerate{
    \item Birthweight Data, 7-d\cr
    Venables, W. N., and Ripley, B. D. (1994),
    "Modern Applied Statistics with S-Plus", Springer-Verlag, New York, 
    pp. 193-195.
    
    \item Bio Oxygen Data, 2-d\cr
    Bates, D. M., and Watts, D. G. (1988),
    Nonlinear Regression and Its Applications, Wiley, New York,
    pp. 40-43, 270.
    theta_1  = 60exp(par(1))/( 1 + exp(par(1) )\cr
    theta_2  = 6exp(par(2))/( 1 + exp(par(2) )
    
    \item Contingency Table Data, 9-d\cr
    Example 2 from M. Evans and T. Swartz 
    "Methods for Approxmating Integrals with Applications to Statistics
    with Special Emphasis on Bayesian Integration Problems"
    Statistical Science, 10 (1995), 254-272.\cr
    theta = exp(par(1))/( 2 + 2*exp(par(1)) ),\cr
    alpha_{1,1} = exp(par(2))/( 1 + exp(par(2)) + exp(par(3)) ),\cr
    alpha_{2,1} = exp(par(3))/( 1 + exp(par(2)) + exp(par(3)) ),\cr
    alpha_{1,2} = exp(par(4))/( 1 + exp(par(4)) + exp(par(5)) ),\cr
    alpha_{2,2} = exp(par(5))/( 1 + exp(par(4)) + exp(par(5)) ),\cr
    beta_{1,1} = exp(par(6))/( 1 + exp(par(6)) + exp(par(8)) ),\cr
    beta_{2,1} = exp(par(7))/( 1 + exp(par(6)) + exp(par(8)) ),\cr
    beta_{1,2} = exp(par(8))/( 1 + exp(par(8)) + exp(par(9)) ) and\cr
    beta_{2,2} = exp(par(9))/( 1 + exp(par(8)) + exp(par(9)) ).
    
    \item Extreme Value Distribution, 5-d\cr
    Exact means should all be approximately -.57721
    \item Fuller Data, 3-d\cr
    Fuller, W.A. (1976),
    "Introduction to Statistical Time Series"
    John Wiley and Sons, New York, p. 228.
    \item Stanford Heart Transplant Data, 3-d\cr
    Naylor, J. C. and Smith, A. F. M. (1988),  "Econometric 
    Illustrations of Novel Numerical Integration Strategies for 
    Bayesian Inference", J. Econometrics, 38, pp. 103-125.\cr
    Original variables are tranformed using lambda = exp(par(1)),
    tau = exp(par(2)), and p = exp(par(3)). 
    \item Econometric Data, 3-d\cr
    Kloek, T., and van Dijk, H. K. (1978), 
    ``Bayesian Estimates of Equation System Parameters: 
    An  Application of Integration by Monte Carlo,''
    Econometrica, 46, 1-19.\cr
    Original variables are tranformed using \cr
    t_1 = exp(par(1))/( 1 + exp(par(1)) + exp(par(2)) ),\cr
    t_2 = exp(par(2))/( 1 + exp(par(1)) + exp(par(2)) ) and\cr
    t_3 = exp(par(3))/( 1 + exp(par(3)) ).
    
    \item Econometric Data, 5-d\cr
    Kloek, T., and van Dijk, H. K. (1978), 
    ``Bayesian Estimates of Equation System Parameters: 
    An  Application of Integration by Monte Carlo,''
    Econometrica, 46, 1-19.\cr

    Original variables are tranformed using \cr
    t_3 = exp(par(3))/( 1 + exp(par(3)) + exp(par(4)) ),\cr
    t_4 = exp(par(4))/( 1 + exp(par(3)) + exp(par(4)) ) and\cr
    t_5 = exp(par(5))/( 1 + exp(par(5)) ).
    
    \item Econometric Data, 8-d\cr
    Kloek, T., and van Dijk, H. K. (1978), 
    ``Bayesian Estimates of Equation System Parameters: 
    An  Application of Integration by Monte Carlo,''
    Econometrica, 46, 1-19.\cr

    Original variables are tranformed using \cr
    t_3 = exp(par(3))/( 1 + exp(par(3)) + exp(par(4)) ),\cr
    t_4 = exp(par(4))/( 1 + exp(par(3)) + exp(par(4)) ),\cr
    t_5 = exp(par(5))/( 1 + exp(par(5)) ),\cr
    t_6 = exp(par(6)) and\cr
    t_8 = exp(par(8)).

    \item Lubricant Data, 10-d\cr
    Bates, D. M., and Watts, D. G. (1988),
    Nonlinear Regression and Its Applications, Wiley, New York,
    pp. 87-89, 275. \cr
    Original variable sigma = exp(theta(10)) and original variables
    theta_1 and theta_2 are scaled by 1/100 and 1/10, respectively.
    
    \item Multivariate Logistic Distribution, 5-d\cr
    \item Linear Model Data, 10-d\cr
    Evans, M. and Swartz, T. (1995), 
    "Methods for Apprixmating Integrals with Applications to Statistics
    with Special Emphasis on Bayesian Integration Problems"
    Statistical Science, 10, 254-272.\cr
    Original variable sigma is transformed using sigma = exp(par(10)).
    
    \item Motorette Data, 3-d\cr
    Tanner, M. A. (1993), Tools for Statistical Inference, 2nd Ed.,
    Springer-Verlag, New York, p. 41.\cr
    Original variables are transformed using
    sigma  = exp( par(3) )

    \item Nonlinear Regression Model, 3-d \cr
    P. M. Reilly, "The Numerical Computation of Posterior
    Distributions iin Bayesian Statistical Inference",
    Appl Statist. 25 (1976), pp. 201-209.\cr
    Original variables are transformed using \cr
    alpha = exp(par(1))\cr
    beta  = exp(par(2))\cr
    sigma = exp(par(3))\cr

    \item Multivariate Normal Distribution, 6-d\cr
    \item Porgi Data, 14-d\cr
    Shaw, J. E. H. (1988), Aspects of Numerical Integration and 
    Summarisation, in Bayesian Statistics 3, J. Bernardo, H. H. Degroot, 
    D. V. Lindley and A. F. M. Smith (Eds.), Oxford University Press,   
    pp. 411-428.
    
    \item Photocarcinogen Data, 5-d\cr
    Dellaportas, P. amd Smith, A. F. M. (1993), Bayesian Inference for 
    Generalized Linear and Proportional Hazards Models via Gibbs
    Sampling, Appl. Statist. 42, 443-459.\cr
    Original variable p is transformed p  = exp(par(1))

    \item Proportional Hazards Data, 7-d\cr
    Dellaportas, P., and Wright, D. (1992), 
    A Numerical Integration Strategy in Bayesian Analysis,
    in Bayesian Statistics 4, J. M. Bernardo, et al. eds, 
    Oxford University Press, 601-606.\cr
    Original variables p is tranformed using p = exp(par(1)). 
    
    \item Pearson Type IV, 1-d\cr
     Genz, A., and Kass, R. (1997),
     Subregion Adaptive Integration of Functions Having a Dominant Peak,
     J. Comp. Graph. Stat., 6, 92-111\cr
     Mode at X = 32, mean at X = 160/3.

     \item Radiotherapy Data, 2-d\cr
     Tanner, M. A. (1993), Tools for Statistical Inference, 2nd Ed.,
     Springer-Verlag, New York, p. 14.

     \item Tornado Data, 11-d\cr
     Monahan, J. F., Schrab, K. J., and Anderson, C. E. (1993), 
     ``Statistical Methods for Forecasting Tornado Intensity,'' 
      in Statistical Sciences and Data Analysis, 
      K. Matushita, et al., eds., VSP, Utrecht, The Netherlands, 13-24.
  }
}
\value{
  The function returns the value of the log posterior density for the
  corresponding problem.
}
\author{
  Alan Genz (R interface Bjoern Bornkamp)
}
\seealso{
  \code{\link{banint}}
}
\examples{
## 21 log posterior examples and their corresp. modes
mode1 <- c(-1.80561, 0.01919,-0.00560, 1.89159, 2.65261,
          0.49013,-1.05303)
exlogPost(mode1, type = 1)

mode2 <- c(-0.786, -2.25)
exlogPost(mode2, type = 2)

mode3 <- c(1.40560,-2.14270,-0.47991, 1.63775,-0.33202,-1.73717,
         -0.03699, 2.63396, 1.66428)
exlogPost(mode3, type = 3)

mode4 <- rep(0,5)
exlogPost(mode4, type = 4)

mode5 <- c(141.6062, -83.9727, 1.34778)
exlogPost(mode5, type = 5)

mode6 <- c(3.38498, -0.09179, -0.72319)
exlogPost(mode6, type = 6)

mode7 <- c(0.01221, -1.61596, -0.53693)
exlogPost(mode7, type = 7)

mode8 <- c(4.0859, -0.7206,  0.1162, -1.4866, -0.5329)
exlogPost(mode8, type = 8)

mode9 <- c(3.6524, -0.8612,  0.2569, -1.3141, -0.5332,
         -1.5009, 0.0709, -3.4271)
exlogPost(mode9, type = 9)

mode10 <- c(10.55, 20.66, 1.46, -0.259, 0.0225, 0.402,
          .0352, 57.4, -0.476, -3.198)
exlogPost(mode10, type = 10)

mode11 <- rep(0,5)
exlogPost(mode11, type = 11)

mode12 <- c(2.01819, 0.11640, 0.88300, 0.02861,-0.40753,-0.17198,-0.27059,
         -0.62908, 0.36961,-0.23244)
exlogPost(mode12, type = 12)

mode13 <- c(-4.9306, 3.7471, -1.8631)
exlogPost(mode13, type = 13)

mode14 <- c(1.3862, 0.8676, -5.969)
exlogPost(mode14, type = 14)

mode15 <- rep(0,6)
exlogPost(mode15, type = 15)

mode16 <- c(0.0172,-0.2626,-0.3833, 1.4509, 0.1196, 0.0934, 1.4998, 
          0.3177, 0.5976, 1.2672, 0.2683, 0.6650, 1.2261, 0.3795)
exlogPost(mode16, type = 16)

mode17 <- c(1.1903, -10.8469, -1.1791, -0.3517, 0.4020)
exlogPost(mode17, type = 17)

mode18 <- c(0.13732,-3.90576, 1.87334,-0.13476,-0.02137,
         -0.03869, 0.13793)
exlogPost(mode18, type = 18)

mode19 <- 32
exlogPost(mode19, type = 19)

mode20 <- c(3.81944, -0.08648)
exlogPost(mode20, type = 20)

mode21 <- c(0.6584,  0.2414,  0.1538, -0.8338,  0.1092, 17.8495, 18.3402,
         20.0268, 21.2174, 23.2027, 28.2810)
exlogPost(mode21, type = 21)
}
\keyword{ distribution }
