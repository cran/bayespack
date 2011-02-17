#######################################################################
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.

exlogPost2 <- function(par, type){
  dm <- c(2, 10, 11, 12, 7, 4)
  dm <- dm[type]
  if(length(par) != dm)
    stop("par of wrong length")
  if(type == 1){ # bivariate skew t
    S <- rbind(c(1,-0.9),c(-0.9,1))
    alpha <- c(0,15)
    out <- sn:::dmst(par, Omega = S, df = 5, alpha = alpha, log = TRUE)
  }
  if(type == 2){ # 10 dim banana shaped
    dim <- 10;b <- 0.03;sigma12 <- 100
    y <- c(par[1], par[2]+b*(par[1]^2-sigma12), par[3:dim])
    cc <- c(1/sqrt(sigma12), rep(1, dim-1))
    out <- -0.5*sum((y*cc)^2)
  }
  if(type == 3){ # climate pressure difference between Darwin and Easter Islands
    x <- 1:168
    y <- c(12.9,11.3,10.6,11.2,10.9,7.5,7.7,11.7,12.9,14.3,10.9,
           13.7,17.1,14,15.3,8.5,5.7,5.5,7.6,8.6,7.3,7.6,12.7,11,
           12.7,12.9,13,10.9,10.4,10.2,8,10.9,13.6,10.5,9.2,12.4,
           12.7,13.3,10.1,7.8,4.8,3,2.5,6.3,9.7,11.6,8.6,12.4,10.5,
           13.3,10.4,8.1,3.7,10.7,5.1,10.4,10.9,11.7,11.4,13.7,
           14.1,14,12.5,6.3,9.6,11.7,5,10.8,12.7,10.8,11.8,12.6,15.7,
           12.6,14.8,7.8,7.1,11.2,8.1,6.4,5.2,12,10.2,12.7,10.2,14.7,
           12.2,7.1,5.7,6.7,3.9,8.5,8.3,10.8,16.7,12.6,12.5,12.5,9.8,
           7.2,4.1,10.6,10.1,10.1,11.9,13.6,16.3,17.6,15.5,16,15.2,11.2,
           14.3,14.5,8.5,12,12.7,11.3,14.5,15.1,10.4,11.5,13.4,7.5,0.6,
           0.3,5.5,5,4.6,8.2,9.9,9.2,12.5,10.9,9.9,8.9,7.6,9.5,8.4,10.7,
           13.6,13.7,13.7,16.5,16.8,17.1,15.4,9.5,6.1,10.1,9.3,5.3,
           11.2,16.6,15.6,12,11.5,8.6,13.8,8.7,8.6,8.6,8.7,12.8,13.2,14,
           13.4,14.8)
    b <- par[1:10]
    sig <- exp(par[11])
    mu <- b[1] + b[5]*cos(6.283185*x/b[2]) + b[6]*sin(6.283185*x/b[2]) +
      b[7]*cos(6.283185*x/b[3]) + b[8]*sin(6.283185*x/b[3]) + 
        b[9]*cos(6.283185*x/b[4]) + b[10]*sin(6.283185*x/b[4])
    out <- sum(dnorm(y, mu, sig, log = TRUE))
    out <- out + dcauchy(b[1], 0, 100, log = TRUE)
    out <- out + sum(dcauchy(b[c(5,6,7,8,9,10)], 0, 10, log = TRUE))
    out <- out + sum(dunif(b[c(2,3,4)], 0, 100, log = TRUE))    
    out <- out + dgamma(sig, 0.1, 0.1, log = TRUE) + par[11]
  }
  if(type == 4){ # pump example
    y <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
    time <- c(94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.05, 
              1.05, 2.1, 10.48)
    thet <- exp(par[1:10])
    alpha <- exp(par[11]);beta <- exp(par[12])
    out <- sum(dpois(y, thet*time, log = TRUE))
    out <- out + sum(dgamma(thet, alpha, beta, log = TRUE))
    out <- out + sum(dexp(alpha, 1, log = TRUE))
    out <- out + sum(dgamma(beta, 0.1, 1, log = TRUE))
    logjakobi <- sum(par[1:12])
    out <- out + logjakobi
  }
  if(type == 5){# parameters in order: theta1-theta4, mu, sigma, tau
    yA <- c(62, 60, 63, 59)
    yB <- c(63, 67, 71, 64, 65, 66)
    yC <- c(68, 66, 71, 67, 68, 68)
    yD <- c(56, 62, 60, 61, 63, 64, 63, 59)
    sigm <- exp(par[6])
    tau <- exp(par[7])
    out <- sum(dnorm(par[1:4], par[5], tau, log = TRUE))
    out <- out + sum(dnorm(yA, par[1], sigm, log = TRUE))
    out <- out + sum(dnorm(yB, par[2], sigm, log = TRUE))
    out <- out + sum(dnorm(yC, par[3], sigm, log = TRUE))
    out <- out + sum(dnorm(yD, par[4], sigm, log = TRUE))
    out <- out + par[7]
  }
  if(type == 6){ # dugong example from WinBUGS manual
    age <- c(1,1.5,1.5,1.5,2.5,4,5,5,7,8,8.5,9,9.5,9.5,10,12,
             12,13,13,14.5,15.5,15.5,16.5,17,22.5,29,31.5)
    lngth <- c(1.80,1.85,1.87,1.77,2.02,2.27,2.15,2.26,2.47,
               2.19,2.26,2.40,2.39,2.41,2.50,2.32,2.32,2.43,
               2.47,2.56,2.65,2.47,2.64,2.56,2.70,2.72,2.57)
    alpha <- par[1]
    beta <- par[2]
    gamma <- (exp(par[3])*1+0.5)/(1+exp(par[3]))
    sigm <- exp(par[4])
    mu <- alpha - beta*gamma^age
    out <- dnorm(alpha, 0, sqrt(1e6), log = TRUE)
    out <- out + dnorm(beta, 0, sqrt(1e6), log = TRUE)
    out <- out + dgamma(1/sigm^2, 0.001, 0.001, log = TRUE)
    out <- out + sum(dnorm(lngth, mu, sigm, log = TRUE))
    logjakobi <- log(exp(-par[3])/(1+exp(-par[3]))^2)+par[4]
    out <- logjakobi + out
  }
  if(is.na(out) | out == -Inf)
    return(-.Machine$double.xmax)
  else
    out
}
