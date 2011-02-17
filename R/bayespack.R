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

banint <- function(logPost, mnFns = function(x, ...) NULL,
                   mode = NULL, start = NULL, method = "Gauss-Hermite",
                   trans = "split-t", control = list(), optimctrl = list(),
                   optimMethod = "BFGS", ...){
  if(is.null(start) & is.null(mode))
    stop("either mode or start argument needed")
  if(is.null(mode)){
    fn <- function(x) -logPost(x, ...)
    optobj <- optim(start, fn, method = optimMethod,
                    control = optimctrl)
    if(optobj$convergence > 0)
      stop("failure in optimization")
    mode <- optobj$par
  }
  if(length(mode) > 20)
    stop("BAYESPACK only works for dimensions <= 20")
  con <- list(relreq = 0.01, maxvls=10000, numtrn = "split-t",
              df = 4, print = FALSE)
  mt <- c("Monte Carlo", "Lattice Rule", "Subregion-Adaptive",
          "Mixed Spherical-Radial", "Gauss-Hermite",
          "Modified Gauss-Hermite", "Stochastic Radial-Spherical")
  con$method <- match.arg(method, mt)
  stopifnot(names(control) %in% names(con))
  con[(namc <- names(control))] <- control
  if(!is.logical(con$print))
    stop("control$print needs to be a logical")
  con$mthd <- switch(con$method,
                     "Monte Carlo" = 0,
                     "Lattice Rule" = 1,
                     "Subregion-Adaptive" = 10,
                     "Mixed Spherical-Radial" = 20,
                     "Gauss-Hermite" = 30,
                     "Modified Gauss-Hermite" = 31,
                     "Stochastic Radial-Spherical" = 32,
                     #"Adaptive Radial-Spherical" = 33,
                     #"Monte Carlo Radial-Spherical" = 34,
                     #"Lattice Rule Radial-Spherical" = 35,
                     NA)
  if(is.na(con$mthd)) # should not happen
    stop("invalid integration method selected.")
  if(con$mthd == 10 | con$mthd == 33){
    d <- length(mode)
    if(con$mthd == 10)
      num <- 1+4*d+6*d^2+4*d*(d-1)*(d-2)/3+2^d
    if(con$mthd == 33)    
      num <- 1+2*d*(d+3)+2^d
    if(con$maxvls < ceiling(num*(2^3-1))){
      mes <- paste("need at least maxvls =", ceiling(num*(2^3-1)),
            "for dim =", d,"and method", con$method)
      stop(mes)
    }
  }
  if(con$mthd == 30){
    d <- length(mode)
    if(con$maxvls < 1+2^d){
      mes <- paste("need at least maxvls =", 1+2^d,
                   "for dim =", d,"and method", con$method)
      stop(mes)
    }
  }
  trans <- match.arg(trans, c("none", "mvt", "mvn", "split-t"))
  con$numtrn <- trans
  if(con$numtrn == "mvt"){
    if(con$df > 9 | con$df < 1)
      stop("invalid df argument selected!")
    trans <- paste("mvt (df=", con$df,")", sep="")
  }
  
  con$numtrn <- switch(con$numtrn,
                       "none" = 0,
                       "mvt" = con$df,
                       "mvn" = 10,
                       "split-t" = 20,
                       NA)
  if(con$mthd >= 20 & con$numtrn < 10){ # methods > 30 do not need standardization
    con$numtrn <- 10
    trans <- "mvn"
  }
  if(is.na(con$numtrn))
    stop("invalid transformation selected.")
  m <- length(mode)
  mn <- length(mnFns(mode))+m

  lgpst <- quote(logPost(.par, ...))
  mns <- quote(mnFns(.par, ...))

  rho <- new.env()
  mode2 <- numeric(m)
  for(i in 1:m)
    mode2[i] <- mode[i]
  assign(".m",      as.integer(m)          ,envir = rho) 
  assign(".mn",     as.integer(mn)         ,envir = rho)
  assign(".mu",     as.double(mode2)       ,envir = rho)
  assign(".par",    double(length(mode))   ,envir = rho)
  assign(".means",  double(mn)             ,envir = rho)
  assign(".errors", double(mn)             ,envir = rho)
  assign(".covrnc", double(m*(m+1)/2)      ,envir = rho)
  assign(".covmod", double(m*(m+1)/2)      ,envir = rho)
  assign(".relreq", as.double(con$relreq)  ,envir = rho)
  assign(".maxvls", as.integer(con$maxvls) ,envir = rho)
  assign(".numtrn", as.integer(con$numtrn) ,envir = rho)
  assign(".method", as.integer(con$mthd)   ,envir = rho)
  assign(".nrmcon", double(1)              ,envir = rho)
  assign(".inform", integer(1)             ,envir = rho)
  assign(".pru", as.integer(con$print)     ,envir = rho)

  .Call("BAYPOST", lgpst, mns, rho)

  nrmcon <- get(".nrmcon", rho)
  if(nrmcon <= 0){
    stop("calculated normalization constant <= 0!")
  } else {
    out <- list()
    out$method <- con$method
    out$trans <- trans
    out$nams <- names(mode)
    out$m <- m;out$mn <- mn
    out$mode <- get(".mu", rho)
    out$max <- logPost(out$mode, ...)
    out$nrmcon <- log(nrmcon)+out$max
    relreq <- get(".relreq", rho)
    log2 <- function(x) ifelse(x>0, log(x), -Inf)
    attr(out$nrmcon, "Int.Error Bounds") <- c(log2(nrmcon-relreq),log2(nrmcon+relreq))+out$max
    out$neval <- get(".maxvls", rho)
    out$inform <- get(".inform", rho)
    out$means <- get(".means", rho)
    out$errors <- get(".errors", rho)
    out$cov <- fillMat(get(".covrnc", rho))
    if(any(diag(out$cov) <= 0)){
      out$cov <- NULL
      warning("negative variances calculated by BAYESPACK")
    }
    mc <- fillMat(get(".covmod", rho))
    mc[upper.tri(mc)] <- 0
    out$modcov <- tcrossprod(mc)
  }
  class(out) <- "bayespack"
  out
}

print.bayespack <- function(x, digits = 4, ...){
  cat("BAYESPACK\n\n")
  cat("Log Norm. Constant:", round(x$nrmcon, 4), "\n")
  out <- c("yes", "no!")[x$inform+1]
  cat("Relative Error tolerance:", out, "\n\n")
  m <- x$m
  cat("Parameters:","\n")
  if(is.null(x$cov)){
    stdv <- rep(NA, m)
  } else {
    stdv <- formatC(sqrt(diag(x$cov)), digits = digits)
  }
  out <- data.frame("Means" = formatC(x$means[1:m], digits = digits),
                    "Stand Dev" = stdv)
  rownames(out) <- x$nams
  print(out)
}

summary.bayespack <- function(object, digits = 4,...){
  class(object) <- "summary.bayespack"
  print(object, digits = digits)
}

print.summary.bayespack <- function(x, digits = 4,...){
  cat("BAYESPACK\n\n")

  cat("Integration algorithm:", x$method, "\n")
  cat("Transformation:", x$trans, "\n")
  cat("No. of function evaluations:", x$neval, "\n")
  nc <- round(x$nrmcon, digits)
  eb <- round(attr(nc, "Int.Error Bounds"), digits)
  cat("Log Norm. Constant: ", nc, "\n", sep="")
  cat("Int.Error Bounds: (", eb[1], ",", eb[2] ,")\n", sep="")
  out <- c("yes", "no!")[x$inform+1]
  cat("Relative error tolerance:", out, "\n\n")
  m <- x$m
  cat("Parameter summary\n")
  if(is.null(x$cov)){
    stdv <- rep(NA, m)
  } else {
    stdv <- formatC(sqrt(diag(x$cov)), digits = digits)
  }
  out <- data.frame("Means" = formatC(x$means[1:m], digits = digits),
                    "Int Errors" = formatC(x$errors[1:m], digits = 2),
                    "Stand Dev" = stdv)
  nams <- x$nams
  if(is.null(nams))
    nams <- 1:m
  rownames(out) <- nams
  print(out)

  if(!is.null(x$cov)){
    cat(" Correlations:\n")
    corMat <- cov2cor(x$cov)
    dimnames(corMat) <- list(nams, nams)
    corMat <- round(corMat, digits)
    cf <- format(corMat, digits = digits)
    cf[row(cf) > col(cf)] <- ""
    print(cf, quote = FALSE)
  }

  mn <- x$mn
  if(mn > m){
    cat("\nAdditional means summary\n")
    out <- data.frame("Add Means" = formatC(x$means[(m+1):mn],
                      digits=digits), "Int Errors" =
                      formatC(x$errors[(m+1):mn], digits=2))
    print(out)
  }
}


fillMat <- function(x){
  d <- length(x)
  dm <- sqrt(8*d-1)/2
  mat <- matrix(0, dm, dm)
  z <- 1
  for(i in 1:dm){
    for(j in 1:i){
      mat[i,j] <- mat[j,i] <- x[z]
      z <- z+1
    }
  }
  mat
}
