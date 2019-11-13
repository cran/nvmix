### d/p/rnvmixcop() ############################################################


##' Density function of a Multivariate Normal Variance Mixture Copula
##' @param u (n,d) matrix of evaluation points. Have to be in (0,1)
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale matrix).
##' @param factor Cholesky factor (lower triangular matrix) of 'scale'
##' @param control see ?pnvmixcop()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @param log logical indicating whether. the logarithmic density is to be computed
##' @param ... see ?pnvmix
##' @author Erik Hintz and Marius Hofert
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)

dnvmixcop <- function(u, qmix, scale = diag(d), factor = NULL, control = list(), 
                      verbose = FALSE, log = FALSE, ...){
  
  ## Most arguments are checked by qnvmix() and pnvmix()
  if(!is.matrix(u)) u <- rbind(u)
  d <- ncol(u) # for 'scale'
  
  ## Change accuracy for logdensity in qnvmix() to the one from dnvmix() here
  ## if the user did not provide a different one 
  ## (The default for abstol.newton.logdensity is chosen somewhat large for
  ## efficiency reasons as the logdensity there is only needed for Newton)
  names.control <- names(control)
  if(!any(names.control == "newton.logdens.abstol")){
    ## 'newton.logdens.abstol' was *not* provided:
    control <- get.set.parameters(control)
    control$newton.logdens.abstol <- control$dnvmix.abstol 
  }
  ## If it was provided, we don't change it. 
  
  ## Obtain quantiles. Note that qnvmix() takes in and returns a vector
  qu <- qnvmix(as.vector(u), qmix = qmix, control = control, 
               verbose = verbose, q.only =  FALSE, ...) # length n*d 
  ## log f_{X, scale} (F_{X1}^{-1}(u_{j1}),...,F_X1 ^{-1}(u_{jd})), j = 1,...,n
  num <- dnvmix(matrix(qu$q, ncol = d), qmix = qmix, scale = scale, factor = factor,
                control = control, verbose = verbose, log = TRUE, ...)# length n
  
  ## sum_{i=1}^d log f_{X1}( F_{X1}^{-1}(u_{ji})), j = 1,..,n
  ## Note that the log-density values are already calculated by qnvmix() 
  denom <- rowSums(matrix(qu$log.density, ncol = d)) # length n
  
  if(log) num - denom else exp(num - denom)
}



##' Distribution function of a Multivariate Normal Variance Mixture Copula
##' @param u (n,d) matrix of evaluation points. Have to be in (0,1)
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale matrix).
##' @param control see ?pnvmixcop()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @author Erik Hintz and Marius Hofert
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)

pnvmixcop <- function(u, qmix, scale = diag(d), control = list(), 
                      verbose = FALSE, ...){
  
  ## Most arguments are checked by qnvmix() and pnvmix()
  if(!is.matrix(u)) u <- rbind(u)
  d <- ncol(u) # for 'scale'
  ## Obtain quantiles. Note that qnvmix() returns a vector
  qu <- matrix(qnvmix(as.vector(u), qmix = qmix, control = control, 
                      verbose = verbose, q.only = TRUE, ...), ncol = d)
  pnvmix(qu, qmix = qmix, scale = scale, control = control, 
         verbose = verbose, ...)
}


##' Random Number Generator for Multivariate Normal Variance Mixtures
##' @param n sample size
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param factor (d, k)-matrix such that factor %*% t(factor) = scale;
##'        internally determined via chol() (and then an upper triangular
##'        matrix) if not provided
##' @param method character string indicating the method to be used:
##'         - "PRNG":    pure Monte Carlo
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         Note: For the methods "sobol" and "ghalotn", qmix() must be provided
##'         and rmix() is ignored. For the method "PRNG", either qmix() or rmix()
##'         needs to be provided. If both are provided, qmix() is ignored and
##'         rmix() is used.
##' @param skip numeric integer. How many points should be skipped when method='sobol'?
##' @param control see ?pnvmixcop()
##' @param verbose indicating whether a warning is given if the required precision
##'        in the underlying pnvmix() is not reached. 
##' @author Erik Hintz and Marius Hofert
##' @return (n, d)-matrix with NVM(0, scale)-copula samples 

rnvmixcop <- function(n, qmix, scale = diag(2), factor = NULL,
                      method = c("PRNG", "sobol", "ghalton"), skip = 0, 
                      control = list(), verbose = FALSE, ...)
{               
  d <- dim(scale)[1]
  stopifnot(dim(scale) == c(d,d))
  scale <- cov2cor(scale) # only need correlation matrix
  
  ## Sample from the nvmix dist'n
  sample.nvmix <- rnvmix(n = n, qmix = qmix, scale = scale, factor = factor,
                method = method, skip = skip, ...)
  ## Apply univariate margins
  ## Need (n,1) matrix as input so that pnvmix() gets the dimension right:
  sample.nvmixcop <- pnvmix(upper = matrix(sample.nvmix, ncol = 1), qmix = qmix, 
                            scale = matrix(1), standardized = TRUE, 
                            control = control, verbose = verbose, ...)
  ## Get dimensions correct and return
  matrix(sample.nvmixcop, ncol = d)
}