### d/p/rNorm() ################################################################

##' @title Density of the Multivariate Normal Distribution
##' @param x (n, d)-matrix of evaluation points
##' @param loc d-vector (location = mean vector here)
##' @param scale (d, d)-covariance matrix, positive definite (scale = covariance matrix here)
##' @param factor *lower triangular* factor R of the covariance matrix 'scale'
##'        such that R^T R = 'scale' here (otherwise det(scale) not computed
##'        correctly!)
##' @param log logical indicating whether the logarithmic density is computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @param ... additional arguments passed to the underlying dnvmix()
##' @return n-vector of N(loc, scale) density values
##' @author Erik Hintz and Marius Hofert
dNorm <- function(x, loc = rep(0, d), scale = diag(d),
                  factor = NULL, # needs to be triangular!
                  log = FALSE, verbose = TRUE, ...)
{
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x) # for 'loc', 'scale'
    dnvmix(x, qmix = "constant", loc = loc, scale = scale,
           factor = factor, log = log, verbose = verbose, ...)
}

##' @title Distribution Function of the Multivariate Normal Distribution
##' @param upper d-vector of upper evaluation limits
##' @param lower d-vector of lower evaluation limits
##' @param loc d-vector (location = mean vector here)
##' @param scale (d, d)-covariance matrix (scale = covariance matrix here)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param method character string indicating the method to be used:
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         - "prng":    pure Monte Carlo
##' @param precond logical; if TRUE (recommended), variable reordering
##'        similar to Genz and Bretz (2002, pp. 955--956) is performed.
##'        Variable reordering can lead to a significant variance reduction
##'        and decrease in computational time.
##' @param abstol numeric >= 0 providing the absolute precision required.
##'        If abstol = 0, algorithm will run until total number of function
##'        evaluations exceeds fun.eval[2].
##' @param CI.factor Monte Carlo confidence interval multiplier. Algorithm runs
##'        CI.factor * (estimated standard error) < abstol. If CI.factor = 3.3
##'        (default), one can expect the actual absolute error to be less than
##'        abstol in 99.9% of the cases
##' @param fun.eval 2-vector giving the initial function evaluations (in the
##'        first loop; typically powers of 2) and the maximal number of
##'        function evaluations
##' @param B number of randomizations to get error estimates.
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)
##' @author Erik Hintz and Marius Hofert
pNorm <- function(upper, lower = rep(-Inf, d),
                  loc = rep(0, d), scale = diag(d), standardized = FALSE,
                  method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
                  abstol = 1e-3, CI.factor = 3.3, fun.eval = c(2^6, 1e8), B = 12,
                  verbose = TRUE)
{
    d <- length(upper) # for 'lower', 'loc', 'scale'
    pnvmix(upper, lower = lower, qmix = "constant", loc = loc, scale = scale,
           standardized = standardized, method = method, precond = precond,
           abstol = abstol, CI.factor = CI.factor, fun.eval = fun.eval, B = B,
           verbose = verbose)
}

##' @title Random Number Generator for the Multivariate Normal Distribution
##' @param n sample size
##' @param loc d-vector (location = mean vector here)
##' @param scale (d, d)-covariance matrix (scale = covariance matrix here)
##' @param factor factor R of the covariance matrix 'scale' with d rows
##'        such that R R^T = 'scale'.
##' @return (n, d)-matrix with N(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rNorm <- function(n, loc = rep(0, d), scale = diag(2),
                  factor = NULL, method = c("PRNG", "sobol", "ghalton"), skip = 0) # needs to be triangular!
{
    d <- if(!is.null(factor)) { # for 'loc', 'scale'
             nrow(factor <- as.matrix(factor))
         } else {
             nrow(scale <- as.matrix(scale))
         }
    rnvmix(n, qmix = "constant", rmix = "constant",
           loc = loc, scale = scale, factor = factor,
           method = method, skip = skip)
}
