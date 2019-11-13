### d/p/rStudent() #############################################################

##' @title Density of the Multivariate Student t Distribution
##' @param x (n, d)-matrix of evaluation points
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix, positive definite (scale != covariance matrix here)
##' @param factor *lower triangular* factor R of the covariance matrix 'scale'
##'        such that R^T R = 'scale' here (otherwise det(scale) not computed
##'        correctly!)
##' @param log logical indicating whether the logarithmic density is computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @param ... additional arguments passed to the underlying dnvmix()
##' @return n-vector of t_nu(loc, scale) density values
##' @author Erik Hintz and Marius Hofert
dStudent <- function(x, df, loc = rep(0, d), scale = diag(d),
                     factor = NULL, # needs to be triangular!
                     log = FALSE, verbose = TRUE, ...)
{
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x) # for 'loc', 'scale'
    dnvmix(x, qmix = "inverse.gamma", loc = loc, scale = scale,
           factor = factor, log = log, verbose = verbose, df = df, ...)
}

##' @title Distribution Function of the Multivariate Student t Distribution
##' @param upper d-vector of upper evaluation limits
##' @param lower d-vector of lower evaluation limits
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param method character string indicating the method to be used:
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         - "prng":    pure Monte Carlo
##' @param precond logical; if TRUE (recommended), variable reordering
##'        as described in Genz and Bretz (2002, pp. 955--956) is performed.
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
##' 


pStudent <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d),
                     df, loc = rep(0, d), scale = diag(d), standardized = FALSE,
                     control = list(), verbose = TRUE)
{
   ## Checks (needed to get the default for 'lower' correctly)
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   
    pnvmix(upper, lower = lower, qmix = "inverse.gamma", loc = loc, scale = scale,
           standardized = standardized, control = control,
           verbose = verbose, df = df)
}

##' @title Random Number Generator for the Multivariate Student t Distribution
##' @param n sample size
##' @param df degrees of freedom > 0; if df = Inf, sample from a Normal dist'n is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param factor factor R of the covariance matrix 'scale' with d rows
##'        such that R R^T = 'scale'.
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rStudent <- function(n, df, loc = rep(0, d), scale = diag(2),
                     factor = NULL, method = c("PRNG", "sobol", "ghalton"), skip = 0)
{
    d <- if(!is.null(factor)) { # for 'loc', 'scale'
             nrow(factor <- as.matrix(factor))
         } else {
             nrow(scale <- as.matrix(scale))
         }
    rnvmix(n, qmix = "inverse.gamma", rmix = "inverse.gamma",
           loc = loc, scale = scale, factor = factor, df = df,
           method = method, skip = skip)
}
