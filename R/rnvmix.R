### rnvmix() ###################################################################

##' @title Random Number Generator for Multivariate Normal Variance Mixtures
##' @param n sample size
##' @param rmix specification of random number generator of the  (mixture) distribution
##'        of W. This can be:
##'        1) character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "r", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) function being interpreted as a random number generator of W.
##'           Additional arguments can be passed via '...'
##'        4) n-vector containing a random sample from W.
##' @param qmix specification of the quantile function of the  (mixture) distribution
##'        of W. This needs to be supplied for the methods "sobol" and "ghalton".This can be:
##'        1) character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "qmix" random number generator.
##'        3) function being interpreted as the quantile function F_W^-.
##'           Additional arguments can be passed via '...'
##' @param loc d-vector (location != mean vector here)
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
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return (n, d)-matrix with NVM(loc,scale, F_W) samples
##' @author Marius Hofert and Erik Hintz
##' @note - For the Student t distribution, W ~ df/rchisq(n, df = df) but
##'         rchisq() simply calls rgamma(); see ./src/nmath/rchisq.c
##'         => W ~ 1/rgamma(n, shape = df/2, rate = df/2)
##'       - For a generalized inverse Gaussian distribution one could use:
##'         + "Runuran": faster if n large and parameters fixed; based on density
##'         + "GIGrvg":  faster if n small and often called with several parameters
##'         see examples of 'GIGrvg' for both methods
rnvmix <- function(n, rmix = NULL, qmix = NULL, loc = rep(0, d), scale = diag(2),
                   factor = NULL, method = c("PRNG", "sobol", "ghalton"),
                   skip = 0, ...)
{
    ## Basic checks
    stopifnot(n >= 1)
    method <- match.arg(method)

    ## Dealing with 'factor' (more general here than in dnvmix() and pnvmix1())
    if(is.null(factor)) { # => let 'factor' (internally here) be an *upper* triangular matrix
        factor <- chol(scale) # *upper* triangular; by this we avoid t() internally here and below around Z
        d <- nrow(factor)
        k <- d # => factor a square matrix here
    } else { # => 'factor' is a provided (d, k)-matrix (factor %*% factor^T = (d, d)-scale)
        d <- nrow(factor <- as.matrix(factor)) # (d, k)-matrix here...
        k <- ncol(factor)
        factor <- t(factor) # ... converted to a (k, d)-matrix to avoid t() below around Z
    }

    ## Determine if inversion is to be used
    inversion <- FALSE
    ## This is the case if the method used is "sobol" or "ghalton"
    if(method != "PRNG") inversion <- TRUE
    ## Or if the method is  "PRNQ" but rmix was not provivded
    if(method == "PRNG" && is.null(rmix)) inversion <- TRUE

    ## Get realizations of W in each case.
    if(inversion) { # work with 'qmix'

        ## In this case, we need qmix to be provided and use inversion
        if(is.null(qmix)) stop("'qmix' needs to be provided for methods 'sobol' and 'ghalton'")

        ## Get low discrepancy pointset
        ## Note that we need a k dimensional normal dist'n later, hence need k+1 uniforms
        ## since the first one is being used to get realizations of W
        U <- switch(method,
                    "sobol" = {
                        qrng::sobol(n, d = k + 1, randomize = TRUE, skip = skip)
                    },
                    "ghalton" = {
                        qrng::ghalton(n, d = k + 1, method = "generalized")
                    },
                    "PRNG" = {
                        matrix(runif( n* (k + 1)), ncol = k + 1)
                    })  # (n, k+1) matrix

        ## get quasi realizations using qmix
        W <- if(is.character(qmix)) { # 'qmix' is a character vector specifying supported mixture distributions (utilizing '...')
                 qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
                 switch(qmix,
                        "constant" = {
                            rep(1, n)
                        },
                        "inverse.gamma" = {
                            if(hasArg(df)) {
                                df <- list(...)$df
                            } else {
                                stop("'qmix = \"inverse.gamma\"' requires 'df' to be provided.")
                            }
                            ## Still allow df = Inf (normal distribution)
                            stopifnot(is.numeric(df), length(df) == 1, df > 0)
                            if(is.finite(df)) {
                                df2 <- df/2
                                1 / qgamma(1 - U[,1], shape = df2, rate = df2)
                            } else {
                                rep(1, n)
                            }
                        },
                        stop("Currently unsupported 'qmix'"))
             } else if(is.list(qmix)) { # 'qmix' is a list of the form (<character string>, <parameters>)
                 stopifnot(length(qmix) >= 1, is.character(distr <- qmix[[1]]))
                 qmix. <- paste0("q", distr)
                 if(!existsFunction(qmix.))
                     stop("No function named '", qmix., "'.")
                 do.call(qmix., append(list(U[,1]), qmix[-1]))
             } else if(is.function(qmix)) { # 'qmix' is interpreted as the quantile function F_W^- of the mixture distribution F_W of W
                 qmix(U[,1],...)
             } else stop("'qmix' must be a character string, list, quantile function or n-vector of non-negative random variates.")

    } else { # work with 'rmix'

        ## In this case, method == "PRNG" and we do not need inversion.
        ## Get realizations from rmix:
        if(is.null(rmix)) stop("'rmix() needs to be provided.")
        W <- if(is.character(rmix)) { # 'rmix' is a character vector specifying supported mixture distributions (utilizing '...')
                 rmix <- match.arg(rmix, choices = c("constant", "inverse.gamma"))
                 switch(rmix,
                        "constant" = {
                            rep(1, n)
                        },
                        "inverse.gamma" = {
                            if(hasArg(df)) {
                                df <- list(...)$df
                            } else {
                                stop("'rmix = \"inverse.gamma\"' requires 'df' to be provided.")
                            }
                            ## Still allow df = Inf (normal distribution)
                            stopifnot(is.numeric(df), length(df) == 1, df > 0)
                            if(is.finite(df)) {
                                df2 <- df/2
                                1 / rgamma(n, shape = df2, rate = df2)
                            } else {
                                rep(1, n)
                            }
                        },
                        stop("Currently unsupported 'rmix'"))
             } else if(is.list(rmix)) { # 'rmix' is a list of the form (<character string>, <parameters>)
                 stopifnot(length(rmix) >= 1, is.character(distr <- rmix[[1]]))
                 mix <- paste0("r", distr)
                 if(!existsFunction(mix))
                     stop("No function named '", rmix, "'.")
                 do.call(mix, c(n, rmix[-1]))
             } else if(is.function(rmix)) {
                 rmix(n, ...)# 'rmix' is interpreted as a random number generator for W
             } else if(is.numeric(rmix) && length(rmix) == n && all(rmix >= 0)) { # 'mix' is the vector of realizations of W
                 rmix
             } else stop("'rmix' must be a character string, list, random number generator or n-vector of non-negative random variates.")

    }

    ## Generate Z ~ N(0, I_k)
    Z <- if(!inversion){
        matrix(rnorm(n * k), ncol = k) # (n, k)-matrix of N(0, 1)
    } else {
        qnorm(U[, 2:(k+1)]) # (n, k)-matrix of N(0, 1)
    }
    ## Generate Y ~ N_d(0, scale)
    ## Recall that factor had been transposed, i.e. factor is (k,d)
    Y <- Z %*% factor # (n, k) %*% (k, d) = (n, d)-matrix of N(0, scale)
    ## Generate X ~ M_d(0, Sigma, LS[F_W])
    X <- sqrt(W) * Y # also fine for different k
    ## Generate X ~ M_d(mu, Sigma, LS[F_W])
    sweep(X, 2, loc, "+")
}
