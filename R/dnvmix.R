### dnvmix() ###################################################################

##' @title Density of a Multivariate Normal Variance Mixture
##' @param x (n, d)-matrix of evaluation points
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param loc d-vector (location vector)
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param factor Cholesky factor (lower triangular matrix) of 'scale';
##'        important here so that det(scale) is computed correctly!
##' @param method character string indicating the method to be used:
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         - "prng":    pure Monte Carlo
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
##' @param log logical indicating whether the logarithmic density is to be computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return n-vector with computed density values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert
dnvmix <- function(x, qmix, loc = rep(0, d), scale = diag(d),
                   factor = NULL, # needs to be lower triangular!
                   method = c("sobol", "ghalton", "PRNG"),
                   abstol = 0.001, CI.factor = 3.3, fun.eval = c(2^6, 1e8), B = 12,
                   log = FALSE, verbose = TRUE, ...)
{
    ## Checks
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x) # dimension
    if(!is.matrix(scale)) scale <- as.matrix(scale)
    stopifnot(length(loc) == d, dim(scale) == c(d, d), # note: 'qmix' is tested later
              abstol >= 0, CI.factor >= 0, length(fun.eval) == 2, fun.eval >= 0, B >= 1,
              is.logical(log))
    method <- match.arg(method)

    ## If factor is not provided, determine it here as a *lower* triangular matrix
    if(is.null(factor)) factor <- t(chol(scale)) # lower triangular

    ## 1 Define the quantile function of the mixing variable ###################

    ## If 'mix' is "constant" or "inverse.gamma", we use the analytical formulas
    is.const.mix <- FALSE # logical indicating whether we have a multivariate normal
    inv.gam <- FALSE # logical indicating whether we have a multivariate t
    qW <- if(is.character(qmix)) { # 'qmix' is a character vector
              qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
              switch(qmix,
                     "constant" = {
                         is.const.mix <- TRUE
                         function(u) 1
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
                             inv.gam <- TRUE
                             df2 <- df / 2
                             mean.sqrt.mix <- sqrt(df) * gamma(df2) / (sqrt(2) * gamma((df+1)/2)) # used for preconditioning
                             function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
                         } else {
                             is.const.mix <- TRUE
                             mean.sqrt.mix <- 1 # used for preconditioning
                             function(u) 1
                         }
                     },
                     stop("Currently unsupported 'qmix'"))
          } else if(is.list(qmix)) { # 'mix' is a list of the form (<character string>, <parameters>)
              stopifnot(length(qmix) >= 1, is.character(distr <- qmix[[1]]))
              qmix. <- paste0("q", distr)
              if(!existsFunction(qmix.))
                  stop("No function named '", qmix., "'.")
              function(u)
                  do.call(qmix., append(list(u), qmix[-1]))
          } else if(is.function(qmix)) { # 'mix' is the quantile function F_W^- of F_W
              function(u)
                  qmix(u, ...)
          } else stop("'qmix' must be a character string, list or quantile function.")

    ## Build result object (log-density)
    n <- nrow(x)
    lres <- rep(-Inf, n) # n-vector of results
    notNA <- rowSums(is.na(x)) == 0
    lres[!notNA] <- NA
    x <- x[notNA,] # non-missing data (rows)

    ## 2 Actual computation ####################################################

    ## Recall that 'scale' is *lower triangular*. For short, let 'scale' = L
    ## Solve L * z = x_i - mu for z, so z = L^{-1} * (x_i - mu)   (d vector)
    ## => z^2 (=> componentwise) = z^T z = (x_i - mu)^T * (L^{-1})^T L^{-1} (x_i - mu)
    ##                           = z^T z = (x_i - mu)^T * (L L^T )^{-1} (x_i - mu)
    ##                           = (x_i - mu)^T * scale^{-1} * (x_i - mu)
    ##                           = quadratic form
    ## Now do this for *all* x_i simultaneously using that L is lower triangular:
    ## Forwardsolve: "right-hand sides" of equation must be in the *columns*, thus t(x)
    z <- forwardsolve(factor, t(x) - loc, transpose = FALSE)
    maha2 <- colSums(z^2) # = sum(z^T z); n-vector of squared Mahalanobis distances from x to mu w.r.t. scale
    ## Note: could probably be done with mahalanobis() but unclear how we would
    ##       get det(scale) then.

    ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
    ##                       = log(prod(diag(R))) = sum(log(diag(R)))
    lrdet <- sum(log(diag(factor)))
    if(!is.finite(lrdet)) stop(paste("Density not defined for singular 'scale' "))

    ## Counters
    total.fun.evals <- 0 # total.fun.evals will count the total number of function evaluations
    numiter <- 0 # initialize counter; this will count the number of iterations in the while loop

    ## Deal with the different distributions
    if(inv.gam) { # multivariate t
        df.d.2 <- (df + d) / 2
        lres[notNA] <- lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) -
            lrdet - df.d.2 * log1p(maha2 / df)
        error <- 0
    } else if(is.const.mix) { # multivariate normal
        lres[notNA] <- -(d/2) * log(2 * pi) - lrdet - maha2/2
        error <- 0
    } else { # general case of a multivariate normal variance mixture (RQMC)

        ## Basics
        CI.factor <- CI.factor / sqrt(B) # instead of dividing sigma by sqrt(B) each time
        current.n <- fun.eval[1] # initial n
        rqmc.estimates <- matrix(0, ncol = n, nrow = B) # matrix to store RQMC estimates
        error <- abstol + 42 # initialize error to > abstol to enter while loop
        useskip <- 0 # needed so that first iteration differs from the others
        denom <- 1

        ## Make sure seed exists for 'method' being "sobol"
        if(method == "sobol") {
            if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
            seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
        }

        ## 3 Main loop #########################################################

        while(error > abstol && total.fun.evals < fun.eval[2]) {

            if(method == "sobol") .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )

            ## Get B RQCM estimates
            for(l in 1:B) {
                ## Get the point set
                U <- switch(method,
                            "sobol"   = {
                                qrng::sobol(current.n, d = 1, randomize = TRUE, skip = (useskip * current.n))
                            },
                            "gHalton" = {
                                qrng::ghalton(current.n, d = 1, method = "generalized")
                            },
                            "prng"    = {
                                cbind(runif(current.n)) # 1-column matrix
                            })
                W <- qW(U) # current.n-vector of W's

                ## exp-log trick
                b <- - (d/2) * log(2 * pi * W) - lrdet - outer(1/W, maha2 / 2) # (current.n, n)-matrix, each column corresponds to "one x"
                bmax <- apply(b, 2, max) # n-vector of maximal b's
                rqmc.estimates[l,] <- (rqmc.estimates[l,] - log(current.n) + bmax +
                                       log(colSums(exp(b - rep(bmax, each = current.n))))) / denom
            }

            ## Update various variables
            total.fun.evals <- total.fun.evals + B * current.n

            ## Change denom and useskip; this is done exactly once, namely in the first iteration.
            if(numiter == 0){
                denom <- 2
                useskip <- 1

            } else {
                ## Increase the sample size n; this is done in all iterations except the first
                current.n <- 2 * current.n
            }

            ## Compute error measures and update counter
            sig <- max(apply(rqmc.estimates, 2, sd)) # get largest standard deviation (over all columns)
            error <- CI.factor * sig # update error; CI.factor here is actually CI.factor/sqrt(N)
            numiter <- numiter + 1 # update counter

            ## Update abserr to ensure that the precision is reached for *both*
            ## the log-density and density
            log.estimates <- colMeans(rqmc.estimates)
            abstol <- abstol / max(1, max(exp(log.estimates)))

        } # while()

        ## Finalize
        lres[notNA] <- log.estimates

        ## If 'log = FALSE' need to adjust error estimate.
        ## This is a conservative error estimate.
        if(!log) error <- error * max( exp(max( lres[notNA])), 1)

        ## Finalize
        if(verbose && (error > abstol))
            warning("'abstol' not reached; consider increasing 'fun.eval[2]'")
    }

    ## Return
    attr(lres, "error")   <- error
    attr(lres, "numiter") <- numiter
    if(log) lres else exp(lres)
}
