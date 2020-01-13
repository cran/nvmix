### fitnvmix() #################################################################

### Functions to estimate weights given 'nu', 'loc', 'scale' ###################

##' @title Estimate Weights for fitnvmix()
##' @param maha2.2 squared maha distances divided by 2 (length n)
##' @param qW see ?fitnvmix() ('qmix' there)
##' @param nu parameter (vector) nu of W
##' @param lrdet log(sqrt(det(scale)))
##' @param d dimension
##' @param special.mix either NA or string. Currently supported are 'inverse.gamma'
##'         and 'pareto' in which case analytical weights are calculated
##' @param control see ?fitnvmix()
##' @param verbose see ?fitnvmix()
##' @return List of three:
##'         $weights n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either relative
##'         error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms
##'         (only if return.all = TRUE)
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
##' @note corresponds to delta_ki in the paper
weights_ <- function(maha2.2, qW, nu, lrdet, d, special.mix, control, verbose)
{
    verbose <- as.logical(verbose) # only logical needed here
    if(!is.na(special.mix)) { # weights are known analytically
        weights <- switch(special.mix,
                          "inverse.gamma" = {
                              (nu + d) / (nu + maha2.2*2)
                          },
                          "pareto" = {
                              pgamma(1, shape = nu+d/2+1, scale = 1/maha2.2)/
                                  pgamma(1, shape = nu+d/2, scale = 1/maha2.2)*
                                  (nu + d/2)/maha2.2
                          })
        numiter <- 0
        error <- rep(0, length(maha2.2))
    } else { # weights need to be estimated
        ## Absolte/relative precision?
        if(is.na(control$weights.reltol)) {
            tol <- control$weights.abstol
            do.reltol <- FALSE
        } else {
            ## Use relative error
            tol <- control$weights.reltol
            do.reltol <- TRUE
        }
        ## lconst for the integrand
        lconst <- rep(-lrdet - d/2*log(2*pi), length(maha2.2))
        ## Call RQMC procedure without any stratification
        rqmc.obj <- weights_rqmc(maha2.2, qW = qW, nu = nu,
                                          lconst = lconst, d = d,
                                          max.iter.rqmc = control$dnvmix.max.iter.rqmc.pilot,
                                          control = control, return.all = TRUE)
        ## Extract results
        weights <- rqmc.obj$weights
        numiter <- rep(rqmc.obj$numiter, length(maha2.2))
        error   <- rqmc.obj$error
        if(any(error > tol)) {
            ## Accuracy not reached for at least one 'maha2.2' value
            ## => Use adaptive approach for those
            notRchd <- which(error > tol)
            qW. <- function(u) qW(u, nu = nu)
            ldens.obj <- densmix_adaptrqmc(qW., maha2.2 = maha2.2[notRchd],
                                           lconst = lconst[notRchd], d = d, 
                                           UsWs = rqmc.obj$UsWs, control = control)
            lcond.obj <- densmix_adaptrqmc(qW., maha2.2 = maha2.2[notRchd],
                                           lconst = lconst[notRchd], d = d, 
                                           k = d + 2, UsWs = rqmc.obj$UsWs,
                                           control = control)
            weights[notRchd]   <- exp(lcond.obj$ldensities - ldens.obj$ldensities)
            ## Which weights cannot be reliably estimated?
            which.errorNA     <- which(is.na(lcond.obj$error) | is.na(ldens.obj$error))
            if(any(which.errorNA)) {
                if(verbose) warning('Some weights cannot be reliably estimated')
                ## Check if 'weights' decreasing in 'maha2.2'
                n <- length(weights)
                for(i in 1:n) {
                    if(i <= 2) {
                        next
                    } else if(weights[i] <= weights[i-1]) {
                        next
                    } else {
                        ## In this case, weights[i] > weights[i-1]
                        ## Case 1: Is there any weight beyond i which is smaller than weights[i-1]?
                        smallerthani <- which(weights[i:n] <= weights[i-1]) + i - 1
                        if(length(smallerthani) > 0) {
                            ## If that's the case, interpolate all weights between i
                            ## and the first one smaller than weights[i-1]
                            firstsmaller <- smallerthani[1]
                            slope <- (weights[firstsmaller] - weights[i-1]) /
                                (maha2.2[firstsmaller] - maha2.2[i-1])
                            weights[i:(firstsmaller-1)] <- weights[i-1] + slope*
                                (maha2.2[i:(firstsmaller-1)] - maha2.2[i-1])
                        } else {
                            ## Behavior of the function suggests log-log extrapolation:
                            slope <- log(weights[i-1] / weights[i-2]) / log(maha2.2[i-1] / maha2.2[i-2])
                            weights[i:n] <- weights[i-1] * exp(slope*log(maha2.2[i:n] / maha2.2[i-1]))
                        }
                    }
                }
            }
        }
    }
    list(weights = weights, numiter = numiter, error = error)
}

##' @title Estimate Weights for fitnvmix()
##' @param maha2.2 squared maha distances divided by 2 (length n)
##' @param qW see ?fitnvmix() ('qmix' there)
##' @param nu parameter (vector) nu of W
##' @param lconst see ?densmix_
##' @param d dimension
##' @param control see ?fitnvmix()
##' @param max.iter.rqmc maximum number of RQMC iterations
##' @param return.all logical; if true, matrix (U, qW(U)) also returned.
##' @return List of three:
##'         $weights n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either relative
##'         error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms
##'         (only if return.all = TRUE)
##' @author Erik Hintz
##' @note corresponds to delta_ki in the paper
weights_rqmc <- function(maha2.2, qW, nu, lconst, d, max.iter.rqmc,
                                  control, return.all)
{
    ## Define various quantites:
    dblng           <- TRUE
    B               <- control$B # number of randomizations
    n               <- length(maha2.2) # sample size
    current.n       <- control$fun.eval[1] #initial sample size
    numiter         <- 0 # counter for the number of iterations
    ZERO            <- .Machine$double.neg.eps
    total.fun.evals <- 0
    ## Absolte/relative precision?
    if(is.na(control$weights.reltol)) {
        ## Use absolute error
        tol <- control$weights.abstol
        do.reltol <- FALSE
    } else {
        ## Use relative error
        tol <- control$weights.reltol
        do.reltol <- TRUE
    }
    ## Store seed if 'sobol' is used to get the same shifts later:
    if(control$method == "sobol") {
        useskip <- 0 # to get correct shifts if method = "sobol"
        if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
        seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
    }
    denom <- 1
    ## Matrix to store RQMC estimates
    rqmc.estimates.lweights <- matrix(-Inf, ncol = n, nrow = B)
    ## Will be needed a lot:
    CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
    ## Initialize 'max.error' to > tol so that we can enter the while loop:
    max.error <- tol + 42
    ## Matrix to store U, W values => nrows = maximal number of funevals
    if(return.all) {
        max.nrow <- current.n*B*2^(max.iter.rqmc-1)
        UsWs <- matrix(NA, ncol = 2, nrow = max.nrow)
        curr.lastrow <- 0 # will count row-index additional points are being inserted after
    }

    ## Main loop
    ## while() runs until precision abstol is reached or the number of function
    ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
    ## the desired log-densities are calculated.
    while(max.error > tol && numiter < max.iter.rqmc &&
          total.fun.evals < control$fun.eval[2])
    {
        ## Reset seed to have the same shifts in sobol(...)
        if(control$method == "sobol" && numiter > 0)
            .Random.seed <<- seed # reset seed to have the same shifts in sobol(...)
        for(b in 1:B) {
            ## Get the point set
            U <- sort(switch(control$method,
                             "sobol" = {
                                 if(dblng) {
                                     qrng::sobol(n = current.n, d = 1,
                                                 randomize = TRUE,
                                                 skip = (useskip * current.n))
                                 } else {
                                     qrng::sobol(n = current.n, d = 1,
                                                 randomize = TRUE,
                                                 skip = (numiter * current.n))
                                 }
                             },
                             "ghalton" = {
                                 qrng::ghalton(n = current.n, d = 1,
                                               method = "generalized")
                             },
                             "PRNG" = {
                                 runif(current.n)
                             })) # sorted for later!
            ## Evaluate the integrand at the (next) point set
            W <- qW(U, nu = nu) # realizations of the mixing variable; sorted!
            ## Need to replace values < ZERO by ZERO. W is *sorted*, so check using
            ## loop instd of 'pmax' (more efficient)
            for(ind in 1:current.n) if(W[ind] < ZERO) W[ind] <- ZERO else break
            ## Update 'UsWs'
            if(return.all) {
                UsWs[(curr.lastrow + 1) : (curr.lastrow + current.n), ] <- cbind(U, W)
                curr.lastrow <- curr.lastrow + current.n
            }
            next.est.condexp <- .Call("eval_densmix_integrand",
                                      W          = as.double(W),
                                      maha2_2    = as.double(maha2.2),
                                      current_n  = as.integer(current.n),
                                      n          = as.integer(n),
                                      d          = as.integer(d),
                                      k          = as.integer(d + 2), # k=d+2 here!
                                      lconst     = as.double(lconst))
            next.est.ldens <- .Call("eval_densmix_integrand",
                                    W          = as.double(W),
                                    maha2_2    = as.double(maha2.2),
                                    current_n  = as.integer(current.n),
                                    n          = as.integer(n),
                                    d          = as.integer(d),
                                    k          = as.integer(d), # k=d+2 here!
                                    lconst     = as.double(lconst))

            ## Update RQMC estimates
            rqmc.estimates.lweights[b, ] <-
                .Call("logsumexp2",
                      a = as.double(rqmc.estimates.lweights[b, ]),
                      b = as.double(next.est.condexp - next.est.ldens),
                      n = as.integer(n)) - log(denom)

        } # end for(b in 1:B)
        ## Update of various variables
        ## Double sample size and adjust denominator in averaging as well as useskip
        if(numiter == 0) {
            ## Change denom and useksip (exactly once, in the first iteration)
            denom <- 2
            useskip <- 1
        } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
        }
        ## Total number of function evaluations:
        total.fun.evals <- total.fun.evals + B * current.n
        numiter <- numiter + 1
        ## Update error for 'weights' (NOT log values!)
        rqmc.estimates.weights <- exp(rqmc.estimates.lweights)
        weights  <- colMeans(rqmc.estimates.weights)
        vars     <- .colMeans((rqmc.estimates.weights -
                               rep(weights, each = B))^2, B, n, 0)
        errors   <- if(!do.reltol) {
                        sqrt(vars)*CI.factor.sqrt.B
                    } else {
                        sqrt(vars)/abs(weights)*CI.factor.sqrt.B
                    }
        max.error <- max(errors)
    } # while()

    ## Return
    ret.obj <- if(return.all) {
                   list(weights = weights, numiter = numiter, error = errors,
                        UsWs = UsWs)
               } else {
                   list(weights = weights, numiter = numiter, error = errors)
               }
    ret.obj
}


### Function to estimate 'nu' given 'loc', 'scale'##############################

##' @title Estimate 'nu' Given 'loc', 'scale' by Maximizing Log-Likelihood
##' @param tx t(x) where x is as in ?fitnmvix()
##' @param qW quantile function of W; must be function(u, nu)
##' @param init.nu initial estimate of nu
##' @param loc current estimate of 'loc'
##' @param scale current estimate of 'scale'
##' @param factor cholesky factor of the 'scale'; if not provided, it's calculated
##' @param mix.param.bounds see ?fitnvmix()
##' @param special.mix string specifying if W has a special distribution for which
##'        weights are known; currently, only 'inverse.gamma' and 'pareto' supported
##' @param control see ?fitnvmix()
##' @param control.optim passed to optim; see ?optim
##' @param verbose see ?fitnvmix
##' @return list of two: $nu.est (scalar of vector of length 'init.nu'; MLE estimate of nu)
##'                      $max.ll (negative log-likelihood at nu.est)
##'                      $ll.counts (total number of calls to likelihood)
##'                      $opt.obj (object returned by underlying 'optim' call)
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
estim_nu <- function(tx, qW, init.nu, loc, scale, factor = NA, mix.param.bounds,
                     special.mix = NA, control, control.optim, verbose)
{
    ## Obtain 'factor' if not provided
    if(is.na(factor)) factor <- t(chol(scale))
    ll.counts <- 0 # counts number of calls to 'neg.log.likelihood.nu'
    if(is.character(special.mix)) {
        ## In this case, dnvmix() uses analytical formula for the density
        neg.log.likelihood.nu <-
            switch(special.mix,
                   "inverse.gamma" = {
                       function(nu) {
                           ll <- -sum(dnvmix(t(tx), qmix = "inverse.gamma",
                                             loc = loc, factor = factor, df = nu,
                                             log = TRUE, verbose = verbose))
                           if(verbose >= 3) cat(".") # print dot after each call to likelihood
                           ll.counts <<- ll.counts + 1 # update 'll.counts' in the parent environment
                           ll # return
                       }
                   },
                   "pareto" = {
                       function(nu) {
                           ll <- -sum(dnvmix(t(tx), qmix = "pareto",
                                             loc = loc, factor = factor, alpha = nu,
                                             log = TRUE, verbose = verbose))
                           if(verbose >= 3) cat(".") # print dot after each call to likelihood
                           ll.counts <<- ll.counts + 1 # update 'll.counts' in the parent environment
                           ll # return
                       }
                   })
    } else {
        ## Get various quantitites passed to 'densmix_'
        z        <- forwardsolve(factor, tx - loc, transpose = FALSE)
        maha2.2  <- sort(colSums(z^2)/2)
        lrdet    <- sum(log(diag(factor)))
        d        <- ncol(factor)
        ## Get and store current seed (=> same shifts in sobol)
        if(!exists(".Random.seed")) runif(1)
        seed <- .Random.seed
        ## Set up -loglikelihood as a function of 'nu'
        lconst <- rep(-lrdet - d/2*log(2*pi), length(maha2.2))
        neg.log.likelihood.nu <- function(nu) {
            .Random.seed <<- seed # reset seed => monotonicity (not bc of sobol shifts!)
            qmix. <- function(u) qW(u, nu = nu) # function of u only
            ## Call 'densmix_()' which by default returns the log-density
            ldens.obj <- densmix_(qW = qmix., maha2.2 = maha2.2, lconst = lconst,
                                          d = d, control = control, verbose = verbose)
            if(verbose >= 3) cat(".") # print dot after each call to likelihood
            ll.counts <<- ll.counts + 1 # update 'll.counts' in parent environment
            ## Return -log-density
            -sum(ldens.obj$ldensities)
        }
    }
    ## Optimize neg.log.likelihood over nu
    opt.obj <- optim(init.nu, fn = neg.log.likelihood.nu,
                     lower = mix.param.bounds[, 1],
                     upper = mix.param.bounds[, 2],
                     method = "L-BFGS-B", control = control.optim)
    list(nu.est    = opt.obj$par,
         max.ll    = opt.obj$value,
         ll.counts = ll.counts,
         opt.obj   = opt.obj) # also return full 'opt.obj'
}


### Main function ##############################################################

##' @title Fitting Multivariate Normal Variance Mixtures
##' @param x (n,d) data matrix
##' @param qmix character string ("constant", "inverse.gamma") or function. If
##'        function, it *has* to be qmix(u, nu)
##' @param mix.param.bounds either a vector of length two (in which case 'nu'
##'        is a scalar) or a matrix with two columns, where element in row i in
##'        1st/2nd column corresponds to lower/upper limits for component i of
##'        'nu'. All elements need to be finite, numeric values.
##'        (eg if nu/W~chi^2_nu, then eg 'mix.param.bounds = c(1, 10)'
##'        Note: The smaller the range, the better.
##' @param nu.init initial estimate for 'nu'; either NA in which case it is
##'        estimated or a vector of length = length of parameter vector 'nu'.
##'        If provided and close to MLE, can speed up 'fitnvmix' significantly
##' @param init.size.subsample if 'is.na(nu.init)', size of subsample of 'x'
##'        used to obtain initial estimate of nu.
##' @param size.subsample numeric, <= nrow(x). Number of rows of 'x' used in ECME
##'        iteration to optimize the log-likelihood. Defaults to n
##'        (all datapoints are used)
##' @param control list of algorithm specific parameters, see ?get_set_param
##'        and ?fitnvmix
##' @param verbose numeric or logical. 0: No warnings; 1: Warnings;
##'        2: Warnings + short tracing; 3: Warnings + complete tracing.
##' @return list of three (if qmix = "constant"), otherwise 5 or 7
##'         $nu: estimate for nu (omitted if qmix = "constant")
##'         $loc: estimate for the location vector
##'         $scale: estimate for scale matrix
##'         $iter: number of ECME iterations
##'         $max.ll: log-likelihood at (nu, loc, scale)
##'         $nu.Ests: matrix of estimates of 'nu' with corresponding loglikelihood
##'                   in each iteration (only if isTRUE(control.addRetruns))
##'         $iter.converged: iterations needed until convergence
##'                          (only if isTRUE(control.addRetruns))
##' TODO    include option to give names to parameters etc
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
fitnvmix <- function(x, qmix, mix.param.bounds, nu.init = NA,
                     init.size.subsample = min(n, 100),
                     size.subsample = n, control = list(),
                     verbose = TRUE)
{
    ## Basics
    if(!is.matrix(x))
        x <- rbind(x)
    ## Initialize various quantities
    control <- get_set_param(control)
    ## Get quantile function
    special.mix <- NA  # to record if we have a special dist'n (normal, t,...)
    ## Set up qW as function(u, nu)
    qW <- if(is.character(qmix)) {# 'qmix' is a character vector
              qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma", "pareto"))
              switch(qmix,
                     "constant" = {
                         special.mix <- "constant"
                         function(u, nu) rep(1, length(u))
                     },
                     "inverse.gamma" = {
                         special.mix <- "inverse.gamma"
                         function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
                     },
                     "pareto" = {
                         special.mix <- "pareto"
                         function(u, nu) (1-u)^(-1/nu)
                     },
                     stop("Currently unsupported 'qmix'"))
          } else if(is.list(qmix)) { # 'mix' is a list of the form (<character string>, <parameters>)
              ## Not supported yet
              stop("'qmix' cannot be a list when passed to 'fitnvmix'")
          } else if(is.function(qmix)) { # 'mix' is the quantile function F_W^- of F_W
              function(u, nu)
                  qmix(u, nu)
          } else stop("'qmix' must be a character string, list or quantile function.")

    ## Case of MVN: MLEs are sample mean and sample cov matrix
    if(is.character(special.mix) && special.mix == "constant") {
        loc.est <- colMeans(x)
        ## TODO: Do this better (as.matrix to remove attributes, can be done better)
        scale.est <- as.matrix(nearPD(cov(x))$mat) # sample covariance matrix
        return(list(loc = loc.est, scale = scale.est, iter = 0))
    }

    ## Check 'verbose' argument:
    if(is.logical(verbose)) {
        verbose <- as.integer(verbose)
    } else if(!(verbose %in% c(0, 1, 2, 3))) {
        stop("'verbose' has to be either logical or an integer between 0 and 3")
    }

    ## Check inputs, get dimensions
    notNA <- rowSums(is.na(x)) == 0
    x     <- x[notNA,, drop = FALSE] # non-missing data (rows)
    tx    <- t(x)
    n     <- nrow(x)
    d     <- ncol(x)

    ## Use only sub-sample to estimate 'nu'?
    if(size.subsample < n) {
        sampled.ind <- sample(n, size.subsample)
        x.sub       <- x[sampled.ind,, drop = FALSE]
        tx.sub      <- tx[, sampled.ind, drop = FALSE]
    } else {
        x.sub  <- x
        tx.sub <- tx
    }

    ## Check parameter bounds on 'nu' and get 'mix.param.length'
    mix.param.length <- if(is.vector(mix.param.bounds)) {
                            stopifnot(length(mix.param.bounds) == 2)
                            mix.param.bounds <- matrix(mix.param.bounds, nrow = 1)
                            1
                        } else if(is.matrix(mix.param.bounds)) {
                            stopifnot(dim(mix.param.bounds)[2] == 2)
                            dim(mix.param.bounds)[1]
                        } else stop("'mix.param.bounds' has to be either a vector of length 2 or a matrix with 2 columns")

    if(control$addReturns) {
        ## Matrix storing all 'nu' estimates with corresponding likelihoods
        nu.Ests <- matrix(NA, ncol = mix.param.length + 1,
                          nrow = (nrow <- if(control$laststep.do.nu) control$ECME.maxiter+2
                                          else control$ECME.maxiter + 1))
        rownames(nu.Ests) <- if(control$laststep.do.nu) {
                                 c("Initial",
                                   sapply(1:control$ECME.maxiter,
                                          function(i) paste0("ECME-iteration ", i)),
                                   "Laststep")
                             } else {
                                 c("Initial",
                                   sapply(1:control$ECME.maxiter,
                                          function(i) paste0("ECME-iteration ", i)))
                             }
        colnames(nu.Ests) <- c(sapply(1:mix.param.length,
                                      function(i) paste0("nu[", i, "]")),
                               "Log-likelihood")
        current.iter.total <- 1
    }

    ## 1 Initial estimates for nu, loc, scale ##################################

    ## Unbiased estimator for 'loc' based on full sample:
    loc.est <- colMeans(x)
    ## Sample covariance matrix based on full sample:
    SCov    <- as.matrix(nearPD(cov(x))$mat) # TODO do this smarter
    ## Check if 'init.nu' was provided. If so, calculate 'scale' as (1/E(W))*SCov
    if(!is.na(nu.init)) {
        stopifnot(length(nu.init) == mix.param.length)
        if(verbose >= 2) cat("Step 1: Initial estimate for 'nu': Was provided")
        nu.est <- nu.init
        scale.est <- 1/mean(qW(runif(1e4), nu.est))* SCov
        scale.est <- SCov
    } else {
        if(verbose >= 2) cat("Step 1: Initial estimate for 'nu' by optimizing log-likelihood ")
        ## Optimize log-likelihood:
        ll.counts <- 0 # counts number of calls to 'neg.log.likelihood.init'
        ## Get and store current seed (=> same shifts in sobol)
        if(!exists(".Random.seed")) runif(1)
        seed <- .Random.seed
        ## -loglikelihood as function of param=(nu, c) of length mix.param.length + 1

        neg.log.likelihood.init <-
            if(is.character(special.mix)) {
                switch(special.mix,
                       "inverse.gamma" = {
                           function(param) {
                               ll <- -sum(dnvmix(x[sample(n, init.size.subsample),, drop = FALSE],
                                                 qmix = "inverse.gamma", loc = loc.est,
                                                 scale = param[2] * SCov, df = param[1],
                                                 log = TRUE))
                               if(verbose >= 3) cat(".") # print dot after each call to likelihood
                               ll.counts <<- ll.counts + 1 # update 'll.counts' in the parent environment
                               ll # return
                           }
                       },
                       "pareto" = {
                           function(param) {
                               ll <- -sum(dnvmix(x[sample(n, init.size.subsample),, drop = FALSE],
                                                 qmix = "pareto", loc = loc.est,
                                                 scale = param[2] * SCov, alpha = param[1],
                                                 log = TRUE))
                               if(verbose >= 3) cat(".") # print dot after each call to likelihood
                               ll.counts <<- ll.counts + 1 # update 'll.counts' in the parent environment
                               ll # return
                           }
                       })
            } else {
                function(param) {
                    ## Define a 'qmix.' function of u only that can be passed to dnvmix():
                    .Random.seed <<- seed # for monotonicity
                    ## Return - loglikelihood
                    ll <- -sum(dnvmix(x[sample(n, init.size.subsample),, drop = FALSE],
                                      qmix = qW, loc = loc.est,
                                      scale = param[mix.param.length + 1] * SCov,
                                      control = control, verbose = verbose, log = TRUE,
                                      nu = param[1:mix.param.length]))
                    if(verbose >= 3) cat(".") # print dot after each call to likelihood
                    ll.counts <<- ll.counts + 1 # update ll.counts from parent environment
                    ll
                }
            }
        ## Optimize -log.likelihood over (nu = nu, scale = c*SCov), 'c' scalar
        ## Initial parameter: For 'nu', midpoint of feasible parameter set;
        ##  for 'c', 1/E(W) where E(W) estimated via RQMC taking 'nu = nu.init'
        init.param <- c(rowMeans(mix.param.bounds),
                        1/mean(qW(runif(1e4), nu = rowMeans(mix.param.bounds))))
        opt.obj <- optim(init.param, fn = neg.log.likelihood.init,
                         lower = c(mix.param.bounds[, 1], 0.1),
                         upper = c(mix.param.bounds[, 2], NA),
                         method = "L-BFGS-B", control = control$control.optim)
        ## Grab estimate for 'nu' as well as for the 'scale' matrix
        nu.est    <- opt.obj$par[1:mix.param.length]
        scale.est <- opt.obj$par[mix.param.length + 1] * SCov
        max.ll    <- opt.obj$value
        if(verbose >= 3)
            cat(paste0(" DONE (", ll.counts, " calls to likelihood needed)", '\n'))
    }
    ## Store if needed
    if(control$addReturns) {
        ## Matrix storing all 'nu' estimates with corresponding likelihoods
        ll <- sum(dnvmix(x, qmix = qW, loc  = loc.est, scale = scale.est, nu = nu.est,
                         control = control, verbose = verbose, log = TRUE))
        nu.Ests[current.iter.total, ] <- c(nu.est, ll)
        current.iter.total <- current.iter.total + 1
        iter.converged <- control$ECME.maxiter + 2
    }

    ## 2 ECME step #############################################################

    if(control$ECMEstep) {
        if(verbose == 2) cat(paste0('\n')) # if 'verbose==3' linebreak already happened
        if(verbose >= 2) cat(paste0("Step 2: ECME iteration.", '\n'))
        ## Initialize various quantities
        iter.ECME            <- 0
        converged            <- FALSE
        ## Main loop:
        while(iter.ECME < control$ECME.maxiter & !converged) {
            if(verbose >= 2) cat(paste0("  Iteration ",iter.ECME + 1, '\n'))

            ## 'loc.est' and 'scale.est' updates
            converged.locscale   <- FALSE
            iter.locscaleupdate  <- 1
            ## Update 'scale.est' and 'loc.est' given current estimate of 'nu.est'
            ## until convergence.
            if(verbose >= 3) cat(paste0("    Estimating weights and updating 'loc' and 'scale'"))

            ## Inner loop (iterating over 'loc' and 'scale' with 'nu.est' held fixed)
            while(!converged.locscale &&
                  iter.locscaleupdate < control$max.iter.locscaleupdate)
            {
                ## Get new 'weights'
                ## First, get new maha distances (with current 'loc.est' and 'scale.est')
                factor               <- t(chol(scale.est))
                lrdet                <- sum(log(diag(factor)))
                z                    <- forwardsolve(factor, tx - loc.est, transpose = FALSE) # use the full sample!
                maha2.2.new          <- colSums(z^2)/2
                order.maha2.2.new    <- order(maha2.2.new)
                maha2.2.new          <- maha2.2.new[order.maha2.2.new] # sorted increasingly
                if(iter.locscaleupdate == 1) {
                    ## Only in the first iteration do we approximate *all* weights by RQMC.
                    ## Get weights:
                    weights <- weights_(maha2.2.new, qW = qW, nu = nu.est, 
                                        lrdet = lrdet, d = d, 
                                        special.mix = special.mix, 
                                        control = control, verbose = verbose)$weights
                    weights.new <- weights[order(order.maha2.2.new)]
                    maha2.2     <- maha2.2.new # need to store maha-distances for interpolation
                    length.maha <- n # store length of 'maha2.2' and 'weights'
                    if(verbose >= 3) cat(".") # print dot after estimation of weights.
                } else {
                    ## Linearly interpolate 'weights' to get new weights
                    weights.new          <- rep(NA, n)
                    curr.index           <- 1 # index to look for close values in 'maha2.2'
                    notInterpol          <- rep(NA, n)
                    notInterpolcounter   <- 1
                    for(ind in 1:n) {
                        curr.maha2.2 <- maha2.2.new[ind]
                        if(curr.maha2.2 < maha2.2[1] || curr.maha2.2 > maha2.2[length.maha]) {
                            notInterpol[notInterpolcounter] <- ind
                            notInterpolcounter <- notInterpolcounter + 1
                        } else {
                            ## Start looking for close maha values in 'maha2.2'
                            found <- FALSE
                            while(!found && curr.index < length.maha) {
                                if(maha2.2[curr.index] <= curr.maha2.2 &&
                                   curr.maha2.2 <= maha2.2[curr.index+1]) {
                                    found <- TRUE
                                    ## Now check if we can interpolate (ie rel.error small)
                                    if(abs(weights[curr.index+1] - weights[curr.index])/
                                       weights[curr.index+1] < control$weights.interpol.reltol) {
                                        weights.new[ind] <- weights[curr.index] +
                                            (curr.maha2.2 - maha2.2[curr.index])*
                                            (weights[curr.index+1] - weights[curr.index])/
                                            (maha2.2[curr.index+1]-maha2.2[curr.index])
                                    } else {
                                        ## if not, will use 'get.weights' for this maha.
                                        notInterpol[notInterpolcounter] <- ind
                                        notInterpolcounter <- notInterpolcounter + 1
                                    }
                                } else {
                                    curr.index <- curr.index + 1
                                }
                            }
                        }
                    }
                    ## Now need to approximate weights for those maha in 'notInterpol'
                    if(notInterpolcounter > 1) {
                        notInterpol <- notInterpol[1:(notInterpolcounter-1)]
                        weights.new[notInterpol] <- weights_(maha2.2.new[notInterpol],
                                                                     qW = qW, nu = nu.est,
                                                                     lrdet = lrdet, d = d,
                                                                     special.mix = special.mix,
                                                                     control = control,
                                                                     verbose = verbose)$weights
                        ## Add estimated weights to 'weights' and corresponding
                        ## 'maha2.2.new' to 'maha2.2' so that they can be reused
                        maha2.2 <- c(maha2.2, maha2.2.new[notInterpol])
                        temp.ordering <- order(maha2.2) # only needed here
                        maha2.2 <- maha2.2[temp.ordering] # sort 'maha2.2' again
                        weights <- c(weights, weights.new[notInterpol])[temp.ordering] # and weights accordingly
                        length.maha <- length.maha + notInterpolcounter - 1
                    }
                    ## Recover original ordering and set negative weights to zero.
                    weights.new <- weights.new[order(order.maha2.2.new)] ## TODO: Omit?
                    if(verbose >= 3) cat(".") # print dot after estimation of weights.
                } # done estimating 'weights.new'

                ## Get new 'scale.est': 1/n * sum_{i=1}^n weights_i (x_i-mu)(x_i-mu)^T
                ## where 'mu' corresponds to current 'loc.est'
                scale.est.new <- crossprod(sqrt(weights.new)*sweep(x, 2, loc.est,
                                                                   check.margin = FALSE))/n
                ## Get new 'loc.est': sum_{i=1}^n weights_i x_i / (sum weights)
                ## as.vector because we need 'loc.est' as a vector, not (d, 1) matrix
                loc.est.new <- as.vector(crossprod(x, weights.new)/sum(weights.new))
                ## Did we converge?
                scale.est.rel.diff <- abs((scale.est - scale.est.new)/scale.est)
                loc.est.rel.diff   <- abs((loc.est - loc.est.new)/loc.est)
                converged.locscale <- (max(loc.est.rel.diff) < control$ECME.rel.conv.tol[2]) &&
                    (max(scale.est.rel.diff) < control$ECME.rel.conv.tol[2])
                ## Update counter
                iter.locscaleupdate <- iter.locscaleupdate + 1
                ## Update 'loc.est' and 'scale.est'
                loc.est     <- loc.est.new
                scale.est   <- scale.est.new
            } # done updating 'loc.est' and 'scale.est'
            if(verbose >= 3) cat(paste0(".DONE (", iter.locscaleupdate, " iterations needed)", '\n'))

            ## Update 'nu.est', if desired/necessary:
            if(control$ECMEstep.do.nu) {
                ## New subsample used for this 'nu' update?
                if(control$resample && size.subsample < n) {
                    if(exists(".Random.seed")) rm(".Random.seed") # destroy the reseted seed
                    runif(1) # get a new seed
                    sampled.ind <- sample(n, size.subsample)
                    tx.sub      <- tx[,sampled.ind, drop = FALSE]
                    seed        <- .Random.seed
                }
                ## Optimize neg.log.likelihood over 'nu'
                if(verbose >= 3) cat(paste0("    Optimizing likelihood over 'nu' with new 'loc' and 'scale'"))
                est.obj <- estim_nu(tx, qW = qW, init.nu = nu.est,
                                    loc = loc.est, scale = scale.est,
                                    mix.param.bounds = mix.param.bounds,
                                    special.mix = special.mix, control = control,
                                    control.optim = control$control.optim,
                                    verbose = verbose)
                nu.est.new        <- est.obj$nu.est
                nu.est.rel.diff   <- abs((nu.est.new - nu.est)/nu.est)
                nu.est            <- nu.est.new
                max.ll            <- est.obj$max.ll
                if(verbose >= 3) cat(paste0("DONE (", est.obj$ll.counts, " calls to likelihood needed)", '\n'))
            } else {
                nu.est.rel.diff <- 0
            }
            ## Did we converge?
            converged <- if(iter.ECME >= control$ECME.miniter) {
                             prod(abs(nu.est.rel.diff) < control$ECME.rel.conv.tol[3])
                         } else FALSE
            ## Update counter and 'nu.Ests'
            iter.ECME <- iter.ECME + 1
            if(control$addReturns) {
                ## Store new 'nu.est' along with log-likelihood
                nu.Ests[current.iter.total, ] <- c(nu.est, -max.ll)
                ## If 'converged', set all future iteration values to current one
                if(converged) {
                    nu.Ests[current.iter.total:(dim(nu.Ests)[1]), ] <-
                        matrix(c(nu.est, -max.ll), ncol = mix.param.length+1,
                               nrow = dim(nu.Ests)[1]+1-current.iter.total,
                               byrow = TRUE)
                    iter.converged <- current.iter.total
                }
                current.iter.total <- current.iter.total + 1
            }
        } # end while()

        iter.ECME < control$ECME.maxiter & !converged

        if(iter.ECME == control$ECME.maxiter & !converged & verbose >= 1)
            warning("Maximum number of ECME iterations exhausted, consider increasing 'ECME.maxiter' in the 'control' argument.")
    } #end if(control$ECMEstep)

    ## 3 Another last 'nu.est' with *full* sample? #############################

    if(control$laststep.do.nu) {
        if(verbose >= 2)
            cat(paste0("Step 3: One last 'nu' update", '\n'))
        if(verbose >= 3)
            cat(paste0("  Optimizing likelihood over 'nu' with new 'loc' and 'scale'"))
        ## One last nu update with the *full* sample and 'control.optim.laststep' as control
        ## (as oppsed to 'control.optim')
        est.obj <- estim_nu(tx, qW = qW, init.nu = nu.est, loc = loc.est,
                            scale = scale.est, mix.param.bounds = mix.param.bounds,
                            special.mix = special.mix, control = control,
                            control.optim = control$control.optim.laststep,
                            verbose = verbose)
        nu.est <- est.obj$nu.est
        max.ll <- est.obj$max.ll
        if(control$addReturns) {
            ## Store new 'nu.est' along with log-likelihood
            nu.Ests[dim(nu.Ests)[1], ] <- c(nu.est, -max.ll)
            current.iter.total <- current.iter.total + 1
        }
        if(verbose >= 3) cat(paste0("DONE (", est.obj$ll.counts, " calls to likelihood needed)", '\n'))
    }
    if(verbose >= 2) cat(paste0("RETURN.", '\n'))

    ## Return
    if(control$addReturns) {
        list(nu = nu.est, loc = loc.est, scale = scale.est, iter = iter.ECME,
             max.ll = -max.ll, nu.Ests = nu.Ests, iter.converged = iter.converged)
    } else {
        list(nu = nu.est, loc = loc.est, scale = scale.est, iter = iter.ECME,
                    max.ll = -max.ll)
    }
}
