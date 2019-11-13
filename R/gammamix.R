### pgammamix() ###################################################################

##' @title Distribution function of the mahalanobis distance of a normal variance mixture
##' @param x n-vector of evaluation points
##' @param qmix see ?pnvmix()
##' @param lower.tail logical if P(X <= m) or P(X>m) shd be returned       
##' @param d dimension of the underlying normal variance mixture
##' @param control see ?pnvmix()
##' @param verbose see ?pnvmix()
##' @param ... asee ?pnvmix()
##' @return n-vector with computed cdf values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert

pgammamix <- function(x, qmix, d, lower.tail = TRUE, 
                      control = list(), verbose = TRUE, ...)
{
   ## Checks
   if(!is.vector(x)) x <- as.vector(x) 
   n <- length(x) # length of input
   ## Deal with algorithm parameters, see also get.set.parameters():
   ## get.set.parameters() also does argument checking, so not needed here.
   control <- get.set.parameters(control)
   
   ## 1 Define the quantile function of the mixing variable ###################
   special.mix <- NA 
   qW <- if(is.character(qmix)) { # 'qmix' is a character vector
      qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma", "pareto"))
      switch(qmix,
             "constant" = {
                special.mix <- "constant"
                function(u) rep(1, length(u))
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
                   special.mix <- "inverse.gamma"
                   df2 <- df / 2
                   function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
                } else {
                   special.mix <- "constant"
                   function(u) rep(1, length(u))
                }
             },
             "pareto"= {
                if(hasArg(alpha)){
                   alpha <- list(...)$alpha
                } else {
                   stop("'qmix = \"pareto\"' requires 'alpha' to be provided.")
                }
                special.mix <- "pareto"
                function(u) (1-u)^(-1/alpha)
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
   
   ## Build result object 
   pres <- rep(0, n) # n-vector of results
   notNA <- which(!is.na(x)) 
   pres[!notNA] <- NA
   x <- x[notNA, drop = FALSE] # non-missing data (rows)
   ## Counter
   numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
   
   ## Deal with special distributions
   if(!is.na(special.mix)){
      if(!(special.mix == "pareto")){
         ## Only for "inverse.gamma" and "constant" do we have analytical forms:
         pres[notNA] <- switch(special.mix,
                               "inverse.gamma" = {
                                  ## D^2 ~ d* F(d, nu)
                                  pf(x/d, df1 = d, df2 = df, lower.tail = lower.tail)
                               },
                               "constant" = {
                                  ## D^2 ~ chi^2_d = Gamma(shape = d/2, scale = 2)
                                  pgamma(x, shape = d/2, scale = 2, lower.tail = lower.tail)
                               })
         ## Return in those cases
         attr(pres, "error")   <- rep(0, n)
         attr(pres, "numiter") <- numiter
         return(pres)
      }
   }
   
   ## Define various quantites for the RQMC procedure
   dblng           <- (control$increment == "doubling")
   B               <- control$B # number of randomizations
   current.n       <- control$fun.eval[1] #initial sample size
   total.fun.evals <- 0
   ZERO            <- .Machine$double.neg.eps
   ## Absolte/relative precision?
   if(is.na(control$pgammamix.reltol)){
      ## Use absolute error
      tol <- control$pgammamix.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error
      tol <- control$pgammamix.reltol
      do.reltol <- TRUE
   }
   ## Store seed if 'sobol' is used to get the same shifts later:
   if(control$method == "sobol") {
      if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
      seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
   }
   ## Additional variables needed if the increment chosen is "dblng"
   if(dblng) {
      if(control$method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Matrix to store RQMC estimates
   rqmc.estimates <- matrix(0, ncol = n, nrow = B)
   ## Will be needed a lot:
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ## Initialize 'max.error' to > tol so that we can enter the while loop:
   max.error <- tol + 42
   
   ## 2 Main loop ###############################################################
   ## while() runs until precision abstol is reached or the number of function
   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
   ## the desired probability are calculated.
   while(max.error > tol && numiter < control$max.iter.rqmc &&
         total.fun.evals < control$fun.eval[2])
   {
      ## Reset seed to have the same shifts in sobol( ... )
      if(control$method == "sobol" && numiter > 0)
         .Random.seed <<- seed # reset seed to have the same shifts in sobol( ... )
      for(b in 1:B){
         
         ## 2.1 Get the point set ###########################################
         U <- switch(control$method,
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
                     }) # sorted for later!
         
         ## 2.2 Evaluate the integrand at the (next) point set #############
         W <- pmax(qW(U), ZERO) # realizations of the mixing variable
         next.estimate <- .colMeans(pgamma( 
            sapply(x, function(i) i/W), shape = d/2, scale = 2, lower.tail = lower.tail), 
            current.n, n, 0)
         
         ## 2.3 Update RQMC estimates #######################################
         rqmc.estimates[b,] <-
            if(dblng) {
               ## In this case both, rqmc.estimates[b,] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b,] + next.estimate)/denom
            } else {
               ## In this case, rqmc.estimates[b,] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b,] + next.estimate)/(numiter + 1)
            }
      } # end for(b in 1:B)
      
      ## Update of various variables
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(dblng) {
         ## Change denom and useksip (exactly once, in the first iteration)
         if(numiter == 0){
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Total number of function evaluations:
      total.fun.evals <- total.fun.evals + B * current.n
      numiter <- numiter + 1
      ## Update error. The following is slightly faster than 'apply(..., 2, var)' 
      pres[notNA] <- .colMeans(rqmc.estimates, B, n , 0)
      vars <- .colMeans((rqmc.estimates - rep(pres, each = B))^2, B, n, 0)
      errors <- if(!do.reltol){
         sqrt(vars)*CI.factor.sqrt.B
      } else {
         sqrt(vars)/pres*CI.factor.sqrt.B
      }
      max.error <- max(errors)
   } # while()
   if(verbose && max.error > tol)
      warning("Tolerance not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
   
   ## 3 Return ##################################################################
   attr(pres, "error")   <- errors 
   attr(pres, "numiter") <- numiter
   pres
}


### qgammamix() ###################################################################

##' @title Quantile function of the mahalanobis distance of a normal variance mixture
##' @param u see ?qnvmix()
##' @param qmix see ?qnvmix()
##' @param d dimension of the underlying normal variance mixture
##' @param control see ?qnvmix()
##' @param verbose see ?qnvmix()
##' @param q.only see ?qnvmix()
##' @param stored.values see ?qnvmix()
##' @param ... see ?qnvmix()
##' @return see ?qnvmix()  
qgammamix <- function(u, qmix, d, control = list(), verbose = TRUE, q.only = TRUE,
                      stored.values = NULL, ...)
{
   quantile.internal(u, qmix = qmix, which = "maha2", d = d, control = control,
                     verbose = verbose, q.only = q.only, 
                     stored.values = stored.values, ...)
}


### rgammamix() ###################################################################

##' @title Random Number Generator for Gamma mixtures ( = Mahalanobis distances for 
##'        normal mixture models)
##' @param n see ?rnvmix()
##' @param rmix ?rnvmix()
##' @param qmix ?rnvmix()
##' @param d dimension of the underlying normal variance mixture, see below under 'return'
##' @param method ?rnvmix()
##' @param skip ?rnvmix()
##' @param ... n-vector or realizations
##' @return n-vector of samples of W*X where W~qmix/rmix and X~chi^2_d = Gamma(shape = d/2, scale = 2)
##' @author Marius Hofert and Erik Hintz
##' @note - For the Student t distribution, W ~ df/rchisq(n, df = df) but
##'         rchisq() simply calls rgamma(); see ./src/nmath/rchisq.c
##'         => W ~ 1/rgamma(n, shape = df/2, rate = df/2)
##'       - For a generalized inverse Gaussian distribution one could use:
##'         + "Runuran": faster if n large and parameters fixed; based on density
##'         + "GIGrvg":  faster if n small and often called with several parameters
##'         see examples of 'GIGrvg' for both methods
##'       - user friendly wrappers are provided in 'rnvmix()' and 'rgammamix()'    
rgammamix <- function(n, rmix = NULL, qmix = NULL, d, method = c("PRNG", "sobol", "ghalton"),
                      skip = 0, ...)
{
   rmix.internal(n, rmix = rmix, qmix = qmix, loc = rep(0, d), scale = diag(d),
                 factor = diag(d), method = method, skip = skip, which = "maha2",
                 ...)
   
}


### dgammamix() ###################################################################

##' @title Density function of the mahalanobis distance of a normal variance mixture
##' @param x n-vector of evaluation points
##' @param qmix see ?pnvmix()
##' @param d dimension of the underlying normal variance mixture
##' @param control see ?pnvmix()
##' @param verbose see ?pnvmix()
##' @param log logical if log-density shall be returned
##' @param ... see ?pnvmix()
##' @return n-vector with computed density values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert

dgammamix <- function(x, qmix, d, control = list(), verbose = TRUE, log = FALSE, 
                      ...){
   ## Checks
   if(!is.vector(x)) x <- as.vector(x)
   n <- length(x) # dimension
   stopifnot(d >= 1)
   verbose <- as.logical(verbose) 
   ## Deal with algorithm parameters, see also get.set.parameters():
   ## get.set.parameters() also does argument checking, so not needed here.
   control <- get.set.parameters(control)
   
   ## Build result object (log-density)
   lres  <- rep(-Inf, n) # n-vector of results
   notNA <- which(!is.na(x))
   lres[!notNA] <- NA
   x <- x[notNA, drop = FALSE] # non-missing data (rows)
   
   ## 1 Define the quantile function of the mixing variable ###################
   ## If 'mix' is "constant", "inverse.gamma" or "pareto", we use the analytical formulas
   special.mix <- NA 
   qW <- if(is.character(qmix)) { # 'qmix' is a character vector
      qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma", "pareto"))
      switch(qmix,
             "constant" = {
                special.mix <- "constant"
                function(u) rep(1, length(u))
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
                   special.mix <- "inverse.gamma"
                   df2 <- df / 2
                   mean.sqrt.mix <- sqrt(df) * gamma(df2) / (sqrt(2) * gamma((df+1)/2)) # used for preconditioning
                   function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
                } else {
                   special.mix <- "constant"
                   mean.sqrt.mix <- 1 # used for preconditioning
                   function(u) rep(1, length(u))
                }
             },
             "pareto"= {
                if(hasArg(alpha)){
                   alpha <- list(...)$alpha
                } else {
                   stop("'qmix = \"pareto\"' requires 'alpha' to be provided.")
                }
                special.mix <- "pareto"
                mean.sqrt.mix <- if(alpha > 0.5) alpha/(alpha-0.5) else NULL
                function(u) (1-u)^(-1/alpha)
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
   
   ## Counter
   numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
   ## Deal with the different distributions
   
   ## Deal with special distributions
   if(!is.na(special.mix)){
      if(!(special.mix == "pareto")){
         ## Only for "inverse.gamma" and "constant" do we have analytical forms:
         lres[notNA] <- switch(special.mix,
                               "inverse.gamma" = {
                                  ## D^2 ~ d* F(d, nu)
                                  df(x/d, df1 = d, df2 = df, log = TRUE) - log(d)
                               },
                               "constant" = {
                                  ## D^2 ~ chi^2_d = Gamma(shape = d/2, scale = 2)
                                  dgamma(x, shape = d/2, scale = 2, log = TRUE)
                               })
         ## Return in those cases
         attr(lres, "error")   <- rep(0, n)
         attr(lres, "numiter") <- numiter
         if(log) return(lres) else return(exp(lres))
      }
   }
   ## General case of a multivariate normal variance mixture (RQMC)
   ## Prepare inputs for densmix.internal.RQMC
   ## Sort input 'x' increasingly and store ordering for later
   ordering.x <- order(x)
   x.ordered  <- x[ordering.x]
   ## Define log-constant for the integration
   lconst <- -lgamma(d/2) - d/2 * log(2) + (d/2-1)*log(x.ordered)
   ## Call internal dnvmix (which itself calls C-Code and handles warnings):
   ests <- densmix.internal(qW, maha2.2 = x.ordered/2, lconst = lconst,
                            d = d, control = control, verbose = verbose)
   ## Grab results, correct 'error' and 'lres' if 'log = FALSE'
   lres[notNA] <- ests$ldensities[order(ordering.x)]
   error <- if(log){
      ests$error[order(ordering.x)]
   } else {
      lres <- exp(lres)
      ests$error[order(ordering.x)]*pmax(lres[notNA], 1)
   }
   numiter <- ests$numiter
   
   ## Return
   ## Note that 'lres' was exponentiated already if necessary.
   attr(lres, "error")   <- error # these are absolute errors, no matter what!
   attr(lres, "numiter") <- numiter
   lres
}