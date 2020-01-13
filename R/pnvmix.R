### pnvmix() ###################################################################

##' @title Evaluate Integrand of pnvmix()
##' @param U (n, d) matrix of uniforms (evaluation points)
##' @param qmix see ?pnvmix
##' @param lower see ?pnvmix
##' @param upper see ?pnvmix
##' @param scale see ?pnvmix. Has to be full rank here though!
##' @param precond see ?get_set_param()
##' @param mean.sqrt.mix see ?get_set_param()
##' @param return.all logical if all function evaluations should be returned.
##' @param ... additional parameters passed to qmix()
##' @return if return.all = TRUE a n-vector with g(U) values, otherwise a 2-vector
##' with estimated mean(g(U)) and estimated var(g(U))
##' @note if return.all = TRUE, all computations are done in R, so can be slow.
##'       if return.all = FALSE, this function calls underlying C code.
##' @author Erik Hintz
##' @note - This function is *only* needed for numerical experiments.
##'       - It corresponds to 'g' in the paper.
pnvmix_g <- function(U, qmix, upper, lower = rep(-Inf, d), scale, precond,
                     mean.sqrt.mix = NULL, return.all = FALSE, verbose = TRUE, ...)
{
   ## Define the quantile function of the mixing variable
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
                } else if(hasArg(nu)) {
                   nu <- list(...)$nu
                   df <- nu
                } else {
                   stop("'qmix = \"inverse.gamma\"' requires 'df' to be provided.")
                }
                ## Still allow df = Inf (normal distribution)
                stopifnot(is.numeric(df), length(df) == 1, df > 0)
                if(is.finite(df)) {
                   special.mix <- "inverse.gamma"
                   df2 <- df / 2
                   mean.sqrt.mix <- sqrt(df) * gamma(df2) /
                      (sqrt(2) * gamma((df+1)/2)) # used for preconditioning
                   function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
                } else {
                   special.mix <- "constant"
                   mean.sqrt.mix <- 1 # used for preconditioning
                   function(u) rep(1, length(u))
                }
             },
             "pareto"= {
                if(hasArg(alpha)) {
                   alpha <- list(...)$alpha
                } else if(hasArg(nu)) {
                   nu <- list(...)$nu
                   alpha <- nu
                } else { 
                   stop("'qmix = \"pareto\"' requires 'alpha' to be provided.")
                }
                special.mix   <- "pareto"
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
   ## Dimension of the problem and number of evaluations
   d <- dim(scale)[1]
   n <- dim(U)[1]
   ## Factor (lower triangular)
   C        <- t(chol(scale))
   rank     <- d # only consider full rank case here
   k.factor <- rep(1, d)
   
   ## Precondition?
   if(precond && d > 2) {
      if(is.null(mean.sqrt.mix))
         mean.sqrt.mix <- mean(sqrt(qW(qrng::sobol(n = 2^12, d = 1, randomize = TRUE))))
      ## Check if provided/approximated' mean.sqrt.mix' is strictly positive
      if(mean.sqrt.mix <= 0)
         stop("'mean.sqrt.mix' has to be positive (possibly after being generated in pnvmix())")
      temp <- precondition(lower, upper = upper, scale = scale, factor = C,
                           mean.sqrt.mix = mean.sqrt.mix)
      if(is.null(temp)) {
         ## Preconditioning did not work, continue with original inputs
         if(verbose) warning("Preconditioning led to (numerically) singular 'scale',
                             continuing with original input.")
      } else {
         a <- temp$lower
         b <- temp$upper
         R <- temp$scale
         C <- temp$factor
      }
   } else {
      a <- lower
      b <- upper
      R <- scale # C did not change
   }
   ## For evaluating qnorm() close to 0 and 1
   ONE <- 1-.Machine$double.neg.eps
   ZERO <- .Machine$double.eps
   ## Transform inputs to realizations of the mixing variable
   ## U will be a matrix with *d+1* columns: Column 1: Realizations of sqrt(mix),
   ## Columns 2 to d: uniforms, Column d+1: Antithetic realization of Column 1
   if(!is.matrix(U)) U <- as.matrix(U)
   U <- cbind(sqrt(qW(U[, 1])), U[, 2:d], sqrt(qW(1 - U[, 1])))
   
   if(return.all) {
      ## Matrix to store results (y_i from paper)
      Yorg <- matrix(NA, ncol = d - 1, nrow = dim(U)[1])
      Yant <- matrix(NA, ncol = d - 1, nrow = dim(U)[1])
      ## First 'iteration' (d1, e1 in the paper)
      dorg <- pnorm(a[1] / (U[, 1] * C[1, 1]))
      dant <- pnorm(a[1] / (U[, d+1] * C[1, 1]))
      eorg <- pnorm(b[1] / (U[, 1] * C[1, 1]))
      eant <- pnorm(b[1] / (U[, d+1] * C[1, 1]))
      forg <- eorg - dorg
      fant <- eant - dant
      ## Recursively calculate (e_i - d_i)
      for(i in 2:d) {
         ## Store realization:
         Yorg[,(i-1)] <- qnorm( pmax( pmin(dorg + U[, i]*(eorg-dorg), ONE), ZERO))
         Yant[,(i-1)] <- qnorm( pmax( pmin(dant + (1-U[, i])*(eant-dant), ONE), ZERO))
         ## Update d__, e__, f___:
         dorg <- pnorm((a[i]/U[,1]   - Yorg[,1: (i-1)] %*% as.matrix(C[i, 1:(i-1)]))/C[i,i])
         dant <- pnorm((a[i]/U[,d+1] - Yant[,1: (i-1)] %*% as.matrix(C[i, 1:(i-1)]))/C[i,i])
         eorg <- pnorm((b[i]/U[,1]   - Yorg[,1: (i-1)] %*% as.matrix(C[i, 1:(i-1)]))/C[i,i])
         eant <- pnorm((b[i]/U[,d+1] - Yant[,1: (i-1)] %*% as.matrix(C[i, 1:(i-1)]))/C[i,i])
         forg <- forg*(eorg-dorg)
         fant <- fant*(eant-dant)
      }
      ## Return:
      (forg+fant)/2
   } else {
      .Call("eval_nvmix_integral",
            lower    = as.double(a),
            upper    = as.double(b),
            U        = as.double(U),
            n        = as.integer(n),
            d        = as.integer(d),
            r        = as.integer(rank),
            kfactor  = as.integer(k.factor),
            factor   = as.double(C),
            ZERO     = as.double(ZERO),
            ONE      = as.double(ONE))
   }
}


##' @title Cholesky Factor for Positive-Semidefinite Matrices
##' @param mat (n,n) symmetric, positive-semidefinite matrix.
##' @return list of length 2:
##'         - C: (n,n) lower triangular matrix such that C % * % t(C) = mat
##'         - D: n-vector with diagonal elements of 'C'
##' @author Erik Hintz and Marius Hofert
##' @note Internal function being called by pnvmix().
cholesky_ <- function(mat, tol = 1e-12)
{
   n <- dim(mat)[1] # dimension
   stopifnot(dim(mat)[2] == n)
   ## First try 'chol()' (=>fast)
   C <- tryCatch(t(chol(mat)), error = function(e) e)
   if(is.matrix(C) && all.equal(dim(C), rep(n, 2))) {
      ## C is the desired Cholesky factor
      ## Grab diagonal
      diag.elements <- diag(C)
   } else {
      ## In this case, 't(chol(scale))' returned an error so that we manually
      ## compute the Cholesky factor of the *singular* matrix 'mat'
      C <- matrix(0, ncol = n, nrow = n) # initialize Cholesky factor
      diag.elements <- rep(NA, n)
      for(col in 1:n) {
         dsq <- mat[col, col] - sum(C[col, 1:(col-1)]^2) # C(col,col)^2
         if(dsq < 0) stop("Matrix not positive semi-definite")
         d <- if(dsq < abs(mat[col, col] * tol)) 0 else sqrt(dsq) # set 0 if too small
         C[col, col] <- d
         diag.elements[col] <- d # store diagnonal element
         if(col < n && d > 0) { # calculate the remaining elements in column 'col'
            for(row in (col+1):n) {
               C[row, col] <- (mat[row, col] -
                                  sum(C[row, 1:(row-1)]*C[col, 1:(row-1)]))/d
            }
         }
      }
   }
   list(C = C, D = diag.elements)
}


##' @title Swap Variables i and j in a, b and R
##' @param i variable to be switched with j
##' @param j variable to be switched with i
##' @param lower d-vector of lower evaluation limits
##' @param upper d-vector of upper evaluation limits
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @return list with lower, upper and scale after components/rows/columns
##'         i and j have been switched
##' @author Erik Hintz and Marius Hofert
swap <- function(i, j, lower, upper, scale)
{
   ## Build vector
   ij <- c(i,j)
   ji <- c(j,i)
   ## Reorder lower, upper
   lower[ij] <- lower[ji]
   upper[ij] <- upper[ji]
   ## Reorder scale
   wo.ij <- setdiff(seq_len(nrow(scale)), ij)
   temp_i <- scale[wo.ij,i, drop = FALSE]
   temp_j <- scale[wo.ij,j, drop = FALSE]
   temp_ii <- scale[i,i]
   scale[wo.ij,i] <- temp_j
   scale[wo.ij,j] <- temp_i
   scale[i,wo.ij] <- temp_j
   scale[j,wo.ij] <- temp_i
   scale[i,i] <- scale[j,j]
   scale[j,j] <- temp_ii
   ## Return
   list(lower = lower, upper = upper, scale = scale)
}


##' @title Reorder Limits and Scale Matrix According to a Permutation
##' @param perm d-vector giving the desired permutation
##' @param upper see ?pnvmix
##' @param lower see ?pnvmix
##' @param scale see ?pnvmix
##' @return list with lower, upper and scale after components/rows/columns
##'         have been switched according to perm
##' @author Erik Hintz
##' @note No input checking is done. This function is mostly for experimenting.
reorder_limits_scale <- function(perm, upper, lower = rep(-Inf, d), scale = diag(d))
{
   d <- length(upper)
   ## Vector to save current positions of original variables
   org.perm <- 1:d
   for(i in 1:d) {
      ## Variable at position i before swapping
      curr.i <- org.perm[i]
      ## Variable at position i after swapping
      next.i <- perm[i]
      ## Do we need to swap?
      if(curr.i != next.i) {
         ## Index of the variable to be put to position i
         ind.next.i <- which(org.perm == next.i)
         ## Swap positions i and ind.next.i
         tmp <- swap(i, ind.next.i, lower = lower, upper = upper,
                     scale = scale)
         lower <- tmp$lower
         upper <- tmp$upper
         scale <- tmp$scale
         ## Update position vector
         org.perm[i] <- next.i
         org.perm[ind.next.i] <- curr.i
      }
   }
   list(lower = lower, upper = upper, scale = scale)
}


##' @title Variance of a Normal Rv over a Truncated Interval (a, b)
##' @param a l-vector
##' @param b l-vector
##' @return l-vector
##' @author Erik Hintz
##' @note formula in Genz and Bretz (2009, p. 38)
trunc_var <- function(a, b)
{
   p.diff <- pnorm(b) - pnorm(a)
   ## Cases -Inf * 0 and Inf * 0:
   adnorma <- a*dnorm(a)
   adnorma[which(is.nan(adnorma))] <- 0
   bdnormb <- b*dnorm(b)
   bdnormb[which(is.nan(bdnormb))] <- 0
   ## Return
   1 + (adnorma - bdnormb)/p.diff - ((dnorm(a)-dnorm(b))/p.diff)^2
}


##' @title Preconditioning (Reordering Variables According to their Expected
##'        Integration Limits)
##' @param lower see ?pnvmix
##' @param upper see ?pnvmix
##' @param scale (d,d) positive definite 'scale' matrix.
##' @param factor Cholesky factor (lower triangular matrix) of 'scale'
##' @param mean.sqrt.mix E(sqrt(W)) where W is the rv corresponding to 'qmix'
##' @param precond.method character, "ExpLength" (=> sorted by expected length),
##'        otherwise sorted by 'trunc_var()'.
##' @param tol if a calculated diagonal element of factor is < sqrt(tol),
##'        factor is deemed singular and preconditioning fails.
##' @return list of length 4 with reordered integration limits, scale matrix and
##'         Cholesky factor as well as a d-vector 'perm' giving the ordering obtained.
##'         If preconditioning was unsuccessful, 'NULL' is returned
##' @author Erik Hintz and Marius Hofert
##' @note See Genz and Bretz (2002, p. 957).
##' It may happen that the original 'scale' admits a (numerically stable full rank)
##' 'factor' and that the reordered 'scale' is numerically singular so that a
##' full rank 'factor' cannot be found. This is detected by 'precondition'
##' and if this happens, NULL is returned. pnvmix1() then throws a warning
##' and estimation is carried out with un-preconditioned inputs.
precondition <- function(lower, upper, scale, factor, mean.sqrt.mix,
                         precond.method = "ExpLength", tol = 1e-16)
{
   d <- length(lower)
   y <- rep(0, d - 1)
   perm <- 1:d
   exp.lengths <- rep(NA, d)
   ## Main
   for(j in 1:(d-1)) {
      ## Case j = 1 somewhat special
      if(j == 1) {
         denom <- sqrt(diag(scale))
         c <- 0
      } else {
         denom <- diag(scale)[j:d] - .rowSums(factor[j:d, 1:(j-1), drop = FALSE]^2,
                                              m = d-j+1, n = j-1)
         if(any(denom <= 0)) return(NULL) else denom <- sqrt(denom)
         c <- factor[j:d, 1:(j-1), drop = FALSE] %*% y[1:(j-1)]
      }
      
      ## Transformed limits with indices greater than j:
      next.uppers <- (upper[j:d] / mean.sqrt.mix - c) / denom
      next.lowers <- (lower[j:d] / mean.sqrt.mix - c) / denom
      
      ## Find i = argmin { <expected length of interval j> }
      i <- if(precond.method == "ExpLength") {
         which.min(pnorm(next.uppers) - pnorm(next.lowers)) + j - 1
      } else {
         ## Find i = argmin { <truncated variance of variable j> }
         which.min(trunc_var(next.lowers, next.uppers)) + j - 1
      }
      
      ## Swap i and j if they are different
      if(i != j) {
         tmp   <- swap(i = i, j = j, lower = lower, upper = upper, scale = scale)
         lower <- tmp$lower
         upper <- tmp$upper
         scale <- tmp$scale
         perm[c(i, j)] <- perm[c(j, i)]
         ## If j>1 and an actual swap has occured, need to reorder Cholesky factor:
         if(j > 1) {
            factor[c(i,j),]   <- factor[c(j,i),, drop = FALSE]
            factor[j,(j+1):i] <- matrix(0, ncol = i - j, nrow = 1)
         }
      }
      ## Update Cholesky factor
      if(j == 1) {
         factor[1, 1] <- sqrt(scale[1, 1])
         factor[2:d, 1] <- scale[2:d, 1, drop = FALSE] / factor[1, 1]
         ## Store y1
         y[1] <- -(dnorm(upper[1]/mean.sqrt.mix) - dnorm(lower[1]/mean.sqrt.mix)) /
            (max(pnorm(upper[1]/mean.sqrt.mix) - pnorm(lower[1]/mean.sqrt.mix), tol)) # avoid division by zero
      } else {
         factorjjsq <- scale[j,j] - sum(factor[j,1:(j-1)]^2)
         if(any(factorjjsq <= tol)) return(NULL) else factor[j,j] <- sqrt(factorjjsq)
         factor[(j+1):d, j] <-
            if(j < d-1) {
               (scale[(j+1):d, j] - factor[(j+1):d, 1:(j-1), drop = FALSE] %*%
                   factor[j, 1:(j-1)]) / factor[j, j]
            } else {
               (scale[(j+1):d, j] - factor[(j+1):d, 1:(j-1)] %*%
                   factor[j, 1:(j-1)]) / factor[j, j]
            }
         ## Get yj
         scprod     <- sum(factor[j, 1:(j-1)] * y[1:(j-1)]) # needed twice
         low.j.up.j <- c(lower[j] / mean.sqrt.mix - scprod,
                         upper[j] / mean.sqrt.mix - scprod) / factor[j, j]
         y[j] <- (dnorm(low.j.up.j[1]) - dnorm(low.j.up.j[2])) /
            (max(pnorm(low.j.up.j[2]) - pnorm(low.j.up.j[1]), tol))
      }
   } # for()
   factorddsq <- scale[d, d] - sum(factor[d, 1:(d-1)]^2)
   if(factorddsq > tol) {
      factor[d,d] <- sqrt(factorddsq)
   } else {
      ## Try 'chol()' for reordered 'scale' (more accurate)
      factor <- tryCatch(t(chol(scale)), error = function(e) e)
      if(!is.matrix(factor)) return(NULL) # else 'factor' is correct => return
   }
   ## Return
   list(lower = lower, upper = upper, scale = scale, factor = factor, perm = perm)
}


##' @title Distribution Function of a Multivariate Normal Variance Mixture
##'        for a Single Observation
##' @param upper d vector
##' @param lower d vector (<= upper)
##' @param qW quantile function of the mixture distribution; built inside pnvmix();
##'        note: different from (the more general) 'qmix'
##' @param is.const.mix logical, TRUE if qmix is constant (=> normal distutions)
##' @param mean.sqrt.mix see details in ?pnvmix
##' @param loc see details in ?pnvmix
##' @param scale see details in ?pnvmix
##' @param factor see details in ?pnvmix
##' @param k.factor vector of length rank(scale) with heights of each step in
##'        'factor'
##' @param method see details in ?pnvmix
##' @param precond see details in ?pnvmix
##' @param tol absolute/relative error tolerance depending on 'do.reltol'
##' @param do.reltol logical if relative precision is required
##' @param CI.factor see details in ?pnvmix
##' @param fun.eval see details in ?pnvmix
##' @param increment see details in ?pnvmix
##' @param B see details in ?pnvmix
##' @param verbose see ?pnvmix
##' @param ... see details in ?pnvmix
##' @return list of length 3:
##'         - value: computed probability
##'         - error: error estimate
##'         - numiter: number of iterations needed
##' @author Erik Hintz and Marius Hofert
##' @note Internal function being called by pnvmix.
pnvmix1 <- function(upper, lower = rep(-Inf, d),
                    qW = NULL, is.const.mix = FALSE, mean.sqrt.mix,
                    loc = rep(0, d), scale = diag(d), factor = NULL,
                    k.factor = rep(1, d),
                    method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
                    tol = 1e-3, do.reltol = FALSE, CI.factor = 3.3,
                    fun.eval = c(2^6, 1e8),
                    max.iter.rqmc = 15, increment = c("doubling", "num.init"),
                    B = 15, verbose = TRUE, ...)
{
   rank <- length(k.factor)
   d    <- sum(k.factor)
   stopifnot(length(lower) == d, lower <= upper, is.function(qW))
   if(any(lower == upper))
      return(list(value = 0, error = 0, numiter = 0))
   
   ## Get (lower triangular) Cholesky factor if not provided
   ## This is only needed if the internal function pnvmix1() is called directly,
   ## if pnvmix() is used, it was determined there.
   ## Note: This will only work if 'scale' is full rank. In the singular case,
   ## pnvmix() needs to be called and handles the singularity issues
   ##  (more efficient, as done only once)
   if(is.null(factor)) factor <- t(chol(scale)) # lower triangular Cholesky factor
   
   ## Preconditioning (resorting the limits; only for d > 2 and non-singular case)
   if(precond && d > 2 && rank == d) {
      ## Note that 'mean.sqrt.mix' has already been calculated in pnvmix()
      temp <- precondition(lower = lower, upper = upper, scale = scale,
                           factor = factor, mean.sqrt.mix = mean.sqrt.mix)
      if(is.null(temp)) {
         ## Preconditioning did not work, continue with original inputs
         if(verbose)
            warning("Preconditioning led to (numerically) singular 'scale',
                           continuing with original input.")
      } else {
         lower <- temp$lower
         upper <- temp$upper
         factor <- temp$factor
      }
   }
   
   ## 1 Basics for while loop below ###########################################
   
   ## Error is calculated as CI.factor * sd( estimates) / sqrt(B)
   ## For performance:
   CI.factor.sqrt.B <- CI.factor / sqrt(B)
   ## Grab the number of sample points for the first iteration
   current.n <- fun.eval[1] #
   ## Vector to store the B RQMC estimates
   rqmc.estimates <- rep(0, B)
   ## Initialize error to > 'tol' so that we can enter the while loop below
   error <- tol + 42
   ## Initialize the total number of function evaluations
   total.fun.evals <- 0
   ## Initialize counter for the number of iterations in the while loop below
   numiter <- 0
   
   ## It may happen that qnorm(u) for u too close to 1 (or 0) is evaluated;
   ## in those cases, u will be replaced by ONE and ZERO which is the largest
   ## (smallest) number different from 1 (0) such that qnorm(u) is not +/- Inf
   ZERO <- .Machine$double.xmin
   ONE <- 1-.Machine$double.neg.eps
   
   ## If method == sobol, we want the same random shifts in each iteration below,
   ## this is accomplished by reseting to the "original" seed
   if(method == "sobol") {
      if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
      seed <- .Random.seed
   }
   
   ## Additional variables needed if the increment chosen is "doubling"
   if(increment == "doubling") {
      if(method == "sobol") useskip <- 0
      denom <- 1
   }
   
   ## 2 Major while() loop ####################################################
   
   ## while() runs until precision 'tol' is reached or the number of function
   ## evaluations exceed fun.eval[2] or until 'max.iter.rqmc' is exhausted.
   ## In each iteration, B RQMC estimates of the desired probability are calculated.
   while(error > tol && total.fun.evals < fun.eval[2] && numiter < max.iter.rqmc)
   {
      if(method == "sobol" && numiter > 0)
         .Random.seed <<- seed # reset seed to have the same shifts in sobol(...)
      
      ## Get B RQCM estimates
      for(b in 1:B)
      {
         ## 2.1 Get the point set ###########################################
         
         ## If is.const.mix = TRUE, we only need (rank - 1) (quasi) random numbers
         ## (is.const.mix = TRUE and rank = 1 has already been dealt with)
         U <- if(is.const.mix) {
            U <- switch(method, # same 'U' to possibly avoid copying
                        "sobol" = {
                           if(increment == "doubling") {
                              qrng::sobol(n = current.n, d = rank - 1,
                                          randomize = TRUE,
                                          skip = (useskip * current.n))
                           } else {
                              qrng::sobol(n = current.n, d = rank - 1,
                                          randomize = TRUE,
                                          skip = (numiter * current.n))
                           }
                        },
                        "ghalton" = {
                           qrng::ghalton(n = current.n, d = rank - 1,
                                         method = "generalized")
                        },
                        "PRNG" = {
                           matrix(runif( current.n * (rank - 1)), ncol = rank - 1)
                        })
            ## First/last column contain 1s corresponding to "simulated"
            ## values from sqrt(mix)
            cbind(rep(1, current.n), U, rep(1, current.n))
         } else {
            U <- switch(method,
                        "sobol" = {
                           if(increment == "doubling") {
                              qrng::sobol(n = current.n, d = rank,
                                          randomize = TRUE,
                                          skip = (useskip * current.n))
                           } else {
                              qrng::sobol(n = current.n, d = rank,
                                          randomize = TRUE,
                                          skip = (numiter * current.n))
                           }
                        },
                        "ghalton" = {
                           qrng::ghalton(n = current.n, d = rank,
                                         method = "generalized")
                        },
                        "PRNG" = {
                           matrix(runif( current.n * rank), ncol = rank)
                        })
            
            ## Case d = 1 somewhat special again:
            if(d == 1) {
               cbind( sqrt(qW(U)), sqrt(qW(1 - U)) )
            } else {
               ## Column 1:sqrt(mix), Columns 2--r: unchanged (still uniforms),
               ## Column r + 1: antithetic realization of sqrt(mix)
               cbind(sqrt(qW(U[, 1])), U[, 2:rank], sqrt(qW(1 - U[, 1])))
            }
         }
         
         ## 2.2 Evaluate the integrand at the (next) point set ##############
         
         next.estimate <-
            if(d == 1) {
               ## Case of dimension 1: Don't need to approximate the
               ## multivariate normal df and can just use pnorm()
               ## Note that d = 1 for a pure normal or t df has already been addressed
               mean((pnorm(upper/U[,1])   - pnorm(lower/U[,1]) +
                        pnorm(upper/U[,d+1]) - pnorm(lower/U[,d+1])) / 2)
            } else {
               .Call("eval_nvmix_integral",
                     lower    = as.double(lower),
                     upper    = as.double(upper),
                     U        = as.double(U),
                     n        = as.integer(current.n),
                     d        = as.integer(d),
                     r        = as.integer(rank),
                     kfactor  = as.integer(k.factor),
                     factor   = as.double(factor),
                     ZERO     = as.double(ZERO),
                     ONE      = as.double(ONE))[1]
            }
         
         ## 2.3 Update RQMC estimates #######################################
         
         rqmc.estimates[b] <-
            if(increment == "doubling") {
               ## In this case both, rqmc.estimates[b] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b] + next.estimate) / denom
            } else {
               ## In this case, rqmc.estimates[b] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b] + next.estimate) / (numiter + 1)
            }
         
      } # end for(b in 1:B)
      ## Update of various variables
      ## Number of function evaluations
      ## (* 2 since antithetic variates are used in eval_nvmix_integral())
      total.fun.evals <- total.fun.evals + 2 * B * current.n
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(increment == "doubling") {
         ## Change denom and useksip (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Update error depending on 'do.reltol'
      error <- if(!do.reltol) {
         CI.factor.sqrt.B * sd(rqmc.estimates)
      } else {
         CI.factor.sqrt.B * sd(rqmc.estimates)/mean(rqmc.estimates)
      }
      numiter <- numiter + 1 # update counter
   } # while()
   ## Finalize
   value <- mean(rqmc.estimates) # calculate the RQMC estimator
   ## Return
   list(value = value, error = error, numiter = numiter)
}


##' @title Distribution Function of a Multivariate Normal Variance Mixture
##' @param upper (n,d) matrix of upper evaluation limits
##' @param lower (n,d) matrix of lower evaluation limits (<= upper)
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param mean.sqrt.mix E(sqrt(W)) or NULL; the latter if not available
##'        in which case it is estimated by QMC
##' @param loc d-vector (location vector)
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param control list() of algorithm parameters; see details in ?pnvmix
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'control$pnvmix.abstol'/'control$pnvmix.reltol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pnvmix <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d), qmix,
                   loc = rep(0, d), scale = diag(d), standardized = FALSE,
                   control = list(), verbose = TRUE, ...)
{
   ## Checks
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
   
   if(!is.matrix(scale)) scale <- as.matrix(scale)
   stopifnot(dim(lower) == c(n, d), length(loc) == d, # 'mean.sqrt.mix' is tested in pnvmix1()
             dim(scale) == c(d, d))
   
   ## Deal with algorithm parameters, see also ?get_set_param()
   control <- get_set_param(control)
   
   ## Grab method, increment and mean.sqrt.mix
   method        <- control$method
   increment     <- control$increment
   mean.sqrt.mix <- control$mean.sqrt.mix
   
   ## Absolute or relative precision required?
   tol <- if(is.null(control$pnvmix.abstol)) {
      ## Set tol to <0 so that algorithm runs until 'fun.eval[2]'
      ## function evaluations or 'max.iter.rqmc' iterations are exhausted
      do.reltol <- FALSE
      -42
   } else if(is.na(control$pnvmix.abstol)) { # if 'abstol = NA' use relative precision
      do.reltol <- TRUE
      control$pnvmix.reltol
   } else { # otherwise use absolute precision (default)
      do.reltol <- FALSE
      control$pnvmix.abstol
   }
   
   ## Define the quantile function of the mixing variable #######################
   special.mix <- NA
   qW <- if(is.character(qmix)) { # 'qmix' is a character vector
      qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma", "pareto"))
      switch(qmix,
             "constant" = {
                special.mix <- "constant"
                mean.sqrt.mix <- 1
                function(u) rep(1, length(u))
             },
             "inverse.gamma" = {
                if(hasArg(df)) {
                   df <- list(...)$df
                } else if(hasArg(nu)) {
                   nu <- list(...)$nu
                   df <- nu
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
                if(hasArg(alpha)) {
                   alpha <- list(...)$alpha
                } else if(hasArg(nu)){
                   nu <- list(...)$nu
                   alpha <- nu
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
   
   ## In the special case d = 1 we call pnvmix1d which is truly *vectorized*
   ## as there is no conditioning to be done. The case of a normal / t distribution
   ## will be handled correctly below (using pnorm() and pt())
   if(d == 1 && (is.na(special.mix) || special.mix == "pareto")) {
      return(pnvmix1d(upper = as.vector(upper), lower = as.vector(lower),
                      qW = qW, loc = loc, scale = as.numeric(scale),
                      standardized = standardized,
                      method = method,
                      tol = tol, do.reltol = do.reltol,
                      CI.factor = control$CI.factor,
                      fun.eval = control$fun.eval,
                      max.iter.rqmc = control$max.iter.rqmc,
                      increment = increment, B = control$B,
                      verbose = verbose))
   }
   ## Grab / approximate mean.sqrt.mix, which will be needed for preconditioning
   ## in pnvmix1(). This only depends on 'qmix', hence it is done (once) here in pnvmix.
   if(control$precond && d > 2) {
      if(is.null(mean.sqrt.mix))
         mean.sqrt.mix <- mean(sqrt(qW(qrng::sobol(n = 2^12, d = 1,
                                                   randomize = TRUE))))
      ## Check if provided/approximated mean.sqrt.mix is strictly positive
      if(mean.sqrt.mix <= 0)
         stop("'mean.sqrt.mix' has to be positive (possibly after being generated in pnvmix())")
   }
   ## Build temporary result list
   res1 <- vector("list", length = n) # results from calls of pnvmix1()
   ## Deal with NA
   NAs <- apply(is.na(lower) | is.na(upper), 1, any) # at least one NA => use NA => nothing left to do
   ## Remove 'loc':
   if(any(loc != 0)) {
      lower <- lower - matrix(rep(loc, n), ncol = d, byrow = TRUE)
      upper <- upper - matrix(rep(loc, n), ncol = d, byrow = TRUE)
   }
   ## Get (lower triangular) Cholesky factor of 'scale'
   factor.obj   <- cholesky_(scale, tol = control$cholesky.tol)
   factor       <- factor.obj$C
   factor.diag  <- factor.obj$D
   ## Obtain rank
   D.zero       <- (factor.diag==0)
   rank         <- d - sum(D.zero)
   ## In case of a singular matrix, need to reorder 'factor':
   if(rank < d) {
      if(verbose) warning("Provided 'scale' is singular.")
      ## In each row i, get minimal index j such that factor[i,k]=0 for all k>j
      length.rows <- apply(factor, 2, function(i) which.max( i != 0))
      order.length.rows <- order(length.rows) # needed later to sort 'low' and 'up'
      length.rows.sorted <- length.rows[order.length.rows]
      factor <- factor[order.length.rows, ] # sort factor
      lower <- lower[, order.length.rows, drop = FALSE]
      upper <- upper[, order.length.rows, drop = FALSE]
      ## Now factor has the form ( * non-zero element, ? any element)
      ##    * 0 0 ......... 0
      ##    ...
      ##    * 0 0 ......... 0
      ##    ? * 0 ......... 0
      ##    ...
      ##    ? * 0 ......... 0
      ##    ...
      ##    ? ... ? * 0 ... 0
      ## Need to standardize: Divide each row by its rightmost ' * ' (which is !=0)
      row.scales   <- factor[cbind(1:d, length.rows.sorted)] # vector of ' * ' elements
      factor       <- diag(1/row.scales) %*% factor
      ## Standardize 'lower' and 'upper' accordingly (also works for +/- Inf)
      lower <- matrix(rep(1/row.scales, each = n), ncol = d, nrow = n, byrow = FALSE) *
         lower
      upper <- matrix(rep(1/row.scales, each = n), ncol = d, nrow = n, byrow = FALSE) *
         upper
      ## Need to swap columns in 'lower'/'upper' multiplied by negative 'row.scales':
      if(any(row.scales < 0)) {
         which.row.scales.neg <- which(row.scales < 0)
         temp.lower <- lower[, which.row.scales.neg]
         lower[, which.row.scales.neg] <- upper[, which.row.scales.neg]
         upper[, which.row.scales.neg] <- temp.lower
      }
      ## Now we have the form as above, with '*' replaced by '1'. Remove 0 columns:
      factor <- factor[, -which(D.zero)]
      ## Get 'height' of each step:
      k.factor <- sapply(unique(length.rows.sorted), function(i) sum(length.rows == i))
      ## no update of'scale' as in singular case, no preconditioning is performed!
   } else {
      ## In case of non-singular 'scale', each of the d steps in factor has height 1:
      k.factor <- rep(1, d)
      ## Standardize 'scale', 'lower', 'upper', 'factor' ('loc' was already taken care of)
      if(!standardized) {
         row.scales     <- diag(factor) # diagonal of cholesky factor
         row.scales.inv <- 1/row.scales # d-vector
         ## The following is equivalent but faster than
         ## 'scale <- diag(row.scales.inv) %*% scale %*% diag(row.scales.inv)'
         ## [http://comisef.wikidot.com/tip:diagonalmatrixmultiplication]
         scale  <- outer(row.scales.inv, row.scales.inv) * scale # standardize scale
         factor <- matrix(rep(row.scales.inv, each = d),
                          ncol = d, byrow = TRUE) * factor # standardized cholesky
         ## Now standardize scale as above: Also works for +/- Inf
         lower <- matrix(rep(row.scales.inv, each = n), ncol = d, nrow = n, byrow = FALSE) *
            lower
         upper <- matrix(rep(row.scales.inv, each = n), ncol = d, nrow = n, byrow = FALSE) *
            upper
      }
   }
   
   ## Loop over observations ##################################################
   
   reached <- rep(TRUE, n) # indicating whether 'abstol' has been reached in the ith integration bounds (needs default TRUE)
   for(i in seq_len(n)) {
      if(NAs[i]) {
         res1[[i]] <- list(value = NA, error = 0, numiter = 0)
         next # => 'reached' already has the right value
      }
      ## Pick out ith row of lower and upper
      low <- lower[i,]
      up  <- upper[i,]
      ## Deal with equal bounds (result is 0 as integral over null set)
      if(any(low == up)) {
         res1[[i]] <- list(value = 0, error = 0, numiter = 0)
         next # => 'reached' already has the right value
      }
      ## Deal with case where components of both low and up are Inf
      lowFin <- is.finite(low)
      upFin  <- is.finite(up)
      lowupFin <- lowFin | upFin # at least one finite limit
      ## Only for full rank case for efficiency as ow rank, k.factor etc change
      ## completely and would need to be recalculated
      if(any(!lowupFin) && rank == d) {
         ## Update low, up
         low <- low[lowupFin]
         up  <- up [lowupFin]
         ## Grab (new) dimension. If 0, then all upper are +Inf, all lower are -Inf
         ## => Return 0
         d <- length(low) # Update dimension and 'k.factor'
         k.factor <- rep(1, d)
         if(d == 0) {
            res1[[i]] <- list(value = 1, error = 0, numiter = 0)
            next
         }
         ## Update scale (for precond)
         scale <- scale[lowupFin, lowupFin, drop = FALSE]
         factorFin <- t(chol(scale)) # Cholesky factor changes
      } else {
         ## If no such component exists, set Choleksy factor correctly
         factorFin <- factor
      }
      
      ## If d = 1, deal with multivariate normal, and t via pnorm() and pt()
      ## Note that everything has been standardized.
      if(d == 1 && !is.na(special.mix)) {
         if(special.mix == "constant") {
            value <- pnorm(up) - pnorm(low)
            res1[[i]] <- list(value = value, error = 0, numiter = 0)
            next
         } else if(special.mix == "inverse.gamma") {
            value <- pt(up, df = df) - pt(low, df = df)
            res1[[i]] <- list(value = value, error = 0, numiter = 0)
            next
         }
      }
      ## Compute result for ith row of lower and upper (in essence,
      ## because of the preconditioning, one has to treat each such
      ## row separately)
      res1[[i]] <- pnvmix1(up, lower = low, qW = qW,
                           is.const.mix = (!is.na(special.mix) &&
                                              special.mix == "constant"),
                           mean.sqrt.mix = mean.sqrt.mix,
                           loc = loc, scale = scale, factor = factorFin,
                           k.factor = k.factor,
                           method = method, precond = control$precond,
                           tol = tol, do.reltol = do.reltol,
                           CI.factor = control$CI.factor,
                           fun.eval = control$fun.eval,
                           max.iter.rqmc = control$max.iter.rqmc,
                           increment = increment,
                           B = control$B, verbose = as.logical(verbose),
                           inv.gam = (!is.na(special.mix) && special.mix == "inverse.gamma"),
                           ...)
      
      ## Check if desired precision was reached
      reached[i] <- res1[[i]]$error <= tol
      if(verbose >= 2 && !reached[i])
         warning(sprintf("Tolerance not reached for pair %d of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc'", i))
   } # for()
   if(verbose == 1 && any(!reached)) { # <=> verbose == 1
      ii <- which(!reached) # (beginning of) indices
      strng <- if(length(ii) > 6) {
         paste0(paste(head(ii), collapse = ", "), ",...")
      } else {
         paste(ii, collapse = ", ")
      }
      warning("Tolerance not reached for pair(s) ",strng," of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc'")
   }
   ## Return
   res <- vapply(res1, function(r) r$value, NA_real_)
   attr(res, "error")   <- vapply(res1, function(r) r$error, NA_real_)
   attr(res, "numiter") <- vapply(res1, function(r) r$numiter, NA_real_)
   res
}


##' @title Distribution Function of a 1d Normal Variance Mixture
##' @param upper n vector of upper evaluation limits
##' @param lower n vector of lower evaluation limits (<= upper)
##' @param qW function of one variable; quantile fct of W
##' @param loc numeric (location)
##' @param scale numeric (scale), >0
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param method see details in ?pnvmix
##' @param tol absolute/relative error tolerance
##' @param do.reltol logical if relative precision is required
##' @param CI.factor see details in ?pnvmix
##' @param fun.eval see details in ?pnvmix
##' @param max.iter.rqmc see details in ?pnvmix
##' @param increment see details in ?pnvmix
##' @param B see details in ?pnvmix
##' @param verbose see details in ?pnvmix
##' @param ... see details in ?pnvmix
##' @return numeric vector with the computed probabilities and attributes
##'         "error" (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pnvmix1d <- function(upper, lower = rep(-Inf,n), qW, loc = 0, scale = 1,
                     standardized = FALSE, method = "sobol",
                     tol = 1e-3, do.reltol = FALSE,
                     CI.factor = 3.5, fun.eval = c(2^6, 1e8), max.iter.rqmc = 15,
                     increment = "doubling", B = 15, verbose = FALSE)
{
   ## 1 Setup #################################################################
   
   n <- length(upper)
   if(!is.vector(upper)) upper <- as.vector(upper)
   if(!is.vector(lower)) lower <- as.vector(lower)
   ## Standardize
   if(!standardized) {
      upper <- (upper - loc) / sqrt(scale)
      lower <- (lower - loc) / sqrt(scale)
   }
   ## Error is calculated as CI.factor * sd( estimates) / sqrt(B).
   CI.factor.sqrt.B <- CI.factor / sqrt(B)
   ## Grab the number of sample points for the first iteration
   current.n <- fun.eval[1] #
   ## Matrix to store the B * length(upper) RQMC estimates
   rqmc.estimates <- matrix(0, ncol = n, nrow = B)
   ## Initialize error to > 'tol' so that we can enter the while loop below
   error <- tol + 42
   ## Initialize the total number of function evaluations
   total.fun.evals <- 0
   ## Initialize counter of the number of iterations in the while loop below
   numiter <- 0
   ## It may happen that qnorm(u) for u too close to 1 (or 0) is evaluated; in those
   ## cases, u will be replaced by ONE and ZERO which is the largest (smallest) number
   ## different from 1 (0) such that qnorm(u) is not +/- Inf
   ZERO <- .Machine$double.eps # for symmetry reasons (-8/+8), use this instead of .Machine$double.xmin
   ONE  <- 1-.Machine$double.neg.eps
   ## If method == sobol, we want the same random shifts in each iteration below,
   ## this is accomplished by reseting to the "original" seed
   if(method == "sobol") {
      if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
      seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used.
   }
   ## Additional variables needed if the increment chosen is "doubling"
   if(increment == "doubling") {
      if(method == "sobol") useskip <- 0
      denom <- 1
   }
   
   ## 2 Main while() loop #####################################################
   
   ## while() runs until precision 'tol' is reached or the number of function
   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
   ## the desired probability are calculated.
   while(max(error) > tol && total.fun.evals < fun.eval[2] && numiter < max.iter.rqmc)
   {
      if(method == "sobol" && numiter > 0)
         .Random.seed <<- seed # reset seed to have the same shifts in sobol(...)
      ## Get B RQCM estimates
      for(b in 1:B)
      {
         ## 2.1 Get the point set ###########################################
         
         U <- switch(method,
                     "sobol" = {
                        if(increment == "doubling") {
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
                        matrix(runif( current.n ), ncol = 1)
                     })
         ## U will contain realizations of 1 / sqrt(W):
         U <- 1 / cbind( sqrt(qW(U)), sqrt(qW(1-U)))
         
         ## 2.2 Evaluate the integrand at the (next) point set ##############
         
         next.estimate <- colMeans( (pnorm(outer( U[,1], upper)) -
                                        pnorm(outer(U[,1], lower)) +
                                        pnorm(outer(U[,2], upper)) -
                                        pnorm(outer(U[,2], lower)))/2)
         
         ## 2.3 Update RQMC estimates #######################################
         
         rqmc.estimates[b, ] <-
            if(increment == "doubling") {
               ## In this case both, rqmc.estimates[b] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b, ] + next.estimate) / denom
            } else {
               ## In this case, rqmc.estimates[b] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b, ] + next.estimate) / (numiter + 1)
            }
      } # end for(b in 1:B)
      
      ## Update of various variables
      ## Number of function evaluations
      ## (* 2 since antithetic variates are used in eval_nvmix_integral())
      total.fun.evals <- total.fun.evals + 2 * B * current.n
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(increment == "doubling") {
         ## Change denom and useksip (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      
      ## Update error depending on 'do.reltol'
      error <- if(!do.reltol) { # absolute error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)
      } else { # relative error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)/colMeans(rqmc.estimates)
      }
      numiter <- numiter + 1 # update counter
   } # while()
   
   ## 3 Finalize ##############################################################
   
   ## Check if error tolerance reached and print warnings accordingly
   reached <- (error <= tol)
   if(any(!reached) && verbose > 0) {
      ii <- which(!reached)
      if(verbose == 1) {
         strng <- if(length(ii) > 6) {
            paste0(paste(head(ii), collapse = ", "), ",...")
         } else {
            paste(ii, collapse = ", ")
         }
         warning("Tolerance not reached for pair(s) ",strng," of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument.")
      } else {
         for(i in 1:length(ii)) {
            warning(sprintf("Tolerance not reached for pair %d of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument", ii[i]))
         }
      }
   }
   ## Obtain estimates
   res <- colMeans(rqmc.estimates)
   attr(res, "error") <- error
   attr(res, "numiter") <- numiter
   
   ## Return
   res
}