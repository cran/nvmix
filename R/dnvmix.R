### dnvmix() ###################################################################

##' @title Exp - log trick for log(sum_i exp(a_i))
##' @param M (n1, n2) matrix
##' @return n2-vector log(colSums(exp(M)))
##' @author Erik Hintz
##' @note NO checking is done for efficiency reasons
logsumexp <- function(M)
{
   cmax <- apply(M, 2, max)
   cmax + log(colSums(exp( M - rep(cmax, each = dim(M)[1]))))
}

##' @title Adaptive RQMC Method to Estimate the Log-Density of an NVMIX
##'        Distribution
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lconst see ?densmix_
##' @param d dimension of the Normal Variance Mixture
##' @param control see ?dnvmix
##' @param UsWs matrix; each row is of the form (u, qW(u)) where u in (0,1)
##' @return list of three ($ldensities, $error, $numiter); each vector of
##'         length(maha2.2)
##' @author Erik Hintz and Marius Hofert
densmix_adaptrqmc <- function(qW, maha2.2, lconst, d, k = d, control, UsWs)
{
   ## Initialization
   numiter           <- 0 # counter for the number of iterations
   ONE               <- 1-.Machine$double.neg.eps
   ZERO              <- .Machine$double.neg.eps
   tol.bisec         <- control$dnvmix.tol.bisec # vector of length 3
   ## [1]/[2]/[3] => tolerance on 'u' / 'W' / 'log-integrand' for bisections
   ## Define result object
   n          <- length(maha2.2)
   ldensities <- rep(NA, n)
   errors     <- rep(NA, n)
   numiters   <- rep(NA, n)
   ## Grab 'UsWs'; store them in vectors and sort
   ordering.U <- order(UsWs[, 1, drop = FALSE])
   ## Will store length(U/W/l.integrand) when appending elements in 'numObs'
   numObs     <- dim(UsWs)[1] + 2 # will add 'ZERO' and 'ONE'
   ## Set up matrix of the form  U | qW(U) | l.integrand
   U.W.lint       <- matrix(NA, ncol = 3, nrow = numObs)
   U.W.lint[, 1]  <- c(ZERO, UsWs[ordering.U, 1], ONE)
   U.W.lint[, 2]  <- qW(U.W.lint[,1])
   U.W.lint[1, 2] <- max(U.W.lint[1, 2], ZERO) # ensure qW(ZERO) > 0 
   ## Check if W is bounded
   W.max <- qW(1) # <= inf
   W.min <- qW(0) # >= 0
   isWbounded <- is.finite(W.max) | (W.min > 0)
   
   ## 2 Main loop over Mahalanobis distances ##################################
   
   for(ind in 1:n) {
      curr.maha2.2 <- max(maha2.2[ind], ZERO) # avoid maha2.2 = 0
      curr.lconst  <- lconst[ind]
      ## Initialize various quantities
      error         <- NA
      ldens.right   <- NA
      ldens.left    <- NA
      ldens.stratum <- NA
      uLuR          <- rep(NA, 2) # c('u.left', 'u.right') for later
      ## Determine location of maximum
      int.argmax.u  <- NA # u* such that qmix(u*) = W where integrand is max
      int.argmax.w  <- 2*curr.maha2.2/k # value of W = qmix(u*) at max assuming W is *unbounded*
      peekRight     <- (int.argmax.w > U.W.lint[numObs, 2]) && !isWbounded # u* close too close to 1
      peekLeft      <- (int.argmax.w < U.W.lint[1, 2]) && !isWbounded # u* close too close to 0
      outofreach.w  <- peekLeft ||  peekRight
      ## Deal with the case that W has support on (W.min, W.max) where W.min>0, W.max<inf
      if(maxLeft <- ((int.argmax.w < W.min) && isWbounded)) {
         ## Maximum at the left endpoint (ie at u=0)
         uLuR[1]       <- 0 # u.left = 0
         int.argmax.u  <- 0
         int.argmax.w  <- W.min
         ldens.left    <- -Inf
      }
      if(maxRight <- ((int.argmax.w > W.max) && isWbounded)) {
         ## Maximum at the right endpoint (ie at u=1)
         uLuR[2]       <- 1 # u.right = 1
         int.argmax.u  <- 1
         int.argmax.w  <- W.max
         ldens.right   <- -Inf
      }
      ## Realizations of log-integrand
      U.W.lint[,3]  <- curr.lconst - log(U.W.lint[,2])*k/2 -
         curr.maha2.2/U.W.lint[,2] # sorted according to ordering(U)!
      ## Value of the theoretical maximum of the log-integrand
      l.int.theo.max <- curr.lconst - log(int.argmax.w)*k/2 -
         curr.maha2.2/int.argmax.w
      ## Threshold: Will only use RQMC where l.integrand > l.tol.int.lower
      l.tol.int.lower <- min(max(log(control$dnvmix.tol.int.lower),
                                 l.int.theo.max - control$dnvmix.order.lower*log(10)), 0)
      
      ## 2.1 Find u* = argmax_u{integrand(u)} ################################
      
      ## Only needed if int.argmax.u is NA (otherwise, we already have it)
      if(!outofreach.w && is.na(int.argmax.u)) {
         ## Want fo find u* such that g(u*) ~= g_max where g is the original
         ## integrand. Equivalently, find u* such that qW(U*) = m/k
         ## (approximately) via binary search. The following finds starting
         ## values (u1, u2): qW(u1)<int.argmax.w & qW(u2) > int.argmax.w
         if(int.argmax.w < U.W.lint[1,2]) {
            ## "Peek" is close to u = 0
            u1u2 <- c(0, U.W.lint[1, 1])
         } else if(int.argmax.w > U.W.lint[numObs - 1, 2]) {
            ## "Peek" is close to u = 1
            u1u2 <- c(U.W.lint[numObs - 1, 1], 1)
         } else {
            ## "Peek" is observed
            ## First index i such that U.W.lint[i,2] (ie W[i]) >  int.argmax.w.
            ## Note: i > 1 bc of condition
            ind.first.bigger <- which(U.W.lint[,2] > int.argmax.w)[1]
            u1u2 <- c(U.W.lint[ind.first.bigger - 1, 1], U.W.lint[ind.first.bigger, 1])
         }
         ## Now find u:
         convd   <- FALSE
         numiter <- 0
         ## Matrix to store additional u's and W's generated:
         additionalVals <- matrix(NA, ncol = 2, nrow = control$dnvmix.max.iter.bisec)
         index.first.u  <- which(U.W.lint[,1] >= u1u2[1])[1] # will insert after this index
         while(!convd && numiter < control$dnvmix.max.iter.bisec) {
            numiter <- numiter + 1
            ## Next point to check (midpoint of u1u2):
            u.next <- mean(u1u2)
            diff <- ((w.next <- max(ZERO, qW(u.next))) - int.argmax.w)
            ## Store u.next and quanitle:
            additionalVals[numiter, ] <- c(u.next, w.next)
            ## Update u1u2 depending on sign of 'diff' and check convergence:
            if(diff > 0) u1u2[2] <- u.next else u1u2[1] <- u.next
            convd <- (abs(diff) < tol.bisec[2]) || (diff(u1u2) < tol.bisec[1])
         }
         int.argmax.u <- u.next
         ## Add additional values of 'U', 'W' and 'l.integrand' in 'U.W.lint'
         ## while *preserving the order*.
         ## Order of first column = order of all columns
         order.additional.Us <- order(additionalVals[1:numiter, 1, drop = FALSE])
         WsToAdd <- additionalVals[order.additional.Us, 2, drop = FALSE] # needed twice
         U.W.lint <- rbind( U.W.lint[1:index.first.u,],
                            cbind(additionalVals[order.additional.Us, 1, drop = FALSE],
                                  WsToAdd,
                                  curr.lconst - log(WsToAdd)*k/2 - curr.maha2.2/WsToAdd),
                            U.W.lint[(index.first.u+1):numObs,])
         ## Update length (=nrow) of 'U.W.lint'
         numObs <- numObs + numiter
      } else if(peekLeft && is.na(int.argmax.u)) {
         int.argmax.u <- ZERO
      } else if(peekRight && is.na(int.argmax.u)) {
         int.argmax.u <- ONE
      }
      ## 2.2 Find stratum (u.left, u.right) over which RQMC is applied #######
      
      ## 2.2.1 Find "candidates" for bisection ###############################
      
      ## Need to distinguish multiple special cases:
      greater.thshold <- (U.W.lint[,3] > l.tol.int.lower)
      if(any(greater.thshold)) {
         ## In this case there is at least one obs exceeding the threshold
         ind.greater <- which(greater.thshold)
         ## Will use bisection to find 'u.left' and 'u.right'; starting
         ## values as follows:
         candid.left  <- c(max(0, U.W.lint[ind.greater[1] - 1, 1]),
                           min(int.argmax.u, U.W.lint[ind.greater[1], 1]) )
         ## Note: ind.greater[1] = 1 => U[ind.greater[1] - 1] = numeric(0)
         ## => max(0, U[ind.greater[1] - 1]) = 0
         ## Need several cases for 'candid.right':
         candid.right <- if(ind.greater[length(ind.greater)] < numObs) {
            ## Let u^ be largest u >= u* such that g(u^) > threshold (exists)
            ## Case 1:
            ## There is u' with u'>=u^ and g(u') < threshold
            ## =>  'u.right' in (u^, u')
            c(max(U.W.lint[ind.greater[length(ind.greater)], 1], int.argmax.u),
              U.W.lint[ind.greater[length(ind.greater)] + 1, 1] )
         } else {
            ## Case 2: the last point in U is such that g(u) > threshold, in other words:
            ## g(ONE) > threshold
            ## => u.right = 1 no matter where the peek is, but the handling of the region
            ## (ONE, 1) *does* depend on where the peek is.
            ldens.right <- if(peekRight) {
               ## Area with peek right of ONE => can't simulate there
               ## => stratify from u.left to u.right = 1 (=ONE)
               ## => use *crude* approximation for region (ONE, 1) where max is
               log(.Machine$double.neg.eps) + l.int.theo.max - log(2)
            } else {
               ## Area with the peek is left of ONE but theoretical 'u.right' is
               ## out of reach due to machine precision.
               ## => stratify from u.left to u.right = 1 (=ONE)
               ## => do nothing in  (ONE, 1)
               -Inf
            }
            ## See above: u.right = 1
            c(1, 1)
         }
      } else if(peekRight) {
         ## In this case, none of the observations is > threshold and peek close to 1:
         ## Will use *all* obs to estimate in (0, ONE) and crude approximation in (ONE, 1)
         ## where max is. NO stratification!
         candid.left    <- c(1, 1)
         candid.right   <- c(1, 1)
         ldens.right    <- log(.Machine$double.neg.eps) + l.int.theo.max - log(2)
         ldens.stratum  <- -Inf
      } else if(peekLeft) {
         candid.left    <- c(0, 0)
         candid.right   <- c(0, 0)
         ldens.left     <- log(.Machine$double.neg.eps) + l.int.theo.max - log(2)
         ldens.stratum <- -Inf
      } else {
         candid.left  <- c(0, int.argmax.u)
         candid.right <- c(int.argmax.u, 1)
      }
      
      ## 2.2.2 Bisection to find 'u.left' and 'u.right' ######################
      
      candids <- rbind(candid.left, candid.right)
      for(i in 1:2) {
         curr.candid <- candids[i, ]
         if(!is.na(uLuR[i])) next # we already set u.left or u.right
         uLuR[i] <- if(diff(curr.candid) <= tol.bisec[1]) {
            curr.candid[1]
         } else {
            ## Use bisection similar to the one used to find u*
            convd <- FALSE
            numiter <- 0
            ## Matrix to store additional u's, W's, l.integrand's generated:
            additionalVals <- matrix(NA, ncol = 3,
                                     nrow = control$dnvmix.max.iter.bisec)
            while(!convd && numiter < control$dnvmix.max.iter.bisec) {
               numiter <- numiter + 1
               ## Next point to check:
               u.next <- mean(curr.candid)
               w.next <- max(qW(u.next), ZERO)
               diff   <- (l.int.next <- (curr.lconst -log(w.next)*k/2 -
                                            curr.maha2.2/w.next)) - l.tol.int.lower
               ## Store values generated
               additionalVals[numiter, ] <- c(u.next, w.next, l.int.next)
               ## Update 'curr.candid' depending on sign of 'diff' and check convergence:
               if(diff > 0) curr.candid[2] <- u.next else 
                  curr.candid[1] <- u.next
               convd <- (abs(diff) < tol.bisec[3])  || (diff(curr.candid) < tol.bisec[1])
            }
            ## Update U.W.lint[]:
            ## First, add additional values
            U.W.lint <- rbind(U.W.lint, additionalVals[1:numiter,, drop = FALSE])
            ## Destroyed the ordering => sort again
            ordering.new <- order(U.W.lint[,1]) # all columns have the same ordering
            U.W.lint <- U.W.lint[ordering.new, ]
            ## Update length (=nrow) of 'U.W.lint'
            numObs <- numObs + numiter
            u.next
         }
      }
      ## For readability
      u.left  <- if(uLuR[1] <= ZERO) 0 else uLuR[1]
      u.right <- if(uLuR[2] >= ONE) 1 else uLuR[2]
      
      ## 2.3 Estimate integral in the three regions ##########################
      
      ## 2.3.1 Estimate log of integral over (0, u.left)
      ldens.left <- if(is.na(ldens.left)) {
         ## 'ldens.left' not set yet => peek is *not* in (0, u.left) => Riemann sums
         if(u.left == 0) {
            -Inf
         } else if(u.left == 1) {
            ## => Special case where no obs > threshold
            ## => Use all obs and Riemann for this region
            weights <- c(U.W.lint[1, 1], U.W.lint[2:numObs, 1] -
                            U.W.lint[1:(numObs-1), 1])
            upper.sum <- logsumexp(
               as.matrix(log(weights) + U.W.lint[1:numObs, 3], ncol = 1))
            lower.sum <- logsumexp(
               as.matrix(log(weights) + c(-Inf, U.W.lint[1:(numObs-1), 3]),
                         ncol = 1))
            (lower.sum + upper.sum) / 2
         } else {
            ## 0 < u.left < 1 => Find obs in (0, u.left)
            usUsed <- (U.W.lint[,1] <= u.left)
            sumusUsed <- sum(usUsed)
            if(sumusUsed > 1) {
               ## Case 1: We have >1 observations in (0, u.left)
               lastUsed <- which(usUsed)[sumusUsed]
               weights <- c(U.W.lint[1, 1],
                            U.W.lint[2:lastUsed, 1] -
                               U.W.lint[1:(lastUsed-1), 1])
               upper.sum <- logsumexp(as.matrix(log(weights) +
                                                   U.W.lint[1:lastUsed, 3],
                                                ncol = 1))
               lower.sum <- logsumexp(as.matrix(log(weights) +
                                                   c(-Inf, U.W.lint[1:(lastUsed-1), 3]),
                                                ncol = 1))
               (lower.sum + upper.sum) / 2
            } else {
               ## Case 2: No observations in (u.right, 1) => u.right > ONE
               log(1 - u.left) + l.tol.int.lower - log(2)
            }
         }
      } else {
         ldens.left
      }
      
      ## 2.3.2 Estimate log of integral over (u.right, 1)
      ldens.right <- if(is.na(ldens.right)) {
         ## 'ldens.right' not set yet => peek is *not* in (u.right, 1) => Riemann sums
         if(u.right == 1) {
            -Inf
         } else if(u.right == 0) {
            ## => Special case where no obs > threshold
            ## => Use all obs and Riemann for this region
            weights <- c(U.W.lint[1, 1],
                         U.W.lint[2:numObs, 1] -
                            U.W.lint[1:(numObs-1), 1])
            upper.sum <- logsumexp(
               as.matrix(log(weights) + U.W.lint[1:numObs, 3],
                         ncol = 1))
            lower.sum <- logsumexp(
               as.matrix(log(weights) + c(-Inf, U.W.lint[1:(numObs-1), 3]),
                         ncol = 1))
            (lower.sum + upper.sum) / 2
         } else {
            ## 0 < u.right < 1 => Find obs in (u.right, 1)
            usUsed <- (U.W.lint[,1] >= u.right)
            if(any(usUsed)) { # maybe redundnat, see above
               ## Case 1: We have observations in (u.right, 1)
               firstUsed <- which(usUsed)[1]
               weights <- c(U.W.lint[ (firstUsed+1):numObs, 1] -
                               U.W.lint[firstUsed:(numObs-1), 1],
                            .Machine$double.neg.eps)
               upper.sum <- 
                  logsumexp(as.matrix(log(weights) + U.W.lint[firstUsed:numObs, 3],
                                      ncol = 1))
               lower.sum <- 
                  logsumexp(as.matrix(
                     log(weights) + c(U.W.lint[(firstUsed+1):numObs, 3], -Inf),
                                                ncol = 1))
               (lower.sum + upper.sum) / 2
            } else {
               ## Case 2: No observations in (u.right, 1) => u.right > ONE
               log1p(-u.right) + l.tol.int.lower - log(2)
            }
         }
      } else {
         ldens.right
      }
      
      ## 2.3.3 Estimate log of integral over (u.left, u.right)
      ## We use RQMC in this region, unless 'ldens.strat' was already set.
      stratlength <- u.right - u.left
      rqmc.numiter <- 0
      ldens.stratum <- if(is.na(ldens.stratum)) {
         if(stratlength > control$dnvmix.tol.stratlength) {
            ldens.obj <- densmix_rqmc(qW, maha2.2 = curr.maha2.2,
                                      lconst = curr.lconst,
                                      d = d, k = k, control = control,
                                      lower.q = u.left, upper.q = u.right,
                                      return.all = FALSE,
                                      max.iter.rqmc = control$max.iter.rqmc -
                                         control$dnvmix.max.iter.rqmc.pilot)
            error        <- ldens.obj$error
            rqmc.numiter <- ldens.obj$numiter
            ldens.obj$ldensities + log(stratlength)
         } else {
            log(max(stratlength, ZERO)) + l.int.theo.max - log(2)
         }
      } else {
         ldens.stratum
      }
      
      ## 3 Combine and return ################################################
      
      ## Combine to one estimate:
      ldensities[ind] <- logsumexp(rbind(ldens.left, ldens.right, ldens.stratum,
                                         deparse.level = 0))
      errors[ind]     <- error
      numiters[ind]   <- rqmc.numiter
   }
   list(ldensities = ldensities, error = errors, numiter = numiters)
}


##' @title Density of a Multivariate Normal Variance Mixture for restricted W
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lrdet log(sqrt(det(scale))) where 'scale' is the scale matrix of
##'        the normal variance mixture distribution.
##' @param d dimension of the Normal Variance Mixture
##' @param k power of qW^{-1} in the root (k=d => density, k=d+2 =>unnormalized
##'        conditional expectation of 1/W)
##' @param control see ?dnvmix
##' @param lower.q numeric in (0,1)
##' @param upper.q numeric in (0,1), > lower.q. Density will be estimated
##'         conditional on W being between its lower.q and upper.q quantile.
##' @param max.iter.rqmc maximum number of iterations
##' @param return.all logical; if true, matrix (U, qW(U)) also returned.
##' @return List of three:
##'         $ldensities n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either
##'         relative error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms
##'         (only if return.all = TRUE)
##' @author Erik Hintz and Marius Hofert
densmix_rqmc <- function(qW, maha2.2, lconst, d, k = d, control,
                         lower.q = 0, upper.q = 1,
                         max.iter.rqmc, return.all)
{
   ## Define various quantites:
   dblng           <- (control$increment == "doubling")
   B               <- control$B # number of randomizations
   n               <- length(maha2.2) # sample size
   current.n       <- control$fun.eval[1] #initial sample size
   numiter         <- 0 # counter for the number of iterations
   total.fun.evals <- 0
   ZERO            <- .Machine$double.neg.eps
   ## Absolte/relative precision?
   if(is.na(control$dnvmix.reltol)) {
      ## Use absolute error
      tol <- control$dnvmix.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error
      tol <- control$dnvmix.reltol
      do.reltol <- TRUE
   }
   ## Store seed if 'sobol' is used to get the same shifts later:
   if(control$method == "sobol") {
      if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
      seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is used
   }
   ## Additional variables needed if the increment chosen is "dblng"
   if(dblng) {
      if(control$method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Matrix to store RQMC estimates
   rqmc.estimates <- matrix(-Inf, ncol = n, nrow = B)
   ## Will be needed a lot:
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ## Define trafo-function that maps u to (q,1) or (1,q) depending on 'up'
   trafo <- function(u) lower.q + (upper.q - lower.q)*u
   ## Initialize 'max.error' to > tol so that we can enter the while loop:
   max.error <- tol + 42
   ## Matrix to store U, W values => nrows = maximal number of funevals
   if(return.all) {
      max.nrow <- if(dblng) current.n*B*2^(max.iter.rqmc-1) else 
         max.iter.rqmc*B*current.n
      UsWs <- matrix(NA, ncol = 2, nrow = max.nrow)
      curr.lastrow <- 0 # will count row-index additional points are being inserted after
   }
   
   ## Main loop ###############################################################
   
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
         W <- qW(U <- trafo(U)) # realizations of the mixing variable; sorted!
         ## Need to replace values < ZERO by ZERO. W is *sorted*, so check using
         ## loop instd of 'pmax' (more efficient)
         for(ind in 1:current.n) if(W[ind] < ZERO) W[ind] <- ZERO else break
         ## Update 'UsWs'
         if(return.all) {
            UsWs[(curr.lastrow + 1) : (curr.lastrow + current.n), ] <- cbind(U, W)
            curr.lastrow <- curr.lastrow + current.n
         }
         next.estimate <- .Call("eval_densmix_integrand",
                                W          = as.double(W),
                                maha2_2    = as.double(maha2.2),
                                current_n  = as.integer(current.n),
                                n          = as.integer(n),
                                d          = as.integer(d),
                                k          = as.integer(k),
                                lconst     = as.double(lconst))
         
         ## Update RQMC estimates
         rqmc.estimates[b,] <-
            if(dblng) {
               ## In this case both, rqmc.estimates[b,] and
               ## next.estimate depend on n.current points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b,]),
                     b = as.double(next.estimate),
                     n = as.integer(n)) - log(denom)
            } else {
               ## In this case, rqmc.estimates[b,] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b,] + log(numiter)),
                     b = as.double(next.estimate),
                     n = as.integer(n)) - log(numiter + 1)
            }
      } # end for(b in 1:B)
      ## Update of various variables
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(dblng) {
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
      ## Total number of function evaluations
      total.fun.evals <- total.fun.evals + B * current.n
      numiter <- numiter + 1
      ## Update error. The following is slightly faster than 'apply(..., 2, var)'
      ldensities <- logsumexp(rqmc.estimates) - log(B) # performs better than .colMeans
      vars <- .colMeans((rqmc.estimates - rep(ldensities, each = B))^2, B, n, 0)
      errors <- if(!do.reltol) { # absolute error
         sqrt(vars)*CI.factor.sqrt.B
      } else { # relative error
         sqrt(vars)/abs(ldensities)*CI.factor.sqrt.B
      }
      max.error <- max(errors)
   } # while()
   
   ## Return
   if(return.all) {
      list(ldensities = ldensities, numiter = numiter, error = errors,
           UsWs = UsWs[1:curr.lastrow,])
   } else {
      list(ldensities = ldensities, numiter = numiter, error = errors)
   }
}


##' @title Integration Routine of log h(u) on the unit interval
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lconst vector of same length as 'maha2.2', see above
##' @param d dimension of the Normal Variance Mixture
##' @param control see ?get_set_param()
##' @param verbose see ?dnvmix()
##' @return List of three:
##'         $ldensities n-vector with computed log-density values
##'         $error n-vector of *absolute* error estimates for log-densities
##'         $numiter n-vector of number of iterations needed
##' @author Erik Hintz and Marius Hofert
##' @note Estimates int_0^1 h(u) du where log h(u) = lconst - d/2 log qW(u) - maha2.2/qW(u).
##' For density of 'nvmix': lconst = rep(-lrdet -d/2log(2pi), length(maha2.2));
##' For density of 'gammamix': lconst = -lgamma(d/2) - d/2 log(2) + (d/2-1)log(2*maha2.2)
densmix_ <- function(qW, maha2.2, lconst, d, control, verbose)
{
   ## Absolte/relative precision?
   tol <- if(is.na(control$dnvmix.reltol)) { # if 'reltol = NA' use absolute precision
      do.reltol <- FALSE
      control$dnvmix.abstol
   } else { # otherwise use relative precision (default)
      do.reltol <- TRUE
      control$dnvmix.reltol
   }
   
   ## Call RQMC procedure without any stratification
   rqmc.obj <- densmix_rqmc(qW, maha2.2 = maha2.2, lconst = lconst, d = d,
                            control = control, lower.q = 0, upper.q = 1,
                            max.iter.rqmc = control$dnvmix.max.iter.rqmc.pilot,
                            return.all = TRUE)
   ## Extract results
   ldens   <- rqmc.obj$ldensities
   numiter <- rep(rqmc.obj$numiter, length(maha2.2))
   error   <- rqmc.obj$error
   if(any(error > tol)) {
      ## Accuracy not reached for at least one 'maha2.2' value
      ## => Use adaptive approach for those
      if(control$dnvmix.doAdapt) {
         notRchd <- which(error > tol)
         rqmc.obj <- densmix_adaptrqmc(qW, maha2.2 = maha2.2[notRchd],
                                       lconst = lconst[notRchd],
                                       d = d, UsWs = rqmc.obj$UsWs,
                                       control = control)
         ldens[notRchd]    <- rqmc.obj$ldensities
         numiter[notRchd]  <- numiter[notRchd] + rqmc.obj$numiter
         error[notRchd]    <- rqmc.obj$error
      }
      ## Handle warnings:
      if(as.logical(verbose)) { # verbose can be passed by 'fitnvmix' => between 0 and 3
         if(any(is.na(error))) {
            ## At least one error is NA
            warning("Estimation unreliable, corresponding error estimate NA")
         }
         whichNA <- which(is.na(error))
         if(any(error[setdiff(1:length(error), whichNA)] > tol)) # 'setdiff' needed if 'whichNA' is empty
            warning("Tolerance not reached for all inputs;
                consider increasing 'max.iter.rqmc' in the 'control' argument.")
      }
      ## Transform error back to *absolute* errors:
      if(do.reltol) error <- error * abs(ldens)
   }
   list(ldensities = ldens, numiter = numiter, error = error)
}


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
##'         - "PRNG":    pure Monte Carlo
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
##' @param max.iter.rqmc maximum number of iterations in the RQMC approach
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
                   control = list(), log = FALSE, verbose = TRUE,...)
{
   ## Checks
   if(!is.matrix(x)) x <- rbind(x)
   d <- ncol(x) # dimension
   if(!is.matrix(scale)) scale <- as.matrix(scale)
   stopifnot(length(loc) == d, dim(scale) == c(d, d))
   verbose <- as.logical(verbose)
   ## Deal with algorithm parameters, see also get_set_param():
   ## get_set_param() also does argument checking, so not needed here.
   control <- get_set_param(control)
   ## If factor is not provided, determine it here as a *lower* triangular matrix
   if(is.null(factor)) factor <- t(chol(scale)) # lower triangular
   
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
                } else if(hasArg(nu)){
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
   ## Build result object (log-density)
   n <- nrow(x)
   lres <- rep(-Inf, n) # n-vector of results
   notNA <- rowSums(is.na(x)) == 0
   lres[!notNA] <- NA
   x <- x[notNA,, drop = FALSE] # non-missing data (rows)
   
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
   ## Counter
   numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
   ## Deal with the different distributions
   if(!is.na(special.mix)) {
      lres[notNA] <- switch(special.mix,
                            "inverse.gamma" = {
                               lgamma((df + d) / 2) - lgamma(df/2) - (d/2) * 
                                  log(df*pi) - lrdet - (df+d)/2 * log1p(maha2/df)
                            },
                            "constant" = {
                               -(d/2) * log(2 * pi) - lrdet - maha2/2
                            },
                            "pareto" = {
                               log(alpha) - d/2*log(2*pi) - lrdet - 
                                  (alpha+d/2)*log(maha2/2) +
                                  pgamma(maha2/2, scale = 1, shape = alpha+d/2, 
                                         log.p = TRUE) + lgamma(alpha+d/2)
                            })
      if(!log) lres <- exp(lres) # already exponentiate
      error <- rep(0, length(maha2))
   } else {
      ## General case of a multivariate normal variance mixture (RQMC)
      ## Prepare inputs for densmix_rqmc
      ## Sort maha-distance and divide by 2; store ordering to recover original
      ## ordering later:
      ordering.maha <- order(maha2)
      maha2.2 <- maha2[ordering.maha]/2
      ## Define log-constant for the integration
      lconst <- rep(-lrdet - d/2*log(2*pi) , length(maha2.2))
      ## Call internal dnvmix (which itself calls C-Code and handles warnings)
      ests <- densmix_(qW, maha2.2 = maha2.2, lconst = lconst,
                       d = d, control = control, verbose = verbose)
      ## Grab results, correct 'error' and 'lres' if 'log = FALSE'
      lres[notNA] <- ests$ldensities[order(ordering.maha)]
      error <- if(log) {
         ests$error[order(ordering.maha)]
      } else {
         lres <- exp(lres)
         ests$error[order(ordering.maha)]*pmax(lres[notNA], 1)
      }
      numiter <- ests$numiter
   }
   
   ## Return
   ## Note that 'lres' was exponentiated already if necessary.
   attr(lres, "error")   <- error # these are absolute errors, no matter what!
   attr(lres, "numiter") <- numiter
   lres
}
