### qnvmix() ##################################################################
##'
##' @param u vector of probabilities 
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param control see ?get.set.parameters()
##' @param verbose logical if warnings should be thrown
##' @param q.only if TRUE, only quantiles will be returned, ow additional quantites (see return )
##' @param stored.values matrix with 3 columns of the form [x, F(x), logf(x)] where
##'           F and logf are the df and log-density of the distribution specified
##'           in 'qmix'. 
##'           If provided it will be used to determine starting values for 
##'           the internal newton proceudure. Only very basic checking is done.   
##' @param ... see ?pnvmix()
##' @return if q.only is TRUE, vector of same length as u with entries F^{-1}(u)
##'         if q.only is FALSE, a list of 
##'         - $q: Vector of quantiles
##'         - $log.density: Vector log-density values at q
##'         - $computed.values: matrix with 3 columns [x, F(x), logf(x)] 
##'         - $newton.iterations: Vector, element i gives nb of 
##'            Newton iterations needed for u[i]
##' @author Erik Hintz, Marius Hofert and Christiane Lemieux  
##' @note - If only the quantiles are needed, abstol.logdensity does not need to be as small.


qnvmix <- function(u, qmix, control = list(), verbose = TRUE, q.only = TRUE,
                   stored.values = NULL, ...)
{
   quantile.internal(u, qmix = qmix, which = "nvmix1", d = 1, control = control,
                     verbose = verbose, q.only = q.only, 
                     stored.values = stored.values, ...)
}



### quantile.internal() ########################################################
##'
##' @param title Internal function used to estimate quantiles for either 
##'        univariate normal mixtures or for gamma mixtures 
##' @param u see ?qnvmix()
##' @param qmix see ?qnvmix()
##' @param which character string; either 'nvmix1' in which case the quantile
##'           of a 1dim nvmix dist'n is returned or 'maha2' in which
##'           case the quantile of the dist'n 'W * chi^2_d' is returned 
##' @param control see ?qnvmix()
##' @param verbose see ?qnvmix()
##' @param q.only see ?qnvmix()
##' @param stored.values see ?qnvmix()
##' @param ... see ?qnvmix()
##' @return see ?qnvmix()
##' @author Erik Hintz, Marius Hofert and Christiane Lemieux  

quantile.internal <- function(u, qmix, which = c('nvmix1', 'maha2'), d = 1,
                              control = list(), verbose = TRUE, 
                              q.only = FALSE, stored.values = NULL, ...) {
   
   ## 1 Checking and definitions ################################################
   
   ## Basic input checking (stored.values is checked below)
   stopifnot(!any(u>=1), !any(u<=0), is.logical(verbose), is.logical(q.only))
   if(!is.vector(u)) u <- as.vector(u) 
   ## Deal with algorithm parameters, see also get.set.parameters():
   ## get.set.parameters() also does argument checking, so not needed here. 
   control <- get.set.parameters(control)
   
   ## Grab method, B and n0
   method    <- control$method
   B         <- control$B
   n0        <- control$fun.eval[1]
   
   ## Define the quantile function of the mixing variable: 
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
   
   ## Build result vectors
   n               <- length(u)
   quantiles       <- rep(NA, n)
   log.density     <- rep(NA, n)
   num.iter.newton <- rep(0,  n)
   
   ## Deal with special distributions
   if(!is.na(special.mix)){
      if(!(special.mix == "pareto")){
         ## Only for "inverse.gamma" and "constant" do we have analytical forms:
         quantiles <- switch(special.mix,
                             "inverse.gamma" = {
                                if(which == "maha2"){
                                   ## D^2 ~ d* F(d, df)
                                   q <- qf(u, df1 = d, df2 = df) * d
                                   if(!q.only) log.density <- df(q/d, df1 = d, df2 = df, log = TRUE)-log(d)
                                   q
                                } else {
                                   ## X ~ t_df
                                   q <- qt(u, df = df)
                                   if(!q.only) log.density <- dt(q, df = df, log = TRUE)
                                   q
                                }
                             },
                             "constant" = {
                                if(which == "maha2"){
                                   ## D^2 ~ chi^2_d = Gamma(shape = d/2, scale = 2)
                                   q <- qgamma(u, shape = d/2, scale = 2)  
                                   if(!q.only) log.density <- dgamma(q, shape = d/2, scale = 2, log = TRUE)
                                   q
                                } else {
                                   ## X ~ N(0,1)
                                   q <- qnorm(u)
                                   if(!q.only) log.density <- dnorm(q, log = TRUE)
                                   q
                                }
                             })
         ## Return in those cases
         if(q.only) return(quantiles) else 
            return(list(q = quantiles, 
                        log.density = log.density, 
                        computed.values = cbind(quantiles, u, log.density, deparse.level = 0),
                        newton.iterations = rep(0, length(u))))
      }
   }
   
   ## 2 Set up est.cdf.dens() ## ################################################
   ## Initialize first pointset needed for rqmc approach:
   if(method == "sobol"){
      if(!exists(".Random.seed")) runif(1)
      seed <- .Random.seed
   } 
   
   ## Get realizations of W and sqrt(W)
   ## Initial point-set with B columns (each column = one shift)
   U0 <- switch(method,
                "sobol"   = {
                   sapply(1:B, function(i) 
                      sobol(n0, d = 1, randomize = TRUE))
                },
                "gHalton" = {
                   sapply(1:B, function(i) 
                      ghalton(n0, d = 1, method = "generalized"))
                },
                "prng"    = {
                   matrix(runif(B*n0), ncol = B)
                }) # (n0, B) matrix
   
   mixings      <- apply(U0, 2, qW) # qW() may be not 'matricized' 
   sqrt.mixings <- sqrt(mixings) # (n0, B) matrix
   
   ## Set up various quantities for est.cdf.dens():
   CI.factor <- control$CI.factor/sqrt(B) # instead of dividing by sqrt(B) all the time
   current.n <- n0 # will give ncol(mixings) 
   ldensity.constant <- switch(which,
                               "nvmix1" = {
                                  -1/2 * log(2* pi)
                               },
                               "maha2" = {
                                  -d/2*log(2) - lgamma(d/2)
                               }) # to avoid recomputation
   
   ## Function that estimates F(x) and logf(x) for *scalar* input x. Similar to
   ## pnvmix() and dnvmix() with increment = "num.init", tailored to the 
   ## one - dimensional case. Previous realizations of the mixing variable are
   ## reused (=> arguments mixings and sqrt.mixings)
   est.cdf.dens <- function(x, mixings, sqrt.mixings){
      ## Define various quantities:
      xx2 <- x*x/2 
      rqmc.estimates.log.density <- rep(NA, B) 
      rqmc.estimates.cdf         <- rep(NA, B)   
      current.n <- dim(mixings)[1]
      
      ## First we use 'mixings' and 'sqrt.mixings' that are already available.
      for(l in 1:B){
         ## Grab realizations corresponding to l'th shift and use exp-log trick 
         log.dens <- switch(which,
                            "nvmix1" = {
                               ldensity.constant - log(mixings[,l])/2 - 
                                  xx2 / mixings[,l] # length current.n
                            },
                            "maha2" = {
                               ldensity.constant - d/2*log(mixings[,l]) +
                                  (d/2-1)*log(x) - x/(2*mixings[,l])
                            })
         log.dens.max <- max(log.dens)
         rqmc.estimates.log.density[l] <- -log(current.n) + log.dens.max + 
            log(sum(exp(log.dens - log.dens.max)))
         rqmc.estimates.cdf[l] <- switch(which,
                                         "nvmix1" = {
                                            mean( pnorm(x/sqrt.mixings[,l]) )
                                         },
                                         "maha2" = {
                                            mean( pgamma(x/mixings[,l], 
                                                         shape = d/2, scale = 2))
                                         })
      }
      ## Check if precisions are reached
      error <- c(sd(rqmc.estimates.cdf)/mean(rqmc.estimates.cdf), 
                 sd(rqmc.estimates.log.density) ) *CI.factor
      precision.reached <- (error[1] <= control$newton.df.reltol & 
                               error[2] <= control$newton.logdens.abstol)
      if(!is.logical(precision.reached)) stop("Problem in est.cdf.dens(): error NA or NaN for x = ", paste(x))
      if(!precision.reached){
         ## Set up while loop
         iter.rqmc <- 1
         while(!precision.reached && iter.rqmc < control$max.iter.rqmc){
            ## Reset seed and get n0 realizations:
            if(method == "sobol") .Random.seed <<- seed
            
            U.next <- switch(method,
                             "sobol"   = {
                                sapply(1:B, function(i) 
                                   sobol(n0, d = 1, randomize = TRUE, skip = current.n))
                             },
                             "gHalton" = {
                                sapply(1:B, function(i) 
                                   ghalton(n0, d = 1, method = "generalized"))
                             },
                             "prng"    = {
                                matrix(runif(B*n0), ncol = B)
                             })
            mixings.next <- apply(U.next, 2, qW) # (n0, B) matrix
            sqrt.mixings.next <- sqrt(mixings.next ) # (n0, B) matrix
            
            ## Update RQMC estimators
            ## Note: Using C here did not really help, as x is only a scalar. 
            for (l in 1:B) {
               ## Grab realizations corresponding to l'th shift and use exp-log trick
               log.dens <- switch(which,
                                  "nvmix1" = {
                                     ldensity.constant - log(mixings.next[, l])/2 - 
                                        xx2 / mixings.next[, l] # length n0
                                  },
                                  "maha2" = {
                                     ldensity.constant - d/2*log(mixings.next[, l]) +
                                        (d/2-1)*log(x) - x/(2*mixings.next[, l])
                                  })
               log.dens.max <- max(log.dens)
               ## Previos estimate based on current.n samples, new one based on n0 samples
               
               rqmc.estimates.log.density[l] <- 
                  (current.n*rqmc.estimates.log.density[l] + 
                      n0*(-log(n0) + log.dens.max + log(sum(exp(log.dens - log.dens.max))))
                  )/(current.n + n0) 
               rqmc.estimates.cdf[l] <- 
                  switch(which,
                         "nvmix1" = {
                            (current.n * rqmc.estimates.cdf[l] + 
                                sum(pnorm(x/sqrt.mixings.next[,l])))/(current.n + n0)
                         },
                         "maha2" = {
                            (current.n * rqmc.estimates.cdf[l] + 
                                sum(pgamma(x/mixings.next[,l], 
                                           shape = d/2, scale = 2)))/(current.n + n0)
                         })
            }
            ## Update mixings and sqrt.mixings so that they can be reused:
            mixings <- rbind(mixings, mixings.next)
            sqrt.mixings <- rbind(sqrt.mixings, sqrt.mixings.next)
            current.n <- current.n + n0
            ## Update iteration number and error(s)
            iter.rqmc <- iter.rqmc + 1
            est.cdf   <- mean(rqmc.estimates.cdf)
            error <- c(sd(rqmc.estimates.cdf)/est.cdf, 
                       sd(rqmc.estimates.log.density) ) * CI.factor
            precision.reached <- (error[1] <= control$newton.df.reltol
                                  && error[2] <= control$newton.logdens.abstol)
         }
      }
      if(verbose){
         if(error[1] > control$newton.df.reltol) 
            warning("'newton.df.reltol' not reached; consider increasing 'max.iter.rqmc'
                          in the 'control' argument.")
         if(error[2] > control$newton.logdens.abstol && !q.only) 
            warning("'abstol.logdensity' not reached; consider increasing 'max.iter.rqmc'
                          in the 'control' argument.")
      }
      ## Return
      return(list(estimates = c(mean(rqmc.estimates.cdf), mean(rqmc.estimates.log.density)),
                  mixings = mixings, # return 'new' mxinings (= old mixings and new mixings)
                  sqrt.mixings = sqrt.mixings))
   }
   
   
   ## 3 Actual computation of quantiles (Newton's method) #######################
   
   if(which == "nvmix1"){
      ## Only compute quantiles for u>=0.5 (use symmetry in the other case)
      lower <- (u < 0.5)
      u[lower] <- 1 - u[lower]
   }
   ## Sort u. Ordering needed to return the quantiles in the correct order later
   ordering <- order(u)
   u.sorted <- u[ordering]
   ZERO     <- .Machine$double.neg.eps
   if(is.null(stored.values)){
      ## Matrix of the form [x, F(x), logf(x)] to store df and density evaluations
      ## Sample some mixings, transform them to realizations
      mixings. <- sample(mixings, n. <- min(100, max(n, 5)))
      x <- if(which == "nvmix1") c(0, mixings. * abs(rnorm(n.))) else 
         mixings. * rgamma(n., shape = d/2, scale = 2)
      x <- as.matrix(x)
      ## Evaluate df and density of the sample and store
      stored.values <- 
         switch(which,
                "nvmix1" = {
                   cbind(x, 
                         pnvmix(x, qmix = qW, scale = as.matrix(1), control = control),
                         dnvmix(x, qmix = qW, scale = as.matrix(1), control = control, log = TRUE))
                },
                "maha2" = {
                   cbind(x, 
                         pgammamix(x, qmix = qW, d = d, control = control),
                         dgammamix(x, qmix = qW, d = d, control = control, log = TRUE))
                })
      # stored.values <-  
      #    switch(which,
      #        "nvmix1" = {
      #           ## First element corresponds x = 0 so that F(x) = 0.5
      #           ## Note: est.cdf.dens(0)[1] = 0.5 (i.e. no error here), but need log-density
      #           cdf.dens.mixings <- est.cdf.dens(0, mixings, sqrt.mixings)
      #           matrix(c(0, cdf.dens.mixings$estimates), nrow = 1)
      #        },
      #        "maha2" = {
      #           ## First element corresponds x = d, need to estimate df and pdf
      #           ## (larger starting value to avoid numerical issues as density
      #           ## is very close to 0 for small x when d is large)
      #           cdf.dens.mixings <- est.cdf.dens(d, mixings, sqrt.mixings)
      #           matrix(c(d, cdf.dens.mixings$estimates), nrow = 1)
      #        })
      
   } else {
      ## Some very basic checking if stored.values was provided
      stopifnot(is.matrix(stored.values), dim(stored.values)[2] == 3,
                !any(stored.values[,2] > 1 | stored.values[,2] < 0))
   }
   
   ## Main loop for Newton's method:
   for (i in 1:n) {
      ## Initialize error and counter for Newton
      error <- control$newton.conv.abstol+ 42
      iter.newton <- 0
      ## Grab current u
      current.u <- u.sorted[i]
      ## Did we already have that u? (could happen if, e.g., u = c(0.4,0.6))
      if(i > 1 && current.u == u.sorted[i-1]){
         quantiles[i]   <- quantiles[i-1]
         log.density[i] <- log.density[i-1]
         next
      }
      ## Get starting value: x in stored.values sth F(x) is close to u
      closest.ind     <- which.min(abs(stored.values[, 2] - current.u))
      current.qu      <- stored.values[closest.ind, 1]
      current.funvals <- stored.values[closest.ind, 2:3] # (F(x), log f(x))
      
      ## Main loop for Newton procedure
      while (error > control$newton.conv.abstol && iter.newton < control$max.iter.newton) {
         ## Update quantile:
         next.qu <- current.qu - sign(current.funvals[1]-current.u)*
            exp(log( abs(current.funvals[1]-current.u)) - current.funvals[2])
         ## Quantiles > 0 here (since u>=0.5 for 'nvmix1')
         if(next.qu < ZERO) next.qu <- current.qu/2
         diff.qu <- (current.qu - (current.qu <- next.qu))
         ## Call est.cdf.dens:
         cdf.dens.mixings <- est.cdf.dens(current.qu, mixings, sqrt.mixings)
         current.funvals  <-  cdf.dens.mixings$estimates
         ## Store these values in stored.values
         stored.values    <- rbind( stored.values, c(current.qu, current.funvals),
                                    deparse.level = 0)
         ## Update mixings and sqrt.mixings
         mixings <- cdf.dens.mixings$mixings
         sqrt.mixings <- cdf.dens.mixings$sqrt.mixings
         ## Update error and increase counter
         error <- abs(diff.qu)
         iter.newton <- iter.newton + 1
      }
      ## Safe result
      quantiles[i]       <- current.qu
      log.density[i]     <- current.funvals[2]
      num.iter.newton[i] <- iter.newton
      if(verbose && error > control$newton.conv.abstol) 
         warning("'newton.conv.abstol' not reached; consider increasing 'max.iter.newton'
                 in the 'control' argument.")
   }
   
   ## 4 Clean-up and return   ###################################################  
   
   ## Order results according to original ordering of input u
   quantiles <- quantiles[order(ordering)]
   if(which == "nvmix1"){
      ## Use symmetry for those u which were < 0.5
      quantiles[lower] <- -quantiles[lower]
   }
   ## Return
   if(q.only) return(quantiles) else{
      num.iter.newton <- num.iter.newton[order(ordering)]
      log.density     <- log.density[order(ordering)]
      return(list(q = quantiles, 
                  log.density = log.density,
                  computed.values = stored.values,
                  newton.iterations = num.iter.newton))
   }
}