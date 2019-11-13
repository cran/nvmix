### get.set.parameters() #######################################################

##' @title  Retrieve algorithm specific default parameters and overwrite them
##' @return list with default values for all functions in the 'nvmix' package
##' @author Erik Hintz and Marius Hofert

get.set.parameters <- function(control = list()){
  ## Set up default controls:
  ctrl <- list(
    ## For pnvmix(): 
    mean.sqrt.mix = NULL, 
    precond = TRUE, 
    pnvmix.abstol = 1e-3, 
    pnvmix.reltol = NA, 
    cholesky.tol = 1e-9, 
    ## For dnvmix():
    dnvmix.abstol = 1e-3, 
    dnvmix.reltol = 1e-2, # If !NA, 'reltol' is used instead of 'abstol'
    dnvmix.max.iter.rqmc.pilot = 4,
    dnvmix.doAdapt = TRUE, 
    dnvmix.tol.int.lower = 1e-30,
    dnvmix.order.lower = 10,
    dnvmix.tol.bisec = c(1e-16, 1e-1, 1e-1),
    dnvmix.tol.stratlength = 1e-20,
    dnvmix.max.iter.bisec = 55,
    ## For pgammamix:
    pgammamix.reltol = NA,
    pgammamix.abstol = 1e-3,
    ## For qnvmix():
    max.iter.newton = 45, 
    newton.conv.abstol = 5e-4,
    newton.df.reltol = 5e-3,
    newton.logdens.abstol = 1e-2, 
    ## For fitnvmix():
    ### Algorithm specifications:
    ECMEstep = TRUE,
    ECMEstep.do.nu = TRUE,
    laststep.do.nu = FALSE,
    resample = FALSE, 
    ### Tolerances:
    ECME.maxiter = 20,
    ECME.miniter = 3,
    max.iter.locscaleupdate = 25,
    weights.abstol = 1e-1, # currently not used
    weights.reltol = 5e-2,
    weights.interpol.reltol = 1e-2,
    ECME.rel.conv.tol = c(1e-2, 1e-2, 5e-3), # [1] => 'loc'; [2] => 'scale'; [3] => 'nu'
    ### For the underlying 'optim':
    control.optim = list(maxit = 10),
    control.optim.laststep = list(), 
    ## For all (randomized) algorithms:
    method = "sobol", 
    increment = "doubling", # "doubling" or "num.init" 
    max.iter.rqmc = NA, # defined below, depending on 'increment'
    CI.factor = 3.5,
    fun.eval = c(2^7, 1e12), 
    B = 15,
    ## Additional returns for testing? (eg estimates after each iteration in
    ## 'fitnvmix')
    addReturns = FALSE)
  if(length(control) > 0){
    ## If input provided, grab input controls and overwrite:
    names.control <- names(ctrl)
    ctrl[(names.provided <- names(control))] <- control
    ## Did the user provide something that is not used?
    if (length(unmatched <- names.provided[!names.provided %in% names.control])) 
      warning("unknown names in control: ", paste(unmatched, collapse = ", "))
    ## If 'pnvmix.reltol' was provided, set 'pnvmix.abstol' to NA:
    if(any(names.provided == 'pnvmix.reltol') && !is.na(control$pnvmix.reltol))
       ctrl$pnvmix.abstol <- NA
    ## Check if 'method' and 'increment' were provided correctly
    ctrl$method     <- match.arg(ctrl$method, 
                                 choices = c("sobol", "ghalton", "PRNG"))
    ctrl$increment  <- match.arg(ctrl$increment, 
                                 choices = c("doubling", "num.init"))
    ## Now some more checkings: ('max.iter.rqmc' checked at the end)
    stopifnot(is.logical(ctrl$precond),
              ctrl$dnvmix.max.iter.rqmc.pilot >= 1,
              ctrl$dnvmix.tol.int.lower > 0,
              ctrl$dnvmix.tol.bisec.w >0,
              ctrl$dnvmix.tol.stratlength > 0,
              ctrl$dnvmix.max.iter.bisec.w > 0, 
              ctrl$max.iter.newton >= 0,
              ctrl$newton.conv.abstol >= 0,
              ctrl$newton.df.abstol >= 0,
              ctrl$newton.logdens.abstol >= 0,
              ctrl$weights.abstol >= 0,
              is.logical(ctrl$ECMEstep.do.nu),
              is.logical(ctrl$laststep.do.nu),
              ctrl$ECME.maxiter >= 0,
              length(ctrl$ECME.rel.conv.tol) == 3, ctrl$ECME.rel.conv.tol >= 0, 
              ctrl$CI.factor >= 0,
              length(ctrl$fun.eval) == 2, ctrl$fun.eval >= 0,
              ctrl$B > 1) # If B = 1 error estimates are NA => need B > 1
  }
  ## Define 'max.iter.rqmc': If it was not provided (=> NA), set defaults
  if(is.na(ctrl$max.iter.rqmc)){
    ctrl$max.iter.rqmc <- if(ctrl$increment == "doubling") 12 else 100
  } else {
    ## If it was provided (=> not NA), check if it's reasonable 
    stopifnot(ctrl$max.iter.rqmc > 1)
  }
  ## Return
  ctrl
}

