### QQ-plot of Mahalanobis distances for visual GOF ############################

##' @title QQ Plot of Mahalanobis distances versus their theoretical quantiles
##' @param x (n,d) data matrix
##' @param qmix see ?pnvmix()
##' @param loc see ?pnvmix()
##' @param scale see ?pnvmix()
##' @param control see ?pnvmix()
##' @param plot.diag logical; if TRUE the curve f(x) = x is plotted additionally
##' @param verbose logical if warnings from underlying 'qgammamix()' shall be
##'        thrown
##' @param control list of algorithm specific parameters, see ?get_set_param()
##'        and ?fitnvmix
##' @param ... see ?pnvmix()
##' @return invisibly returns a list of two: 'maha2' (sorted squared mahalanobis
##'         distances obtained from 'x', sorted) and 'q' (theoretical quantiles)
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
qqplot_maha <- function(x, qmix, loc, scale, plot.diag = TRUE, verbose = TRUE,
                        control = list(), ...)
{
   ## Initialize and check inputs
   control <- get_set_param(control)
   control$newton.df.reltol <- control$qqplot.df.reltol
   if(!is.matrix(x)) x <- rbind(x)
   notNA <- rowSums(is.na(x)) == 0
   x     <- x[notNA,, drop = FALSE] # non-missing data (rows)
   n     <- nrow(x)
   d     <- ncol(x)
   ## Obtain sorted Mahalanobis distances (X-loc)^T scale^{-1} (X-loc)
   maha2 <- sort(mahalanobis(x, center = loc, cov = scale))
   ## Obtain theoretical quantiles
   theoretical.quantiles <-
      qgammamix(ppoints(n), qmix = qmix, d = d, control = control,
                verbose = verbose, q.only = TRUE, stored.values = NULL, ...)
   ## Plot
   plot(theoretical.quantiles, maha2, xlab = "Theoretical quantiles",
        ylab = "Sample quantiles", main = "")
   ## Add a diagonal
   if(plot.diag) lines(theoretical.quantiles, theoretical.quantiles, lty = 2)
   invisible(list(maha2 = maha2, q = theoretical.quantiles))
}
