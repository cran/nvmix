\name{get.set.parameters}
\alias{get.set.parameters}
\title{Algorithm specific parameters for functionalities in the nvmix package}
\description{
  Algorithm specific parameters for functionalities in the \code{nvmix} package,
  notably for \code{\link{fitnvmix}}, \code{\link{dnvmix}} and \code{\link{pnvmix}}.
}
\usage{
get.set.parameters(control = list())
}
\arguments{
  \item{control}{\code{\link{list}} specifying algorithm specific parameters to be set;
  see below under details.}
}
\value{
  \code{get.set.parameters()} returns a \code{\link{list}} with more than 30 elements 
  specifying algorithm specific parameters for the functions \code{\link{fitnvmix}}, 
  \code{\link{dnvmix}} and \code{\link{pnvmix}}. 
  Parameter values passed to \code{get.set.parameters()} via the
  \code{control} argument overwrite the defaults; for parameters not specified in the
  \code{control} argument, the default values are being returned. 
}
\details{
  For most functions in the \code{nvmix} package, internally, an iterative randomized 
  Quasi-Monte Carlo (RQMC) approach is used to estimate probabilities, weights and 
  (log-)densities. There are various parameters than can be individiually
  changed. 
  
  Algorithm specific parameters can be passed as a list via \code{control}. 
  It can contain any of the following:
  
  \describe{
  \item{For all algorithms:}{
      \describe{
         \item{\code{method}}{\code{\link{character}} string indicating the method 
         to be used to compute the integral. Available are:
            \describe{
               \item{\code{"sobol"}:}{Sobol' sequence (default).}
               \item{\code{"ghalton"}:}{generalized Halton sequence.}
               \item{\code{"PRNG"}:}{plain Monte Carlo based on a pseudo-random 
               number generator.}
            }
         }
         \item{\code{increment}}{\code{\link{character}} string indicating how the 
         sample size should be increased in each iteration. Available are:
            \describe{
               \item{\code{"doubling"}:}{next iteration has as many sample points 
               as all the previous iterations combined.}
               \item{\code{"num.init"}:}{all iterations use an additional 
               \code{fun.eval[1]}-many points.}
            }
         }
         \item{\code{CI.factor}}{multiplier of the Monte Carlo confidence interval 
         bounds. The algorithm runs until \code{CI.factor} times the estimated 
         standard error is less than \code{abstol} or \code{reltol} (whichever is 
         provided). If \code{CI.factor = 3.3} (the default), one can expect the 
         actual absolute error to be less than \code{abstol} in 99.9\% of the cases.}
         \item{\code{fun.eval}}{\code{\link{numeric}(2)} providing the size of the 
         first point set to be used to estimate integrals  
         (typically a power of 2) and the maximal number of function 
         evaluations. \code{fun.eval} defaults to \code{c(2^7, 1e8)}.}  
         \item{\code{max.iter.rqmc}}{\code{\link{numeric}}, providing the maximum 
         number of iterations allowed in the RQMC approach; the default is 15 if
         \code{increment = "doubling"} and 100 otherwise.}  
         \item{\code{B}}{number of randomizations for obtaining an error estimate 
         in the randomized quasi-Monte Carlo (RQMC) approach; the default is 15.}
      }   
  }
  \item{For \code{pnvmix()}:}{
      \describe{
         \item{\code{pnvmix.abstol}, \code{pnvmix.reltol}}{non-negative numeric 
         providing the relative/absolute precision required for the distribution 
         function. Relative precision via \code{pnvmix.reltol} is only used 
         when \code{pnvmix.abstol = NA}; in all other cases, absolute precision
         will be used. \code{pnvmix.abstol} defaults to 1e-3. 
         If \code{pnvmix.abstol = 0} and \code{pnvmix.reltol = 0}, the algorithm 
         will typically run until the total number of function evaluations 
         exceeds \code{fun.eval[2]} or until the total number of iterations exeeds
         \code{max.iter.rqmc}, whichever happens first. 
         If \eqn{n > 1} (so \code{upper} has more than 
         one row), the algorithm runs until the precision requirement is reached 
         for all \eqn{n} probability estimates.}
         \item{\code{mean.sqrt.mix}}{expectation of the square root \eqn{\sqrt(W)} 
         of the mixing variable \eqn{W}. If \code{NULL}, it will be estimated via 
         QMC; this is only needed for determining the reordering of the integration 
         bounds, so a rather crude approximation is fine.}
         \item{\code{precond}}{\code{\link{logical}} indicating whether 
         preconditioning is applied, that is, reordering of the integration 
         variables. If \code{TRUE}, integration limits as well as \code{scale} 
         are internally re-ordered in a way such that the overall variance of the 
         integrand is usually smaller than with the original ordering; this 
         usually leads smaller run-times.}
         \item{\code{cholesky.tol}}{non-negative numeric providing lower threshold
         for non-zero elements in the computation of the cholesky factor: If 
         calculated \eqn{C(i,i)^2 < | cholesky.tol * Scale(i,i)|}, the diagonal
         element (and all other elements in column \eqn{i}) of the cholesky factor
         \eqn{C} is set to zero, yielding a singular matrix. \code{cholesky.tol}
         defaults to 1e-9.}
      }
   } 
  \item{For \code{dnvmix()}:}{
      \describe{
         \item{\code{dnvmix.reltol}, \code{dnvmix.abstol}}{non-negative numeric 
         providing the relative/absolute precision for the *log-* density required. 
         Absolute precision via \code{dnvmix.abstol} is only used when 
         \code{dnvmix.reltol = NA}; in all other cases, relative precision will be 
         used. \code{dnvmix.reltol} defaults to 1e-2. 
         If \code{dnvmix.reltol=0} and \code{dnvmix.abstol=0}, the algorithm will 
         typically run until the total number of function evaluations exceeds 
         \code{fun.eval[2]} or until the total number of iterations exeeds
         \code{max.iter.rqmc}, whichever happens first.
         If \eqn{n > 1} (so \code{x} has more than one row), the algorithm runs until 
         the precision requirement is reached for all \eqn{n} log-density estimates.}
         \item{\code{dnvmix.doAdapt}}{\code{\link{logical}} indicating if an
         adaptive integration procedure shall be used; defaults to \code{TRUE}.}
         \item{\code{dnvmix.max.iter.rqmc.pilot}}{\code{\link{numeric}}, providing
         the maximum number of unstratified pilot runs the internal integration
         procedure performs. Defaults to 4.}
         \item{\code{dnvmix.tol.int.lower}, \code{dnvmix.order.lower}}{
         both \code{\link{numeric}} and nonnegative. RQMC integration 
         is only performed where the integrand is > than the maximum of 
         \code{dnvmix.tol.int.lower} and \eqn{10^{-c} g_{max}}, where \eqn{g_{max}}
         is the theoretical maximum of the integrand and \eqn{c} is the specified
         \code{dnvmix.order.lower}. Default to 1e-30 and 10, respectively.}
         \item{\code{dnvmix.tol.bisec}}{\code{\link{numeric}} \code{vector} of 
         length 3 specifying bisection tolerances in the adaptive RQMC algorithm.
         First/second/third element specify the tolerance on \eqn{u}, \eqn{W} and
         the log-integrand and default to 1e-16, 1e-1 and 1e-1, respectively.}
         \item{\code{dnvmix.max.iter.bisec}}{\code{\link{numeric}}, maximum number
         of iterations in the internal bisection procedure to find good cutting 
         points allowed.}
         \item{\code{dnvmix.tol.stratlength}}{\code{\link{numeric}}, nonnegative.
         If the stratum found by the adaptive integration method has length >
         \code{dnvmix.tol.stratlength} RQMC integration is used there; otherwise
         a crude approximation. Defaults to 1e-20.}
      }
   } 
   \item{For \code{fitnvmix()}:}{
      \describe{
         \item{\code{ECMEstep}}{\code{\link{logical}}, if \code{TRUE} (default), 
         ECME iteration is performed; 
         if \code{FALSE}, no ECME step is performed so that \code{fitnvmix()}
         performs between zero and two optimizations over \eqn{nu}, depending on
         \code{laststep.do.nu} and wether \code{nu.init} was provided.}
         \item{\code{ECMEstep.do.nu}}{\code{\link{logical}}, if \code{TRUE} 
         (default), the likelihood is maximed over \eqn{nu} in each ECME iteration; 
         if \code{FALSE}, this step is omitted.}
         \item{\code{laststep.do.nu}}{\code{\link{logical}}, if \code{TRUE} 
         another last maximization of the likelihood over \eqn{nu} is 
         performed using all  observations after the ECME iterations. 
         Only makes sense if either \code{ECMEstep.do.nu=FALSE} or
         if \code{size.subsample} is smaller than the number of observations.
         Defaults to \code{FALSE}.}
         \item{\code{resample}}{\code{\link{logical}}, if \code{TRUE}, a different
         subsample of \code{x} is taken in each optimization over \eqn{nu} in 
         the ECME iterations. Only relevant when \code{size.subsample} is 
         smaller than the number of observations. Defaults to \code{FALSE}.}
         \item{\code{ECME.maxiter}}{\code{\link{numeric}} positive, maximum number of 
         ECME iterations. Defaults to 25.}
         \item{\code{max.iter.locscaleupdate}}{\code{\link{numeric}} positive. Maximum
         number of loc-scale updates (while helding \eqn{nu} fixed) in each
         individual ECME iteration.}
         \item{\code{weights.reltol}}{\code{numeric} non-negative. Relative tolerance to
         estimate internal weights used to update \eqn{loc} and \eqn{scale} in 
         the ECME iterations. Defaults to 5e-2.}
         \item{\code{weights.interpol.reltol}}{\code{numeric} non-negative. Some weights
         can be obtained by interpolating previously calculated weights; if the
         maximal relative interpolation error is smaller than 
         \code{weights.interpol.reltol}, this is done. Defaults to 1e-2.}
         \item{\code{ECME.rel.conv.tol}}{\code{\link{numeric}(3)} vector specifying relative
         convergence tolerances for \code{loc}, \code{scale} and \code{nu} (in this
         order). Defaults to \code{c(5e-2, 5e-2, 1e-2)}.}
         \item{\code{control.optim}}{\code{\link{list}} of control parameters passed to 
         the underlying \code{optim} in the initial step as well as in the ECME
         iterations. See \code{?optim} for details; defaults to 
         \code{list(maxit=10)}.}
         \item{\code{control.optim.laststep}}{like \code{control.optim}; this control is
         passed to \code{optim} in the last-step. Only relevant when 
         \code{laststep.do.nu = TRUE}.}
         }
      }
   }
  Care should be taken when changing algorithm specific parameters, notably
  tolerances, as the accuracy of the result is heavily influenced by those.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2019),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \url{https://arxiv.org/abs/1911.03017}. 
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pnvmix}()}, \code{\link{fitnvmix}()},
  \code{\link{rnvmix}()}
}
\examples{
get.set.parameters() # obtain defaults
}
\keyword{distribution}