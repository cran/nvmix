\name{qnvmix}
\alias{qnvmix}
\title{Quantile Function of a univariate Normal Variance Mixture Distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t} for non-integer degrees of freedom).
}
\usage{
qnvmix(u, qmix, control = list(),
       verbose = TRUE, q.only = TRUE, stored.values = NULL, ...)
}
\arguments{
  \item{u}{vector of probabilities .}
  \item{qmix}{specification of the mixing variable \eqn{W}; see
  \code{\link{pnvmix}()} for details and examples.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see details below.} 
  \item{verbose}{\link{logical}, if \code{TRUE} a warning is printed if one of the 'abstol' is
  not reached.}
  \item{q.only}{\link{logical}. If \code{TRUE}, only the quantiles are returned; if FALSE, see 
  Section 'value' below.}
  \item{stored.values}{\link{matrix} with 3 columns of the form [x, F(x), logf(x)] where
  F and logf are the df and log-density of the distribution specified
  in 'qmix'. 
  If provided it will be used to determine starting values for 
  the internal newton proceudure. Only very basic checking is done.}
  \item{\dots}{additional arguments containing parameters of
    mixing distributions when \code{mix} is a \code{\link{character}}
    string.}
}
\value{
  If \code{q.only = TRUE} a vector of the same length as \code{u} with entries
  \eqn{q_i} where \eqn{q_i} satisfies \eqn{q_i = inf_x { F(x) \ge u_i}} 
  where \eqn{F(x)} the univariate df of the normal variance mixture specified 
  via \code{qmix}; 
  
  if \code{q.only = FALSE} a list of four:
  \describe{
  \item{\code{$q}:}{Vector of quantiles}
  \item{\code{$log.density}:}{Vector log-density values at \code{q}}
  \item{\code{$computed.values}:}{matrix with 3 columns [x, F(x), logf(x)];
  see details above}
  \item{\code{$newton.iterations}:}{Vector giving the number of Newton iterations 
  needed for \code{u[i]}}}
}
\details{
  This function uses a Newton procedure to estimate the quantile of the specified
  univariate normal variance mixture distribution. Internally, a randomized 
  quasi-Monte Carlo (RQMC) approach is used to estimate the distribution and (log)density
  function; the method is similar to the one in \link{pnvmix} and \link{dnvmix}.
  The result depends slightly on \code{.random.seed}.
  
  Internally, symmetry is used for \eqn{u \le 0.5}. Function values (i.e. cdf and
  log-density values) are stored and reused to get good starting values. These
  values are returned if \code{q.only = FALSE} and can be re-used by passing it to
  \code{qnvmix} via the argument \code{stored.values}; this can significantly
  reduce run-time.
  
  Accuracy and run-time depend on both the magnitude of \eqn{u} and on how heavy
  the tail of the underlying distributions is. 
  Numerical instabilities can occur for values of \eqn{u} close to 0 or 1, especially
  when the tail of the distribution is heavy. 
  
  If \code{q.only = TRUE} the log-density values of the underlying distribution
  evaluated at the estimated quantiles are returned as well: This can be useful
  for copula density evaluations where both quantities are needed. 
  
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
  \item{\code{abstol.cdf}}{abstol to estimate the df F(x) internally. See also ?pnvmix.}
  \item{abstol.logdensity}{abstol to estimate the log-density logf(x) internally. 
    See also ?dnvmix.}  
  \item{\code{abstol.newton}}{Convergence criterion for the internal Newton method.}  
  \item{\code{max.iter.newton}}{Maximum number of iterations for the internal Newton 
    method for *each* entry of \code{u}}
  \item{\code{max.iter.rqmc}}{\code{\link{numeric}}, providing the maximum number of 
    iterations allowed in the RQMC approach; the default is 15.}  
  \item{\code{CI.factor}}{multiplier of the Monte Carlo confidence interval
    bounds.  The algorithm runs until \code{CI.factor} times the estimated
    standard error is less than \code{abstol}. If \code{CI.factor = 3.3}
    (the default), one can expect the actual absolute error to be less
    than \code{abstol} in 99.9\% of the cases.}
  \item{\code{n0}}{Size of initial point-set used to internally approximate the df.}
  \item{\code{B}}{number of randomizations for obtaining an error estimate in the
    randomized quasi-Monte Carlo (RQMC) approach; the default is 8.}  
  }  
  
  Care should be taken when changing the algorithm-specific parameters,
  notably \code{method}, \code{precond}, \code{fun.eval[2]} and \code{B}.
  Error estimates will not be reliable for too small \code{B} and the
  performance of the algorithm depends heavily on the (quasi-)Monte
  Carlo point-set used.

  
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{rnvmix}()}, \code{\link{pnvmix}()}
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2019),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \url{https://arxiv.org/abs/1911.03017}. 
  
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
### Examples for qnvmix() ####################################################

## Evaluation points
u <- seq(from = 0.05, to = 0.95, by = 0.1)
set.seed(271) # for reproducibility

## Evaluate the t_{1.4} quantile function
df = 1.4
qmix. <- function(u) 1/qgamma(1-u, shape = df/2, rate = df/2)
## If qmix = "inverse.gamma", qt() is being called 
qt1  <- qnvmix(u, qmix = "inverse.gamma", df = df)
## Estimate quantiles (without using qt())
qt1. <- qnvmix(u, qmix = qmix., q.only = FALSE)
stopifnot(all.equal(qt1, qt1.$q, tolerance = 1e-3))
## Look at absolute error:
abs.error <- abs(qt1 - qt1.$q)
plot(u, abs.error, type = "l", xlab = "u", ylab = "qt(u)")
## Now do this again but provide qt1.$stored.values, in which case at most
## one Newton iteration will be needed:
qt2 <- qnvmix(u, qmix = qmix., stored.values = qt1.$computed.values, q.only = FALSE)
stopifnot(max(qt2$newton.iterations) <= 1)


## Evaluate quantile function where W~Exp(2)
rate = 2
qexp <- qnvmix(u, qmix = list("exp", rate = rate))
## Check: F( F^{-1}(u)) = u 
stopifnot(all.equal(pnvmix(as.matrix(qexp), qmix = list("exp", rate = rate)), u, 
                    tolerance = 5e-4, check.attributes = FALSE))
}
\keyword{distribution}