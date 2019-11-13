\name{gammamix}
\alias{dgammamix}
\alias{pgammamix}
\alias{qgammamix}
\alias{rgammamix}
\title{Functionalities for Gamma Scale Mixture Models }
\description{
  Evaluating density-, distribution- and quantile-function of Gamma scale
  mixtures as well as random variate generation.
}
\usage{
dgammamix(x, qmix, d, control = list(), verbose = TRUE, log = FALSE, ...)
pgammamix(x, qmix, d, lower.tail = TRUE, control = list(), verbose = TRUE, ...)
qgammamix(u, qmix, d, control = list(), verbose = TRUE, q.only = TRUE,
          stored.values = NULL, ...)
rgammamix(n, rmix = NULL, qmix = NULL, d, method = c("PRNG", "sobol", "ghalton"),
          skip = 0, ...)
}
\arguments{
  \item{x}{\eqn{n}-\code{\link{vector}} of evaluation points.}
  \item{u}{\eqn{n}-\code{\link{vector}} of probabilities.}
  \item{qmix}{see \code{?\link{pnvmix}()}}
  \item{rmix}{see \code{?\link{rnvmix}()}}
  \item{d}{dimension of the underlying normal variance mixture, see also
           details below.}
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{lower.tail}{\code{\link{logical}}; if \code{TRUE} (default), probabilities are
  \eqn{P(X<= x)}, otherwise \eqn{P(X > x)}.}
  \item{log}{\code{\link{logical}} indicating whether the log-density shall be returned.}
  \item{q.only}{see \code{?\link{qnvmix}()}}
  \item{stored.values}{see \code{?\link{qnvmix}()}}
  \item{method}{see \code{?\link{rnvmix}()}}
  \item{skip}{see \code{?\link{rnvmix}()}}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{?\link{get.set.parameters}()}}   
  \item{verbose}{
    \code{\link{logical}}
    indicating whether a warning is given if the required precision
    has not been reached.}
  \item{\dots}{additional arguments (for example, parameters) passed to the
    underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{pgammamix()} and \code{dgammamix()} return
  a \code{\link{numeric}} \eqn{n}-vector with the computed probabilities/densities
  and corresponding attributes \code{"error"} (error estimates of the
  RQMC estimator) and \code{"numiter"} (number of iterations).
  
  If \code{q.only = TRUE}, \code{qgammamix()} a vector of the same length as
  \code{u} with entries \eqn{q_i} where \eqn{q_i} satisfies 
  \eqn{q_i = inf_x { F(x) >= u_i}} where \eqn{F(x)} the df of the Gamma mixture 
  specified via qmix; if \code{q.only = FALSE}, see \code{\link{qnvmix}}.
  
  \code{rgammamix()} returns a \eqn{n}-\code{\link{vector}} containing \eqn{n}
  samples of the specified (via mix) Gamma mixture.
}
\details{
  We define a Gamma mixture as a random variable \eqn{Dsq} satisfying, in distribution,
  \eqn{Dsq = W*Gamma(d/2, 2)} where \eqn{W} is specified via \code{qmix}. 
  If \eqn{X} follows a \eqn{d-}dimensional normal variance mixture, the 
  squared mahalanobis distance \eqn{(X-\mu)^T Sigma^{-1}(X-\mu)} has the same
  distribution as \eqn{Dsq}. 
  
  The functions presented here are similar to the corresponding functions
  for normal variance mixtures (\code{d/p/q/rnvmix()}),
  details can be found in the corresponding help-files there.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2019),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \url{https://arxiv.org/abs/1911.03017}. 
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pnvmix}()}, \code{\link{qnvmix}()},
  \code{\link{rnvmix}()}, \code{\link{get.set.parameters}()},
  \code{\link{qqplot.maha}()}, \code{\link{fitnvmix}()}
}
\examples{
set.seed(271) # for reproducibility
## Specify inverse-gamma mixture => results in F(d, nu) dist'n,
## handled correctly when 'qmix = "inverse.gamma"' is specified
qmix <- function(u, nu) 1/qgamma(1 - u, shape = nu/2, rate = nu/2)

## Example for rgammamix()
n  <- 50
nu <- 3
d  <- 5
x  <- rgammamix(n, qmix = qmix, d = d, nu = nu)

## Evaluate distribution function at 'x'
p.true <- pgammamix(x, qmix = "inverse.gamma", d = d, df = nu) 
p.est  <- pgammamix(x, qmix = qmix, d = d, nu = nu)
stopifnot(all.equal(p.true, p.est, tol = 5e-4, check.attributes = FALSE))

## Evaluate density function at 'x'
d.true <- dgammamix(x, qmix = "inverse.gamma", d = d, df = nu)
d.est  <- dgammamix(x, qmix = qmix, d = d, nu = nu)
stopifnot(all.equal(p.true, p.est, tol = 5e-4, check.attributes = FALSE))

## Evaluate quantile function:
u <- seq(from = 0.05, to = 0.95, by = 0.05)
q.true <- qgammamix(u, qmix = "inverse.gamma", d = d, df = nu)
q.est  <- qgammamix(u, qmix = qmix, d = d, nu = nu)
stopifnot(all.equal(q.true, q.est, tol = 5e-4, check.attributes = FALSE))
}
\keyword{distribution}