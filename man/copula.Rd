\name{copula}
\alias{dnvmixcop}
\alias{pnvmixcop}
\alias{rnvmixcop}
\title{Functionalities for Normal Variance Mixture Copulas}
\description{
  Evaluate the density / distribution function of normal variance mixture
  copulas (including Student \emph{t} and normal copula) and generate vectors
  of random variates from normal variance mixture copulas.
}
\usage{
dnvmixcop(u, qmix, scale = diag(d), factor = NULL, control = list(),
          verbose = FALSE, log = FALSE, ...)
pnvmixcop(u, qmix, scale = diag(d), control = list(),
          verbose = FALSE, ...)
rnvmixcop(n, qmix, scale = diag(2), factor = NULL,
          method = c("PRNG", "sobol", "ghalton"), skip = 0,
          control = list(), verbose = FALSE, ...)
}
\arguments{
  \item{u}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points.
    Have to be in \eqn{(0,1)}.}
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{qmix}{specification of the mixing variable \eqn{W}; see
  \code{\link{pnvmix}()} for details and examples.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)} (defaults to
    \eqn{d = 2});
    this equals the covariance matrix of a random vector following
    the specified normal variance mixture distribution divided by
    the expecation of the mixing variable \eqn{W} if and only if the
    former exists.
    Note that \code{scale} must be positive definite; sampling from
    singular normal variance mixtures can be achieved by providing
    \code{factor}.}
  \item{factor}{% more general here than for dnvmix() or pnvmix()
    \eqn{(d, k)}-\code{\link{matrix}} such that
    \code{factor \%*\% t(factor)} equals \code{scale}; the non-square
    case \eqn{k \neq d}{k != d} can be used to sample from singular
    normal variance mixtures. For \code{dnvmixcop()}, this has to be
    a square matrix. Note that this notation coincides with
    McNeil et al. (2015, Chapter 6). If not provided, \code{factor} is
    internally determined via \code{\link{chol}()} (and multiplied
    from the right to an \eqn{(n, k)}-matrix of independent
    standard normals to obtain a sample from a multivariate normal
    with zero mean vector and covariance matrix \code{scale}).}
  \item{method}{\code{\link{character}} string indicating the method
    to be used to obtain the sample. Available are:
    \describe{
      \item{\code{"PRNG"}:}{pseudo-random numbers,}
      \item{\code{"sobol"}:}{Sobol' sequence,}
      \item{\code{"ghalton"}:}{generalized Halton sequence.}
    }
    If \code{method = "PRNG"}, either \code{qmix} or \code{rmix} can be
    provided. If both are provided, \code{rmix} is used and \code{qmix}
    ignored. For the other two methods, sampling is done via inversion,
    hence \code{qmix} has to be provided and \code{rmix} is ignored.
  }
  \item{skip}{\code{\link{integer}} specifying the number of points
    to be skipped when \code{method = "sobol"}, see also example below.
  }
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see details below.}
  \item{verbose}{
    \code{\link{logical}}
    indicating whether a warning is given if the required precision
    \code{abstol} has not been reached.}
  \item{log}{\code{\link{logical}} indicating whether the logarithmic
    density is to be computed.}
  \item{\dots}{additional arguments (for example, parameters) passed to the
    underlying mixing distribution when \code{rmix} or \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  The values returned by \code{dnvmixcop()}, \code{rnvmixcop()} and
  \code{pnvmixcop()} are similar to the ones returned by their
  non-copula alternatives \code{\link{dnvmix}()}, \code{\link{rnvmix}()}
  and \code{\link{pnvmix}()}.
}
\details{
  Functionalities for normal variance mixture copulas provided here
  essentially call \code{\link{pnvmix}()}, \code{\link{dnvmix}()} and
  \code{\link{rnvmix}()} as well as \code{\link{qnvmix}()}, see their
  documentations for more details.

  We remark that computing normal variance mixtures is a challenging
  task; evaluating normal variance mixture copulas additionally requires
  the approximation of a univariate quantile function so that for large
  dimensions and sample sizes, these procedures can be fairly slow. As
  there are approximations on many levels, reported error estimates for
  the copula versions of \code{pnvmix}() and \code{dnvmix}() can be
  flawed.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2019),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \url{https://arxiv.org/abs/1911.03017}.

  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pnvmix}()}, \code{\link{qnvmix}()},
  \code{\link{rnvmix}()}
}
\examples{
## Generate a random correlation matrix in d dimensions
d <- 2 # dimension
set.seed(42) # for reproducibility
rho <- runif(1, min = -1, max = 1)
P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
diag(P) <- 1
## Generate two random evaluation points:
u <- matrix(runif(2*d), ncol = d)
## We illustrate using a t-copula
df = 2.1
## Define quantile function which is inverse-gamma here:
qmix. <- function(u) 1/qgamma(1-u, shape = df/2, rate = df/2)


### Example for dnvmixcop() ####################################################

## If qmix = "inverse.gamma", dnvmix() calls qt and dt:
d1 <- dnvmixcop(u, qmix = "inverse.gamma", scale = P, df = df)
## Use qmix. to force the algorithm to use a rqmc procedure:
d2 <- dnvmixcop(u, qmix = qmix., scale = P)
stopifnot(all.equal(d1, d2, tol = 1e-3, check.attributes = FALSE))


### Example for pnvmixcop() ####################################################

## Same logic as above:
p1 <- pnvmixcop(u, qmix = "inverse.gamma", scale = P, df = df)
p2 <- pnvmixcop(u, qmix = qmix., scale = P)
stopifnot(all.equal(p1, p2, tol = 1e-3, check.attributes = FALSE))


### Examples for rnvmixcop() ###################################################

## Draw random variates and compare
n <- 100
set.seed(1)
X  <- rnvmixcop(n, qmix = "inverse.gamma", df = df, scale = P) # with scale
set.seed(1)
X. <- rnvmixcop(n, qmix = "inverse.gamma", df = df, factor = t(chol(P))) # with factor
stopifnot(all.equal(X, X.))
}
\keyword{distribution}