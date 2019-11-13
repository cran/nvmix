\name{qqplot.maha}
\alias{qqplot.maha}
\title{QQ Plot for Multivariate Normal Variance Mixtures}
\description{
  Visual goodness-of-fit test for multivariate normal variance mixtures: Plotting
  squared mahalanobis distance against their theoretical quantiles. 
}
\usage{
qqplot.maha(x, qmix, loc, scale, plot.diag = TRUE, verbose = TRUE, 
            control = list(), ...)
}
\arguments{
  \item{x}{\eqn{(n, d)}-data \code{\link{matrix}}.}
  \item{qmix}{see \code{?\link{pnvmix}()}}
  \item{loc}{see \code{?\link{pnvmix}()}}
  \item{scale}{see \code{?\link{pnvmix}()}}
  \item{plot.diag}{\code{\link{logical}} indicating if the diagonal \eqn{y=x}
  shall be included in the plot.}
  \item{verbose}{see \code{?\link{pnvmix}()}}
  \item{control}{\code{\link{list}} specifying algorithm specific parameters; see \code{?\link{get.set.parameters}()}.}
  \item{\dots}{additional arguments (for example, parameters) passed to the
    underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{qqolot.maha()} (invisibly) returns a list of two: \code{maha2}, the
  sorted mahalanobis distances of the data in \code{x} with respect to \code{loc}
  and \code{scale} and \code{q}, the theoretical quantiles evaluated at 
  \code{ppoints(n)} where \code{n=nrow(x)}.
}
\details{
  If \eqn{X} follows a multivariate normal variance mixture, the distribution of
  the mahalanobis distance
  \eqn{D^2 = (X-\mu)^T \Sigma^{-1} (X-\mu)}
  is a gamma mixture whose distribution can be approximated.
  The function \code{qqplot.maha()} plots the empirical mahalanobis distances
  from the data in \code{x} (with \eqn{\mu=}\code{loc} and 
  \eqn{Sigma=}\code{scale}) versus their theoretical quantiles which are 
  internally estimated via the function \code{qgammamix()}.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2019),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \url{https://arxiv.org/abs/1911.03017}. 
  
  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.

  Genz, A. and Bretz, F. (1999).
  Numerical computation of multivariate t-probabilities with
  application to power calculation of multiple contrasts.
  \emph{Journal of Statistical Computation and Simulation} 63(4),
  103--117.

  Genz, A. and Bretz, F. (2002).
  Comparison of methods for the computation of multivariate \emph{t}
  probabilities.
  \emph{Journal of Computational and Graphical Statistics} 11(4),
  950--971.
}
\seealso{
  \code{\link{fitnvmix}()}, \code{\link{rnvmix}()}, \code{\link{pnvmix}()}
}
\examples{
set.seed(1) # for reproducibility
n     <- 21
d     <- 2
nu    <- 4 # degree-of-freedom parameter
loc   <- rep(0, d)
scale <- diag(d)
## Define 'qmix' as the quantile function of an IG(nu/2, nu/2) distribution
qmix <- function(u, nu){
    1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
}
## Sample data using 'rnvmix()':
x <- rnvmix(n, qmix = qmix, nu = nu, loc = loc, scale = scale)
## Call 'fitvnmix()' with 'qmix' as a function (so all densities/weights are estimated)
(MyFit <- fitnvmix(x, qmix = qmix, mix.param.bounds = c(0.5, 9)))
## QQ Plot of empirical quantiles vs estimated, theoretical quantiles:
qqplot.maha(x, qmix = qmix, loc = MyFit$loc, scale = MyFit$scale, nu = MyFit$nu)
}
\keyword{distribution}