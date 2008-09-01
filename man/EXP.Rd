% 2005-09-23, Alberto Viglione
%
\name{EXP}
\alias{EXP}
\alias{f.exp}
\alias{F.exp}
\alias{invF.exp}
\alias{Lmom.exp}
\alias{par.exp}
\alias{rand.exp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Two parameter exponential distribution and L-moments}
\description{
\code{EXP} provides the link between L-moments of a sample and the two parameter
exponential distribution.
}
\usage{
f.exp (x, xi, alfa)
F.exp (x, xi, alfa)
invF.exp (F, xi, alfa)
Lmom.exp (xi, alfa)
par.exp (lambda1, lambda2)
rand.exp (numerosita, xi, alfa)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{xi}{vector of exp location parameters}
 \item{alfa}{vector of exp scale parameters}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
}
\details{
See \url{http://en.wikipedia.org/wiki/Exponential_distribution} for a brief introduction on the Exponential distribution.

\bold{Definition}

Parameters (2): \eqn{\xi} (lower endpoint of the distribution), \eqn{\alpha} (scale).

Range of \eqn{x}: \eqn{\xi \le x < \infty}.

Probability density function:
\deqn{f(x) = \alpha^{-1} \exp\{-(x-\xi)/\alpha\}}

Cumulative distribution function:
\deqn{F(x) = 1 - \exp\{-(x-\xi)/\alpha\}}

Quantile function:
\deqn{x(F) = \xi - \alpha \log(1-F)}


\bold{L-moments}

\deqn{\lambda_1 = \xi + \alpha}
\deqn{\lambda_2 = 1/2 \cdot \alpha}
\deqn{\tau_3 = 1/3}
\deqn{\tau_4 = 1/6}


\bold{Parameters}

If \eqn{\xi} is known, \eqn{\alpha} is given by \eqn{\alpha = \lambda_1 - \xi} and the L-moment, moment, and maximum-likelihood estimators are identical.
If \eqn{\xi} is unknown, the parameters are given by
\deqn{\alpha = 2 \lambda_2}
\deqn{\xi = \lambda_1 - \alpha}
For estimation based on a single sample these estimates are inefficient, but in regional frequency analysis they can give reasonable estimates of upper-tail quantiles.

\code{Lmom.exp} and \code{par.exp} accept input as vectors of equal length. In \code{f.exp}, \code{F.exp}, \code{invF.exp} and \code{rand.exp} parameters (\code{xi}, \code{alfa}) must be atomic.
}
\value{
 \code{f.exp} gives the density \eqn{f}, \code{F.exp} gives the distribution function \eqn{F}, \code{invFexp} gives
the quantile function \eqn{x}, \code{Lmom.exp} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}), \code{par.exp} gives the parameters (\code{xi}, \code{alfa}), and \code{rand.exp} generates random deviates.
}
%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{GENLOGIS}}, \code{\link{GENPAR}}, \code{\link{GEV}}, \code{\link{GUMBEL}}, \code{\link{KAPPA}}, \code{\link{LOGNORM}}, \code{\link{P3}}; \code{\link{DISTPLOTS}}, \code{\link{GOFmontecarlo}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.exp(ll[1],ll[2])
f.exp(1800,parameters$xi,parameters$alfa)
F.exp(1800,parameters$xi,parameters$alfa)
invF.exp(0.7870856,parameters$xi,parameters$alfa)
Lmom.exp(parameters$xi,parameters$alfa)
rand.exp(100,parameters$xi,parameters$alfa)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.exp(Rll[1],Rll[2])
Lmom.exp(parameters$xi,parameters$alfa)
}
\keyword{distribution}
