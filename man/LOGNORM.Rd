% 2005-09-23, Alberto Viglione
%
\name{LOGNORM}
\alias{LOGNORM}
\alias{f.lognorm}
\alias{F.lognorm}
\alias{invF.lognorm}
\alias{Lmom.lognorm}
\alias{par.lognorm}
\alias{rand.lognorm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Three parameter lognormal distribution and L-moments}
\description{
\code{LOGNORM} provides the link between L-moments of a sample and the three parameter
log-normal distribution.
}
\usage{
f.lognorm (x, xi, alfa, k)
F.lognorm (x, xi, alfa, k)
invF.lognorm (F, xi, alfa, k)
Lmom.lognorm (xi, alfa, k)
par.lognorm (lambda1, lambda2, tau3)
rand.lognorm (numerosita, xi, alfa, k)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{xi}{vector of lognorm location parameters}
 \item{alfa}{vector of lognorm scale parameters}
 \item{k}{vector of lognorm shape parameters}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{tau3}{vector of L-CA (or L-skewness)}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
}
\details{
See \url{http://en.wikipedia.org/wiki/Log-normal_distribution} for an introduction to the lognormal distribution.

\bold{Definition}

Parameters (3): \eqn{\xi} (location), \eqn{\alpha} (scale), \eqn{k} (shape).

Range of \eqn{x}: \eqn{-\infty < x \le \xi + \alpha / k} if \eqn{k>0};
 \eqn{-\infty < x < \infty} if \eqn{k=0};
 \eqn{\xi + \alpha / k \le x < \infty} if \eqn{k<0}.

Probability density function:
\deqn{f(x) = \frac{e^{ky-y^2/2}}{\alpha \sqrt{2\pi}}}
where \eqn{y = -k^{-1}\log\{1 - k(x - \xi)/\alpha\}} if \eqn{k \ne 0},
\eqn{y = (x-\xi)/\alpha} if \eqn{k=0}.

Cumulative distribution function:
\deqn{F(x) = \Phi(x)}
where
\eqn{\Phi(x)=\int_{-\infty}^x \phi(t)dt}.

Quantile function:
\eqn{x(F)} has no explicit analytical form.

\eqn{k=0} is the Normal distribution with parameters \eqn{\xi} and \eqn{alpha}.


\bold{L-moments}

L-moments are defined for all values of \eqn{k}.

\deqn{\lambda_1 = \xi + \alpha(1 - e^{k^2/2})/k}
\deqn{\lambda_2 = \alpha/k e^{k^2/2} [1 - 2 \Phi(-k/\sqrt{2})]}

There are no simple expressions for the L-moment ratios \eqn{\tau_r} with \eqn{r \ge 3}.
Here we use the rational-function approximation given in Hosking and Wallis (1997, p. 199).


\bold{Parameters}

The shape parameter \eqn{k} is a function of \eqn{\tau_3} alone.
No explicit solution is possible.
Here we use the approximation given in Hosking and Wallis (1997, p. 199).

Given \eqn{k}, the other parameters are given by
\deqn{\alpha = \frac{\lambda_2 k e^{-k^2/2}}{1-2 \Phi(-k/\sqrt{2})}}
\deqn{\xi = \lambda_1 - \frac{\alpha}{k} (1 - e^{k^2/2})}


\code{Lmom.lognorm} and \code{par.lognorm} accept input as vectors of equal length. In \code{f.lognorm}, \code{F.lognorm}, \code{invF.lognorm} and \code{rand.lognorm} parameters (\code{xi}, \code{alfa}, \code{k}) must be atomic.
}
\value{
 \code{f.lognorm} gives the density \eqn{f}, \code{F.lognorm} gives the distribution function \eqn{F}, \code{invFlognorm} gives the quantile function \eqn{x}, \code{Lmom.lognorm} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}), \code{par.lognorm} gives the parameters (\code{xi}, \code{alfa}, \code{k}), and \code{rand.lognorm} generates random deviates.
}

%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{EXP}}, \code{\link{GENLOGIS}}, \code{\link{GENPAR}}, \code{\link{GEV}}, \code{\link{GUMBEL}}, \code{\link{KAPPA}}, \code{\link{P3}}; \code{\link{DISTPLOTS}}, \code{\link{GOFmontecarlo}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.lognorm(ll[1],ll[2],ll[4])
f.lognorm(1800,parameters$xi,parameters$alfa,parameters$k)
F.lognorm(1800,parameters$xi,parameters$alfa,parameters$k)
invF.lognorm(0.7529877,parameters$xi,parameters$alfa,parameters$k)
Lmom.lognorm(parameters$xi,parameters$alfa,parameters$k)
rand.lognorm(100,parameters$xi,parameters$alfa,parameters$k)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.lognorm(Rll[1],Rll[2],Rll[4])
Lmom.lognorm(parameters$xi,parameters$alfa,parameters$k)
}
\keyword{distribution}
