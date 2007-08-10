% 2005-09-23, Alberto Viglione
%
\name{GENLOGIS}
\alias{GENLOGIS}
\alias{f.genlogis}
\alias{F.genlogis}
\alias{invF.genlogis}
\alias{Lmom.genlogis}
\alias{par.genlogis}
\alias{rand.genlogis}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Three parameter generalized logistic distribution and L-moments}
\description{
\code{GENLOGIS} provides the link between L-moments of a sample and the three parameter
generalized logistic distribution.
}
\usage{
f.genlogis (x, xi, alfa, k)
F.genlogis (x, xi, alfa, k)
invF.genlogis (F, xi, alfa, k)
Lmom.genlogis (xi, alfa, k)
par.genlogis (lambda1, lambda2, tau3)
rand.genlogis (numerosita, xi, alfa, k)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{xi}{vector of genlogis location parameters}
 \item{alfa}{vector of genlogis scale parameters}
 \item{k}{vector of genlogis shape parameters}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{tau3}{vector of L-CA (or L-skewness)}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
}
\details{
See \url{http://en.wikipedia.org/wiki/Logistic_distribution} for an introduction to the Logistic Distribution.

\bold{Definition}

Parameters (3): \eqn{\xi} (location), \eqn{\alpha} (scale), \eqn{k} (shape).

Range of \eqn{x}: \eqn{-\infty < x \le \xi + \alpha / k} if \eqn{k>0};
 \eqn{-\infty < x < \infty} if \eqn{k=0};
 \eqn{\xi + \alpha / k \le x < \infty} if \eqn{k<0}.

Probability density function:
\deqn{f(x) = \frac{\alpha^{-1} e^{-(1-k)y}}{(1+e^{-y})^2}}
where \eqn{y = -k^{-1}\log\{1 - k(x - \xi)/\alpha\}} if \eqn{k \ne 0},
\eqn{y = (x-\xi)/\alpha} if \eqn{k=0}.

Cumulative distribution function:
\deqn{F(x) = 1/(1+e^{-y})}

Quantile function:
\eqn{x(F) = \xi + \alpha[1-\{(1-F)/F\}^k]/k} if \eqn{k \ne 0},
\eqn{x(F) = \xi - \alpha \log\{(1-F)/F\}} if \eqn{k=0}.

\eqn{k=0} is the logistic distribution.


\bold{L-moments}

L-moments are defined for \eqn{-1<k<1}.

\deqn{\lambda_1 = \xi + \alpha[1/k - \pi / \sin (k \pi)]}
\deqn{\lambda_2 = \alpha k \pi / \sin (k \pi)}
\deqn{\tau_3 = -k}
\deqn{\tau_4 = (1+5 k^2)/6}


\bold{Parameters}

\eqn{k=-\tau_3}, \eqn{\alpha = \frac{\lambda_2 \sin (k \pi)}{k \pi}}, 
\eqn{\xi = \lambda_1 - \alpha (\frac{1}{k} - \frac{\pi}{\sin (k \pi)})}.
}
\value{
 \code{f.genlogis} gives the density \eqn{f}, \code{F.genlogis} gives the distribution function \eqn{F}, \code{invF.genlogis} gives the quantile function \eqn{x}, \code{Lmom.genlogis} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}), \code{par.genlogis} gives the parameters (\code{xi}, \code{alfa}, \code{k}), and \code{rand.genlogis} generates random deviates.
}
\references{Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

\note{\code{Lmom.genlogis} and \code{par.genlogis} accept input as vectors of equal length. In \code{f.genlogis}, \code{F.genlogis}, \code{invF.genlogis} and \code{rand.genlogis} parameters (\code{xi}, \code{alfa}, \code{k}) must be atomic.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{EXP}}, \code{\link{GENPAR}}, \code{\link{GEV}}, \code{\link{GUMBEL}}, \code{\link{KAPPA}}, \code{\link{LOGNORM}}, \code{\link{P3}}; \code{\link{DISTPLOTS}}, \code{\link{GOFmontecarlo}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.genlogis(ll[1],ll[2],ll[4])
f.genlogis(1800,parameters$xi,parameters$alfa,parameters$k)
F.genlogis(1800,parameters$xi,parameters$alfa,parameters$k)
invF.genlogis(0.7697433,parameters$xi,parameters$alfa,parameters$k)
Lmom.genlogis(parameters$xi,parameters$alfa,parameters$k)
rand.genlogis(100,parameters$xi,parameters$alfa,parameters$k)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.genlogis(Rll[1],Rll[2],Rll[4])
Lmom.genlogis(parameters$xi,parameters$alfa,parameters$k)
}
\keyword{distribution}
