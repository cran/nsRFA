% 2005-09-23, Alberto Viglione
%
\name{GENPAR}
\alias{GENPAR}
\alias{f.genpar}
\alias{F.genpar}
\alias{invF.genpar}
\alias{Lmom.genpar}
\alias{par.genpar}
\alias{rand.genpar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Three parameter generalized Pareto distribution and L-moments}
\description{
\code{GENPAR} provides the link between L-moments of a sample and the three parameter
generalized Pareto distribution.
}
\usage{
f.genpar (x, xi, alfa, k)
F.genpar (x, xi, alfa, k)
invF.genpar (F, xi, alfa, k)
Lmom.genpar (xi, alfa, k)
par.genpar (lambda1, lambda2, tau3)
rand.genpar (numerosita, xi, alfa, k)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{xi}{vector of genpar location parameters}
 \item{alfa}{vector of genpar scale parameters}
 \item{k}{vector of genpar shape parameters}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{tau3}{vector of L-CA (or L-skewness)}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
}
\details{
See \url{http://en.wikipedia.org/wiki/Pareto_distribution} for an introduction to the Pareto distribution.
 
\bold{Definition}

Parameters (3): \eqn{\xi} (location), \eqn{\alpha} (scale), \eqn{k} (shape).

Range of \eqn{x}: \eqn{\xi < x \le \xi + \alpha / k} if \eqn{k>0};
 \eqn{\xi \le x < \infty} if \eqn{k \le 0}.

Probability density function:
\deqn{f(x) = \alpha^{-1} e^{-(1-k)y}}
where \eqn{y = -k^{-1}\log\{1 - k(x - \xi)/\alpha\}} if \eqn{k \ne 0},
\eqn{y = (x-\xi)/\alpha} if \eqn{k=0}.

Cumulative distribution function:
\deqn{F(x) = 1-e^{-y}}

Quantile function:
\eqn{x(F) = \xi + \alpha[1-(1-F)^k]/k} if \eqn{k \ne 0},
\eqn{x(F) = \xi - \alpha \log(1-F)} if \eqn{k=0}.

\eqn{k=0} is the exponential distribution; \eqn{k=1} is the uniform distribution on the interval \eqn{\xi < x \le \xi + \alpha}.

\bold{L-moments}

L-moments are defined for \eqn{k>-1}.

\deqn{\lambda_1 = \xi + \alpha/(1+k)]}
\deqn{\lambda_2 = \alpha/[(1+k)(2+k)]}
\deqn{\tau_3 = (1-k)/(3+k)}
\deqn{\tau_4 = (1-k)(2-k)/[(3+k)(4+k)]}

The relation between \eqn{\tau_3} and \eqn{\tau_4} is given by
\deqn{\tau_4 = \frac{\tau_3 (1 + 5 \tau_3)}{5+\tau_3}}


\bold{Parameters}

If \eqn{\xi} is known, \eqn{k=(\lambda_1 - \xi)/\lambda_2 - 2} and \eqn{\alpha=(1+k)(\lambda_1 - \xi)};
if \eqn{\xi} is unknown, \eqn{k=(1 - 3 \tau_3)/(1 + \tau_3)}, \eqn{\alpha=(1+k)(2+k)\lambda_2} and
\eqn{\xi=\lambda_1 - (2+k)\lambda_2}.

}
\value{
 \code{f.genpar} gives the density \eqn{f}, \code{F.genpar} gives the distribution function \eqn{F}, \code{invF.genpar} gives
the quantile function \eqn{x}, \code{Lmom.genpar} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}), \code{par.genpar}
gives the parameters (\code{xi}, \code{alfa}, \code{k}), and \code{rand.genpar} generates random deviates.
}
\references{Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

\note{ \code{Lmom.genpar} and \code{par.genpar} accept input as vectors of equal length. In \code{f.genpar}, \code{F.genpar}, \code{invF.genpar} and \code{rand.genpar} parameters (\code{xi}, \code{alfa}, \code{k}) must be atomic.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{EXP}}, \code{\link{GENLOGIS}}, \code{\link{GEV}}, \code{\link{GUMBEL}}, \code{\link{KAPPA}}, \code{\link{LOGNORM}}, \code{\link{P3}}; \code{\link{DISTPLOTS}}, \code{\link{GOFmontecarlo}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.genpar(ll[1],ll[2],ll[4])
f.genpar(1800,parameters$xi,parameters$alfa,parameters$k)
F.genpar(1800,parameters$xi,parameters$alfa,parameters$k)
invF.genpar(0.7161775,parameters$xi,parameters$alfa,parameters$k)
Lmom.genpar(parameters$xi,parameters$alfa,parameters$k)
rand.genpar(100,parameters$xi,parameters$alfa,parameters$k)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.genpar(Rll[1],Rll[2],Rll[4])
Lmom.genpar(parameters$xi,parameters$alfa,parameters$k)
}
\keyword{distribution}
