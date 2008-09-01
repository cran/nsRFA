% 2005-09-23, Alberto Viglione
%
\name{GUMBEL}
\alias{GUMBEL}
\alias{f.gumb}
\alias{F.gumb}
\alias{invF.gumb}
\alias{Lmom.gumb}
\alias{par.gumb}
\alias{rand.gumb}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Two parameter Gumbel distribution and L-moments}
\description{
\code{GUMBEL} provides the link between L-moments of a sample and the two parameter
Gumbel distribution.
}
\usage{
f.gumb (x, xi, alfa)
F.gumb (x, xi, alfa)
invF.gumb (F, xi, alfa)
Lmom.gumb (xi, alfa)
par.gumb (lambda1, lambda2)
rand.gumb (numerosita, xi, alfa)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{xi}{vector of gumb location parameters}
 \item{alfa}{vector of gumb scale parameters}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
}
\details{
See \url{http://en.wikipedia.org/wiki/Fisher-Tippett_distribution} for an introduction to the Gumbel distribution.
 
\bold{Definition}

Parameters (2): \eqn{\xi} (location), \eqn{\alpha} (scale).

Range of \eqn{x}: \eqn{-\infty < x < \infty}.

Probability density function:
\deqn{f(x) = \alpha^{-1} \exp[-(x-\xi)/\alpha] \exp\{- \exp[-(x-\xi)/\alpha]\}}

Cumulative distribution function:
\deqn{F(x) = \exp[-\exp(-(x-\xi)/\alpha)]}

Quantile function:
\eqn{x(F) = \xi - \alpha \log(-\log F)}.


\bold{L-moments}

\deqn{\lambda_1 = \xi + \alpha \gamma}
\deqn{\lambda_2 = \alpha \log 2}
\deqn{\tau_3 = 0.1699 = \log(9/8)/ \log 2}
\deqn{\tau_4 = 0.1504 = (16 \log 2 - 10 \log 3)/ \log 2}

Here \eqn{\gamma} is Euler's constant, 0.5772...


\bold{Parameters}

\deqn{\alpha=\lambda_2 / \log 2}
\deqn{\xi = \lambda_1 - \gamma \alpha}



\code{Lmom.gumb} and \code{par.gumb} accept input as vectors of equal length. In \code{f.gumb}, \code{F.gumb}, \code{invF.gumb} and \code{rand.gumb} parameters (\code{xi}, \code{alfa}) must be atomic.
}
\value{
 \code{f.gumb} gives the density \eqn{f}, \code{F.gumb} gives the distribution function \eqn{F}, \code{invF.gumb} gives
the quantile function \eqn{x}, \code{Lmom.gumb} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4})), \code{par.gumb} gives the parameters (\code{xi}, \code{alfa}), and \code{rand.gumb} generates random deviates.
}

%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{EXP}}, \code{\link{GENLOGIS}}, \code{\link{GENPAR}}, \code{\link{GEV}}, \code{\link{KAPPA}}, \code{\link{LOGNORM}}, \code{\link{P3}}; \code{\link{DISTPLOTS}}, \code{\link{GOFmontecarlo}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows[1:10,]
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.gumb(ll[1],ll[2])
f.gumb(1800,parameters$xi,parameters$alfa)
F.gumb(1800,parameters$xi,parameters$alfa)
invF.gumb(0.7686843,parameters$xi,parameters$alfa)
Lmom.gumb(parameters$xi,parameters$alfa)
rand.gumb(100,parameters$xi,parameters$alfa)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.gumb(Rll[1],Rll[2])
Lmom.gumb(parameters$xi,parameters$alfa)
}
\keyword{distribution}
