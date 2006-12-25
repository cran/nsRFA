% 2005-09-23, Alberto Viglione
%
\name{GEV}
\alias{GEV}
\alias{f.GEV}
\alias{F.GEV}
\alias{invF.GEV}
\alias{Lmom.GEV}
\alias{par.GEV}
\alias{rand.GEV}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Three parameter generalized extreme value distribution and L-moments}
\description{
\code{GEV} provides the link between L-moments of a sample and the three parameter
generalized extreme value distribution.
}
\usage{
f.GEV (x, xi, alfa, k)
F.GEV (x, xi, alfa, k)
invF.GEV (F, xi, alfa, k)
Lmom.GEV (xi, alfa, k)
par.GEV (lambda1, lambda2, tau3)
rand.GEV (numerosita, xi, alfa, k)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{xi}{vector of GEV location parameters}
 \item{alfa}{vector of GEV scale parameters}
 \item{k}{vector of GEV shape parameters}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{tau3}{vector of L-CA (or L-skewness)}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
}
\details{
\bold{Definition}

Parameters (3): \eqn{\xi} (location), \eqn{\alpha} (scale), \eqn{k} (shape).

Range of \eqn{x}: \eqn{-\infty < x \le \xi + \alpha / k} if \eqn{k>0};
 \eqn{-\infty < x < \infty} if \eqn{k=0};
 \eqn{\xi + \alpha / k \le x < \infty} if \eqn{k<0}.

Probability density function:
\deqn{f(x) = \alpha^{-1} e^{-(1-k)y - e^{-y}}}
where \eqn{y = -k^{-1}\log\{1 - k(x - \xi)/\alpha\}} if \eqn{k \ne 0},
\eqn{y = (x-\xi)/\alpha} if \eqn{k=0}.

Cumulative distribution function:
\deqn{F(x) = e^{-e^{-y}}}

Quantile function:
\eqn{x(F) = \xi + \alpha[1-(-\log F)^k]/k} if \eqn{k \ne 0},
\eqn{x(F) = \xi - \alpha \log(-\log F)} if \eqn{k=0}.

\eqn{k=0} is the Gumbel distribution; \eqn{k=1} is the reverse exponential distribution.


\bold{L-moments}

L-moments are defined for \eqn{k>-1}.

\deqn{\lambda_1 = \xi + \alpha[1 - \Gamma (1+k)]/k}
\deqn{\lambda_2 = \alpha (1-2^{-k}) \Gamma (1+k)]/k}
\deqn{\tau_3 = 2(1-3^{-k})/(1-2^{-k})-3}
\deqn{\tau_4 = [5(1-4^{-k})-10(1-3^{-k})+6(1-2^{-k})]/(1-2^{-k})}

Here \eqn{\Gamma} denote the gamma function
\deqn{\Gamma (x) = \int_0^{\infty} t^{x-1} e^{-t} dt}


\bold{Parameters}

To estimate \eqn{k}, no explicit solution is possible, but the following approximation has accurancy better than \eqn{9 \times 10^{-4}} for \eqn{-0.5 \le \tau_3 \le 0.5}:
\deqn{k \approx 7.8590 c + 2.9554 c^2}
where
\deqn{c = \frac{2}{3+\tau_3} - \frac{\log 2}{\log 3}}
The other parameters are then given by
\deqn{\alpha = \frac{\lambda_2 k}{(1-2^{-k})\Gamma(1+k)}}
\deqn{\xi = \lambda_1 - \alpha[1 - \Gamma(1+k)]/k}

}
\value{
 \code{f.GEV} gives the density \eqn{f}, \code{F.GEV} gives the distribution function \eqn{F}, \code{invF.GEV} gives
the quantile function \eqn{x}, \code{Lmom.GEV} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}), \code{par.GEV} gives the parameters (\code{xi}, \code{alfa}, \code{k}), and \code{rand.GEV} generates random deviates.
}
\references{Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

\note{ \code{Lmom.GEV} and \code{par.GEV} accept input as vectors of equal length. In \code{f.GEV}, \code{F.GEV}, \code{invF.GEV} and \code{rand.GEV} parameters (\code{xi}, \code{alfa}, \code{k}) must be atomic.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{EXP}}, \code{\link{GENLOGIS}}, \code{\link{GENPAR}}, \code{\link{GUMBEL}}, \code{\link{KAPPA}}, \code{\link{LOGNORM}}, \code{\link{P3}}; \code{\link{DISTPLOTS}}, \code{\link{GOFTESTS}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.GEV(ll[1],ll[2],ll[4])
f.GEV(1800,parameters$xi,parameters$alfa,parameters$k)
F.GEV(1800,parameters$xi,parameters$alfa,parameters$k)
invF.GEV(0.7518357,parameters$xi,parameters$alfa,parameters$k)
Lmom.GEV(parameters$xi,parameters$alfa,parameters$k)
rand.GEV(100,parameters$xi,parameters$alfa,parameters$k)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.GEV(Rll[1],Rll[2],Rll[4])
Lmom.GEV(parameters$xi,parameters$alfa,parameters$k)
}
\keyword{distribution}
