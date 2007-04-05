% 2005-09-23, Alberto Viglione
%
\name{KAPPA}
\alias{KAPPA}
\alias{f.kappa}
\alias{F.kappa}
\alias{invF.kappa}
\alias{Lmom.kappa}
\alias{par.kappa}
\alias{rand.kappa}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Four parameter kappa distribution and L-moments}
\description{
\code{KAPPA} provides the link between L-moments of a sample and the four parameter
kappa distribution.
}
\usage{
f.kappa (x, xi, alfa, k, h)
F.kappa (x, xi, alfa, k, h)
invF.kappa (F, xi, alfa, k, h)
Lmom.kappa (xi, alfa, k, h)
par.kappa (lambda1, lambda2, tau3, tau4)
rand.kappa (numerosita, xi, alfa, k, h)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{xi}{vector of kappa location parameters}
 \item{alfa}{vector of kappa scale parameters}
 \item{k}{vector of kappa third parameters}
 \item{h}{vector of kappa fourth parameters}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{tau3}{vector of L-CA (or L-skewness)}
 \item{tau4}{vector of L-kurtosis}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
}
\details{
\bold{Definition}

Parameters (4): \eqn{\xi} (location), \eqn{\alpha} (scale), \eqn{k}, \eqn{h}.

Range of \eqn{x}: upper bound is \eqn{\xi + \alpha/k} if \eqn{k>0}, \eqn{\infty}{inf} if \eqn{k \le 0};
lower bound is \eqn{\xi + \alpha(1-h^{-k})/k}{xi + \alpha(1-h^(-k))/k} if \eqn{h>0}, \eqn{\xi + \alpha/k} if \eqn{h \le 0} and \eqn{k<0} and \eqn{-\infty}{-inf} if \eqn{h \le 0} and \eqn{k \ge 0}


Probability density function:
\deqn{f(x)=\alpha^{-1} [1-k(x-\xi)/\alpha]^{1/k-1} [F(x)]^{1-h}}{%
      f(x) = \alpha^(-1) [1 - k (x - \xi)/\alpha]^(1/k-1) [F(x)]^(1-h)}

Cumulative distribution function:
\deqn{F(x)=\{1-h[1-k(x-\xi)/\alpha]^{1/k}\}^{1/h}}{%
      F(x) = \{1 - h[1 - k(x - \xi)/\alpha]^(1/k)\}^(1/h)}

Quantile function:
\deqn{x(F)= \xi + \frac{\alpha}{k} \left[ 1-\left( \frac{1-F^h}{h} \right)^k \right]}{%
      x(F) = \xi + \alpha/k [1 - ((1-F^h)/h)^k]}



\eqn{h=-1} is the generalized logistic distribution; 
\eqn{h=0} is the generalized eztreme value distribution;
\eqn{h=1} is the generalized Pareto distribution.


\bold{L-moments}

L-moments are defined for \eqn{h \ge 0} and \eqn{k>-1}, or if \eqn{h<0} and \eqn{-1<k<-1/h}.

\deqn{\lambda_1 = \xi + \alpha(1-g_1)/k}{\lambda1 = \xi + \alpha(1 - g1)/k}
\deqn{\lambda_2 = \alpha(g_1 - g_2)/k}{\lambda2 = \alpha(g1 - g2)/k}
\deqn{\tau_3 = (-g_1 + 3g_2 - 2g_3)/(g_1 - g_2)}{\tau3 = (-g1 + 3g2 - 2g3)/(g1 - g2)}
\deqn{\tau_4 = (-g_1 + 6g_2 - 10g_3 + 5g_4)/(g_1 - g_2)}{\tau4 = (-g1 + 6g2 - 10g3 + 5g4)/(g1 - g2)}
where
\eqn{g_r = \frac{r\Gamma(1+k)\Gamma(r/h)}{h^{1+k}\Gamma(1+k+r/h)}}{gr = (r \Gamma(1+k)\Gamma(r/h))/(h^(1+k) \Gamma(1+k+r/h))} if \eqn{h>0};
\eqn{g_r = \frac{r \Gamma(1+k)\Gamma(-k-r/h)}{(-h)^{1+k}\Gamma(1-r/h)}}{gr = (r\Gamma(1+k)\Gamma(-k-r/h))/((-h)^(1+k)\Gamma(1-r/h))} if \eqn{h<0};

Here \eqn{\Gamma} denote the gamma function
\deqn{\Gamma (x) = \int_0^{\infty} t^{x-1} e^{-t} dt}{\Gamma(x) = integral[from 0 to inf] t^(x-1) e^(-t) dt}



\bold{Parameters}

There are no simple expressions for the parameters in terms of the L-moments.
However they can be obtained with a numerical algorithm considering the formulations of \eqn{\tau_3}{\tau3} and \eqn{\tau_4}{\tau4} in terms of \eqn{k} and \eqn{h}.
Here we use the function \code{optim} to minimize \eqn{(t_3-\tau_3)^2 + (t_4-\tau_4)^2}{(t3-\tau3)^2 + (t4-\tau4)^2} where \eqn{t_3}{t3} and \eqn{t_4}{t4} are the sample L-moment ratios.
}
\value{
 \code{f.kappa} gives the density \eqn{f}, \code{F.kappa} gives the distribution function \eqn{F}, \code{invFkappa} gives
the quantile function \eqn{x}, \code{Lmom.kappa} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}), \code{par.kappa} gives the parameters (\code{xi}, \code{alfa}, \code{k}, \code{h}), and \code{rand.kappa} generates random deviates.
}
\references{Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

\note{\code{Lmom.kappa} and \code{par.kappa} accept input as vectors of equal length. In \code{f.kappa}, \code{F.kappa}, \code{invF.kappa} and \code{rand.kappa} parameters (\code{xi}, \code{alfa}, \code{k}, \code{h}) must be atomic.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{EXP}}, \code{\link{GENLOGIS}}, \code{\link{GENPAR}}, \code{\link{GEV}}, \code{\link{GUMBEL}}, \code{\link{LOGNORM}}, \code{\link{P3}}; \code{\link{optim}}, \code{\link{DISTPLOTS}}, \code{\link{GOFmontecarlo}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.kappa(ll[1],ll[2],ll[4],ll[5])
f.kappa(1800,parameters$xi,parameters$alfa,parameters$k,parameters$h)
F.kappa(1800,parameters$xi,parameters$alfa,parameters$k,parameters$h)
invF.kappa(0.771088,parameters$xi,parameters$alfa,parameters$k,parameters$h)
Lmom.kappa(parameters$xi,parameters$alfa,parameters$k,parameters$h)
rand.kappa(100,parameters$xi,parameters$alfa,parameters$k,parameters$h)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.kappa(Rll[1],Rll[2],Rll[4],Rll[5])
Lmom.kappa(parameters$xi,parameters$alfa,parameters$k,parameters$h)
}
\keyword{distribution}
