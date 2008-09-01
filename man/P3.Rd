% 2005-09-23, Alberto Viglione
%
\name{P3}
\alias{P3}
\alias{f.gamma}
\alias{F.gamma}
\alias{invF.gamma}
\alias{Lmom.gamma}
\alias{par.gamma}
\alias{rand.gamma}
\alias{mom2par.gamma}
\alias{par2mom.gamma}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Three parameters Pearson type III distribution and L-moments}
\description{
\code{P3} provides the link between L-moments of a sample and the three parameter
Pearson type III distribution.
}
\usage{
f.gamma (x, xi, beta, alfa)
F.gamma (x, xi, beta, alfa)
invF.gamma (F, xi, beta, alfa)
Lmom.gamma (xi, beta, alfa)
par.gamma (lambda1, lambda2, tau3)
rand.gamma (numerosita, xi, beta, alfa)
mom2par.gamma (mu, sigma, gamm)
par2mom.gamma (alfa, beta, xi)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector of quantiles}
 \item{mu}{vector of gamma mean}
 \item{sigma}{vector of gamma standard deviation}
 \item{gamm}{vector of gamma third moment}
 \item{F}{vector of probabilities}
 \item{lambda1}{vector of sample means}
 \item{lambda2}{vector of L-variances}
 \item{tau3}{vector of L-CA (or L-skewness)}
 \item{numerosita}{numeric value indicating the length of the vector to be generated}
 \item{alfa}{vector of gamma shape parameters}
 \item{beta}{vector of gamma scale parameters}
 \item{xi}{vector of gamma location parameters}
}
\details{
See \url{http://en.wikipedia.org/wiki/Pearson_distribution} for an introduction to the Pearson distribution, and \url{http://en.wikipedia.org/wiki/Gamma_distribution} for an introduction to the Gamma distribution (the Pearson type III distribution is, essentially, a Gamma distribution with 3 parameters).

\bold{Definition}

Parameters (3): \eqn{\mu} (location), \eqn{\sigma} (scale), \eqn{\gamma} (shape).

If \eqn{\gamma \ne 0}, let \eqn{\alpha=4/\gamma^2}, \eqn{\beta=\frac{1}{2}\sigma |\gamma|}, and \eqn{\xi= \mu - 2 \sigma/\gamma}.
If \eqn{\gamma > 0}, then the range of \eqn{x} is \eqn{\xi \le x < \infty} and
\deqn{f(x) = \frac{(x - \xi)^{\alpha - 1} e^{-(x-\xi)/\beta}}{\beta^{\alpha} \Gamma(\alpha)}}
\deqn{F(x) = G \left(\alpha, \frac{x-\xi}{\beta}\right)/ \Gamma(\alpha)}

If \eqn{\gamma=0}, then the distribution is Normal, the range of \eqn{x} is \eqn{-\infty < x < \infty} and
\deqn{f(x) = \phi \left(\frac{x-\mu}{\sigma}\right)}
\deqn{F(x) = \Phi \left(\frac{x-\mu}{\sigma}\right)}
where
\eqn{\phi(x)=(2\pi)^{-1/2}\exp(-x^2/2)} and
\eqn{\Phi(x)=\int_{-\infty}^x \phi(t)dt}.

If \eqn{\gamma < 0}, then the range of \eqn{x} is \eqn{-\infty < x \le \xi} and
\deqn{f(x) = \frac{(\xi - x)^{\alpha - 1} e^{-(\xi-x)/\beta}}{\beta^{\alpha} \Gamma(\alpha)}}
\deqn{F(x) = G \left(\alpha, \frac{\xi-x}{\beta}\right)/ \Gamma(\alpha)}

In each case, \eqn{x(F)} has no explicit analytical form.
Here \eqn{\Gamma} is the gamma function, defined as \deqn{\Gamma (x) = \int_0^{\infty} t^{x-1} e^{-t} dt}
and
\deqn{G(\alpha, x) = \int_0^x t^{\alpha-1} e^{-t} dt}
is the incomplete gamma function.

\eqn{\gamma=2} is the exponential distribution; \eqn{\gamma=0} is the Normal distribution; \eqn{\gamma=-2} is the reverse exponential distribution.

The parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\gamma} are the conventional moments of the distribution.



\bold{L-moments}

Assuming \eqn{\gamma>0}, L-moments are defined for \eqn{0<\alpha<\infty}.

\deqn{\lambda_1 = \xi + \alpha \beta}
\deqn{\lambda_2 = \pi^{-1/2} \beta \Gamma(\alpha + 1/2)/\Gamma(\alpha)}
\deqn{\tau_3 = 6 I_{1/3} (\alpha, 2 \alpha)-3}
where \eqn{I_x(p,q)} is the incomplete beta function ratio
\deqn{I_x(p,q) = \frac{\Gamma(p+q)}{\Gamma(p)\Gamma(q)} \int_0^x t^{p-1} (1-t)^{q-1} dt}

There is no simple expression for \eqn{\tau_4}.
Here we use the rational-funcion approximation given by Hosking and Wallis (1997, pp. 201-202).

The corresponding results for \eqn{\gamma <0} are obtained by changing the signs of \eqn{\lambda_1}, \eqn{\tau_3} and \eqn{\xi} wherever they occur above.



\bold{Parameters}

\eqn{alpha} is obtained with an approximation.
If \eqn{0<|\tau_3|<1/3}, let \eqn{z=3 \pi \tau_3^2} and use
\deqn{\alpha \approx \frac{1+0.2906 z}{z + 0.1882 z^2 + 0.0442 z^3}}
if \eqn{1/3<|\tau_3|<1}, let \eqn{z=1-|\tau_3|} and use
\deqn{\alpha \approx \frac{0.36067 z - 0.59567 z^2 + 0.25361 z^3}{1-2.78861 z + 2.56096 z^2 -0.77045 z^3}}

Given \eqn{\alpha}, then
\eqn{\gamma=2 \alpha^{-1/2} sign(\tau_3)},
\eqn{\sigma=\lambda_2 \pi^{1/2} \alpha^{1/2} \Gamma(\alpha)/\Gamma(\alpha+1/2)},
\eqn{\mu=\lambda_1}.



\code{Lmom.gamma} and \code{par.gamma} accept input as vectors of equal length. 
In \code{f.gamma}, \code{F.gamma}, \code{invF.gamma} and \code{rand.gamma} parameters (\code{mu}, \code{sigma}, \code{gamm}) must be atomic.
}
\value{
 \code{f.gamma} gives the density \eqn{f}, \code{F.gamma} gives the distribution function \eqn{F}, \code{invFgamma} gives
the quantile function \eqn{x}, \code{Lmom.gamma} gives the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}), \code{par.gamma} gives the parameters (\code{mu}, \code{sigma}, \code{gamm}), and \code{rand.gamma} generates random deviates.

 \code{mom2par.gamma} returns the parameters \eqn{\alpha}, \eqn{\beta} and \eqn{\xi}, given the parameters (moments) \eqn{\mu}, \eqn{\sigma}, \eqn{\gamma}.
}

%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{rnorm}}, \code{\link{runif}}, \code{\link{EXP}}, \code{\link{GENLOGIS}}, \code{\link{GENPAR}}, \code{\link{GEV}}, \code{\link{GUMBEL}}, \code{\link{KAPPA}}, \code{\link{LOGNORM}}; \code{\link{DISTPLOTS}}, \code{\link{GOFmontecarlo}}, \code{\link{Lmoments}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
fac <- factor(annualflows["cod"][,])
split(x,fac)

camp <- split(x,fac)$"45"
ll <- Lmoments(camp)
parameters <- par.gamma(ll[1],ll[2],ll[4])
f.gamma(1800,parameters$xi,parameters$beta,parameters$alfa)
F.gamma(1800,parameters$xi,parameters$beta,parameters$alfa)
invF.gamma(0.7511627,parameters$xi,parameters$beta,parameters$alfa)
Lmom.gamma(parameters$xi,parameters$beta,parameters$alfa)
rand.gamma(100,parameters$xi,parameters$beta,parameters$alfa)

Rll <- regionalLmoments(x,fac); Rll
parameters <- par.gamma(Rll[1],Rll[2],Rll[4])
Lmom.gamma(parameters$xi,parameters$beta,parameters$alfa)

moments <- par2mom.gamma(parameters$alfa,parameters$beta,parameters$xi); moments
mom2par.gamma(moments$mu,moments$sigma,moments$gamm)
}
\keyword{distribution}
