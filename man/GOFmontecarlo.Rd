% 2005-09-23, Alberto Viglione
%
\name{GOFmontecarlo}
\alias{GOFmontecarlo}
\alias{gofNORMtest}
\alias{gofEXPtest}
\alias{gofGUMBELtest}
\alias{gofGENLOGIStest}
\alias{gofGENPARtest}
\alias{gofGEVtest}
\alias{gofLOGNORMtest}
\alias{gofP3test}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Goodness of fit tests}
\description{
 Anderson-Darling goodness of fit tests for Regional Frequency Analysis: Monte-Carlo method.
}
\usage{
 gofNORMtest (x)
 gofEXPtest (x, Nsim=1000)
 gofGUMBELtest (x, Nsim=1000)
 gofGENLOGIStest (x, Nsim=1000)
 gofGENPARtest (x, Nsim=1000)
 gofGEVtest (x, Nsim=1000)
 gofLOGNORMtest (x, Nsim=1000)
 gofP3test (x, Nsim=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{data sample}
 \item{Nsim}{number of simulated samples from the hypothetical parent distribution}
}
\details{
An introduction, analogous to the following one, on the Anderson-Darling test is available on \url{http://en.wikipedia.org/wiki/Anderson-Darling_test}.

Given a sample \eqn{x_i \  (i=1,\ldots,m)}{xi (i=1,...,m)} of data extracted from a distribution \eqn{F_R(x)}{FR(x)}, the test is used to check the null hypothesis \eqn{H_0 : F_R(x) = F(x,\theta)}{H0 : FR(x) = F(x,\theta)}, where \eqn{F(x,\theta)} is the hypothetical distribution and \eqn{\theta} is an array of parameters estimated from the sample \eqn{x_i}{xi}.

The Anderson-Darling goodness of fit test measures the departure between the hypothetical distribution \eqn{F(x,\theta)} and the cumulative frequency function \eqn{F_m(x)}{Fm(x)} defined as:
\deqn{F_m(x) = 0 \ , \ x < x_{(1)}}{Fm(x)=0, x<x(1)}
\deqn{F_m(x) = i/m \ , \  x_{(i)} \leq x < x_{(i+1)}}{Fm(x)=i/m, x(i)<=x<x(i+1)}
\deqn{F_m(x) = 1 \ , \ x_{(m)} \leq x}{Fm(x)=1, x(m)<=x}
where \eqn{x_{(i)}}{x(i)} is the \eqn{i}-th element of the ordered sample (in increasing order).

The test statistic is:
\deqn{Q^2 = m \! \int_x \left[ F_m(x) - F(x,\theta) \right]^2 \Psi(x) \,dF(x)}{Q2 = m int[Fm(x) - F(x,\theta)]^2 \Psi(x) dF(x)}
where \eqn{\Psi(x)}, in the case of the Anderson-Darling test (Laio, 2004), is \eqn{\Psi(x) = [F(x,\theta) (1 - F(x,\theta))]^{-1}}.
In practice, the statistic is calculated as:
\deqn{A^2 = -m -\frac{1}{m} \sum_{i=1}^m \left\{ (2i-1)\ln[F(x_{(i)},\theta)] + (2m+1-2i)\ln[1 - F(x_{(i)},\theta)] \right\}}{A2 = -m -1/m sum\{(2i-1)\ln[F(x(i),\theta)] + (2m+1-2i)\ln[1 - F(x(i),\theta)]\}}

The statistic \eqn{A^2}{A2}, obtained in this way, may be confronted with the population of the \eqn{A^2}{A2}'s that one obtain if samples effectively belongs to the \eqn{F(x,\theta)} hypothetical distribution.
In the case of the test of normality, this distribution is defined (see Laio, 2004).
In other cases, e.g. the Pearson Type III case, can be derived with a Monte-Carlo procedure.
}
\value{
  \code{gofNORMtest} tests the goodness of fit of a normal (Gauss) distribution with the sample \code{x}.

  \code{gofEXPtest} tests the goodness of fit of a exponential distribution with the sample \code{x}.

  \code{gofGUMBELtest} tests the goodness of fit of a Gumbel (EV1) distribution with the sample \code{x}.

  \code{gofGENLOGIStest} tests the goodness of fit of a Generalized Logistic distribution with the sample \code{x}.

  \code{gofGENPARtest} tests the goodness of fit of a Generalized Pareto distribution with the sample \code{x}.

  \code{gofGEVtest} tests the goodness of fit of a Generalized Extreme Value distribution with the sample \code{x}.

  \code{gofLOGNORMtest} tests the goodness of fit of a 3 parameters Lognormal distribution with the sample \code{x}.

  \code{gofP3test} tests the goodness of fit of a Pearson type III (gamma) distribution with the sample \code{x}.

  They return the value \eqn{A_2}{A2} of the Anderson-Darling statistics and its non exceedence probability \eqn{P}.
  Note that \eqn{P} is the probability of obtaining the test statistic \eqn{A_2}{A2} lower than the one that was actually observed, assuming that the null hypothesis is true, i.e., \eqn{P} is one minus the p-value usually employed in statistical testing (see \url{http://en.wikipedia.org/wiki/P-value}).
  If \eqn{P(A_2)}{P(A2)} is, for example, greater than 0.90, the null hypothesis at significance level \eqn{\alpha=10\%} is rejected.
  %If \eqn{P} is, for example, 0.92, the sample shouldn't be considered extracted from the hypothetical parent distribution with significance level greater than 8%.
}
%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{traceWminim}}, \code{\link{roi}}, \code{\link{HOMTESTS}}.}

\examples{
x <- rnorm(30,10,1)
gofNORMtest(x)

x <- rand.gamma(50, 100, 15, 7)
gofP3test(x, Nsim=200)

x <- rand.GEV(50, 0.907, 0.169, 0.0304)
gofGEVtest(x, Nsim=200)

x <- rand.genlogis(50, 0.907, 0.169, 0.0304)
gofGENLOGIStest(x, Nsim=200)

x <- rand.genpar(50, 0.716, 0.418, 0.476)
gofGENPARtest(x, Nsim=200)

x <- rand.lognorm(50, 0.716, 0.418, 0.476)
gofLOGNORMtest(x, Nsim=200)

}
\keyword{htest}
