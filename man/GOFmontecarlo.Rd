% 2005-09-23, Alberto Viglione
%
\name{GOFmontecarlo}
\alias{GOFmontecarlo}
\alias{gofNORMtest}
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
In other cases, e.g. the Pearson Type III case here, can be derived with a Monte-Carlo procedure.
}
\value{
  \code{gofNORMtest} tests the goodness of fit of a normal (Gauss) distribution with the sample \code{x}.

  \code{gofGENLOGIStest} tests the goodness of fit of a Generalized Logistic distribution with the sample \code{x}.

  \code{gofGENPARtest} tests the goodness of fit of a Generalized Pareto distribution with the sample \code{x}.

  \code{gofGEVtest} tests the goodness of fit of a Generalized Extreme Value distribution with the sample \code{x}.

  \code{gofLOGNORMtest} tests the goodness of fit of a 3 parameters Lognormal distribution with the sample \code{x}.

  \code{gofP3test} tests the goodness of fit of a Pearson type III (gamma) distribution with the sample \code{x}.

  They return the value \eqn{A_2}{A2} of the Anderson-Darling statistics and its probability \eqn{P}.
  If \eqn{P(A_2)}{P(A2)} is, for example, greater than 0.90, the test is not passed at level \eqn{\alpha=10\%}.
  %If \eqn{P} is, for example, 0.92, the sample shouldn't be considered extracted from the hypothetical parent distribution with significance level greater than 8%.
}
\references{
D'Agostino R., Stephens M. (1986) Goodness-of-Fit Techniques, chapter Tests based on EDF statistics. Marcel Dekker, New York.

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.

Laio, F., Cramer-von Mises and Anderson-Darling goodness of fit tests for extreme value distributions with unknown parameters, Water Resour. Res., 40, W09308, doi:10.1029/2004WR003204.

Viglione A., Claps P., Laio F. (2006) Utilizzo di criteri di prossimit\`a nell'analisi regionale del deflusso annuo, XXX Convegno di Idraulica e Costruzioni Idrauliche - IDRA 2006, Roma, 10-15 Settembre 2006.

Viglione A. (2007) Metodi statistici non-supervised per la stima di grandezze idrologiche in siti non strumentati, PhD thesis
, In press.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

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
