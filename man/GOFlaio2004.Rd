% 2005-09-23, Alberto Viglione
%
\name{GOFlaio2004}
\alias{GOFlaio2004}
\alias{A2_GOFlaio}
\alias{A2}
\alias{W2}
\alias{Fx}
%\alias{fw2}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Goodness of fit tests}
\description{
 Anderson-Darling goodness of fit tests for extreme-value distributions, from Laio (2004).
}
\usage{
 A2_GOFlaio (x, dist="NORM")
 A2 (F)
 W2 (F)
 Fx (x, T, dist="NORM")
 %fw2 (w)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{data sample}
 \item{dist}{distribution: normal \code{"NORM"}, log-normal \code{"LN"}, Gumbel \code{"EV1"}, Frechet \code{"EV2"}, Generalized Extreme Value \code{"GEV"}, Pearson type III \code{"GAM"}, log-Pearson type III \code{"LP3"}}
 \item{F}{cumulative distribution function}
 \item{T}{parameters (position, scale, shape, ...)}
 %\item{w}{Transformed test statistic (Laio, 2004)}
}
\details{An introduction on the Anderson-Darling test is available on \url{http://en.wikipedia.org/wiki/Anderson-Darling_test} and in the \code{GOFmontecarlo} help page.
 The oroginal paper of Laio (2004) is available on his web site.
}
\value{
  \code{A2_GOFlaio} tests the goodness of fit of a distribution with the sample \code{x}; it return the value \eqn{A_2} of the Anderson-Darling statistics and its probability \eqn{P(A_2)}{P(A2)}. 
  If \eqn{P(A_2)}{P(A2)} is, for example, greater than 0.90, the test is not passed at level \eqn{\alpha=10\%}.

  \code{A2} is the Anderson-Darling test statistic; it is used by \code{A2_GOFlaio}.

  \code{W2} is the Cramer-von Mises test statistic.

  %\code{fw2} is the approximation of the probability distribution of \code{w} (first 2 terms) when \eqn{H_0} is true (Anderson-Darling, 1952); it is used by \code{A2_GOFlaio}.

 \code{Fx} gives the empirical distribution function \code{F} of a sample \code{x}, extracted from the distribution \code{dist} with parameters \code{T}.

}
\references{
Laio, F., Cramer-von Mises and Anderson-Darling goodness of fit tests for extreme value distributions with unknown parameters, Water Resour. Res., 40, W09308, doi:10.1029/2004WR003204.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{GOFmontecarlo}}, \code{\link{MLlaio2004}}.}

\examples{
sm <- sample_generator(100, c(0,1), dist="EV1")
ML_estimation (sm, dist="GEV")
Fx (sm, ML_estimation (sm, dist="GEV"), dist="GEV")
A2(sort(Fx(sm, ML_estimation (sm, dist="GEV"), dist="GEV")))
A2_GOFlaio(sm, dist="GEV")

ML_estimation (sm, dist="GAM")
Fx(sm, ML_estimation (sm, dist="GAM"), dist="GAM")
A2(sort(Fx(sm, ML_estimation (sm, dist="GAM"), dist="GAM")))
A2_GOFlaio(sm, dist="GAM")
}
\keyword{htest}
