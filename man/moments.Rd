% 2005-09-23, Alberto Viglione
%
\name{moments}
\alias{moments}
\alias{CV}
\alias{skew}
\alias{kurt}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Sample moments}
\description{
 \code{moments} provides the estimate of the first 4 moment-statistics of a sample.
}
\usage{
 moments (x)
 CV (x)
 skew (x)
 kurt (x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing a data-sample}
}
\details{
Skewness and kurtosis are defined as:
\deqn{skew = n^(-1) sd(x)^(-3) sum_i (x_i - mean(x))^3}
\deqn{kurt = n^(-1) sd(x)^(-4) sum_i (x_i - mean(x))^4 - 3}
where \eqn{n} is the size of \code{x}.

}
%\references{}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{mean}}, \code{\link{var}}, \code{\link{sd}}, \code{\link{Lmoments}}.}

\examples{
x <- rnorm(30,10,2)
moments(x)

data(hydroSIMN)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]
sapply(split(x,cod),moments)
}
\keyword{univar}
