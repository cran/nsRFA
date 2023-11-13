% 2007-03-01, Alberto Viglione
%
\name{AD.dist}
\alias{AD.dist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Anderson-Darling distance matrix for growth curves}
\description{
 Distance matrix for growth curves. Every element of the matrix is the Anderson-Darling statistic calculated between two series.
}
\usage{
 AD.dist (x, cod, index=2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing data from many samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{index}{if \code{index}=1 samples are divided by their average value;
              if \code{index}=2 (default) samples are divided by their median value}
}
\details{The Anderson-Darling statistic used here is the one defined in \url{https://en.wikipedia.org/wiki/Anderson-Darling_test} as \eqn{A^2}{A2}.}
\value{
 \code{AD.dist} returns the distance matrix between growth-curves built with the Anderson-Darling statistic.
}
%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{traceWminim}}, \code{\link{roi}}.}

\examples{
data(hydroSIMN)

annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]

# Ad.dist
AD.dist(x,cod)             # it takes some time
}
\keyword{cluster}
