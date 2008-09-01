% 2005-09-23, Alberto Viglione
%
\name{nsRFA-internal}
\alias{nsRFA-internal}
\alias{ksampleA2}
\alias{nonparboot}
\alias{test.GOFmontecarlo}
\alias{typeIerrorA2_GOFlaio}
\alias{fw2}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Internal functions}
\description{
 User-level objects which are for \sQuote{internal} use only

 \code{ksampleA2} provides the Anderson-Darling test statistic;

 \code{nonparboot} provides a non-parametric bootstrap;

 \code{test.GOFmontecarlo} is used to test the type I error of goodness of fit tests \code{GOFmontecarlo}.

 \code{typeIerrorA2_GOFlaio} is used to test the type I error of goodness of fit tests \code{A2_GOFlaio}.

 \code{fw2} is the approximation of the probability distribution of \code{w} (first 2 terms) when \eqn{H_0} is true (Anderson-Darling, 1952); it is used by \code{A2_GOFlaio}.
}
\usage{
 ksampleA2 (x,cod)
 nonparboot (z,n=length(z))
 test.GOFmontecarlo (parameters, type="NORM", alfa=.05, n=30, N=100)
 typeIerrorA2_GOFlaio (n, T, alfa=0.05, dist="NORM", Nsim=1000)
 fw2 (w)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing data from many samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{z}{data sample, used for bootstrap}
 \item{n}{length of sample (extracted in \code{nonparboot})}
 \item{parameters}{vector of distribution parameters}
 \item{type}{distribution type: \code{"NORM"}, \code{"GENLOGIS"}, \code{"GENPAR"}, \code{"GEV"}, \code{"LOGNORM"}, \code{"P3"}}
 \item{alfa}{significance level of the test}
 \item{N}{number of iterations}
 \item{T}{parameters (position, scale, shape, ...)}
 \item{dist}{distribution (\code{"NORM"}, \code{"LN"}, \code{"EV1"}, \code{"EV2"}, \code{"GEV"}, \code{"GAM"}, \code{"LP3"})}
 \item{Nsim}{number of iterations, as \code{N}}
 \item{w}{transformed test statistic (Laio, 2004)}
}
%\details{}

%\note{.}
%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

%\seealso{\code{\link{factor}}, \code{\link{HOMTESTS}}, \code{\link{traceWminim}}, \code{\link{roi}}.}

%\examples{}
\keyword{internal}
