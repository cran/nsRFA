% 2005-09-23, Alberto Viglione
%
\name{SERIESPLOTS}
\alias{SERIESPLOTS}
\alias{consistencyplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Series plots}
\description{
 Plots for time-series investigation.
}
\usage{
 consistencyplot (x, cod, cex.axis=.8, mark.code=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing time (e.g. years) of data-samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{cex.axis}{dimensions of points and labels}
 \item{mark.code}{if TRUE (default) codes of samples are plotted on y axis}
 \item{...}{graphical parameters as \code{xlab}, \code{ylab}, \code{main}, ...}
}
%\details{}
\value{
 \code{consistencyplot} displays time-series consistency.
}
%\references{}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{plot}}, \code{\link{DISTPLOTS}}}

\examples{

data(hydroSIMN)
annualflows[c(1:10),]
y <- annualflows["anno"][,]
cod <- annualflows["cod"][,]
consistencyplot(y,cod)
}
\keyword{hplot}
