% 2005-09-23, Alberto Viglione
%
\name{SERIESPLOTS}
\alias{SERIESPLOTS}
\alias{serieplot}
\alias{consistencyplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Series plots}
\description{
 Plots for time-series investigation.
}
\usage{
 serieplot (x, t, lim.x=c(min(x),max(x)), lim.t=c(min(t),max(t)), 
            ...)
 consistencyplot (t, cod, cex.axis=.8, mark.code=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{data sample}
 \item{t}{vector representing time (e.g. years) of data-samples defined with \code{cod}}
 \item{lim.x, lim.t}{plot limits}
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
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{plot}}, \code{\link{DISTPLOTS}}}

\examples{

data(hydroSIMN)
annualflows[c(1:10),]
x <- annualflows["dato"][,]
y <- annualflows["anno"][,]
cod <- annualflows["cod"][,]
consistencyplot(y,cod)


for (i in unique(cod)) {
 serieplot(x[cod==i], y[cod==i], c(0,max(x)), c(min(y),max(y)),
           xlab="", ylab="D [mm]", main=i)
 readline()
}
}
\keyword{hplot}
