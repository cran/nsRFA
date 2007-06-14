% 2005-09-23, Alberto Viglione
%
\name{HW.original}
\alias{HW.original}
\alias{print.HWorig}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Original Hosking and Wallis Fortran routine}
\description{
 The original Fortran routine by Hosking is here used to analyse a region.
}
\usage{
 HW.original (data, cod, Nsim=500)
 print.HWorig (x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{object of class \code{HWorig}}
 \item{data}{vector representing data from many samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{Nsim}{number of regions simulated with the bootstrap of the original region}
 \item{...}{additional parameter for \code{print}}
}
\details{
}
\value{
  \code{HW.original} returns an object of class \code{HWorig} (what the Fortran subroutine REGTST return).
}
\references{
Hosking J., Wallis J. (1993) Some statistics useful in regional frequency analysis. Water Resources Research, 29 (2), pp. 271-281.

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{HW.tests}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]
split(x,cod)

HW.original(x,cod)

fac <- factor(annualflows["cod"][,],levels=c(34:38))
x2 <- annualflows[!is.na(fac),"dato"]
cod2 <- annualflows[!is.na(fac),"cod"]

HW.original(x2,cod2)
}
\keyword{htest}
