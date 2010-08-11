% 2005-09-23, Alberto Viglione
%
\name{STATICPLOTS}
\alias{STATICPLOTS}
\alias{Lmoment.ratio.diagram}
\alias{Lspace.HWvsAD}
\alias{Lspace.limits}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Static plots}
\description{
 Plots from books and articles.
}
\usage{
 Lmoment.ratio.diagram (grid=TRUE, ...)
 Lspace.HWvsAD (grid=TRUE, ...)
 Lspace.limits (grid=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{grid}{should a grid be plotted?}
 \item{...}{other arguments}
}
%\details{}
\value{
 \code{Lmoment.ratio.diagram} plots points corresponding to two parameters distributions and lines corresponding to three parameters distributions on the 'L-CA - L-kur' plane.
 The distributions are:
 E = exponential,
 G = gumbel,
 L = logistic,
 N = normal,
 U = uniform,
 GLO = generalized logistic,
 GEV = generalized extreme-value,
 GPA = generalized Pareto,
 LN3 = lognormal,
 PE3 = Pearson type III.

 \code{Lspace.HWvsAD} separate regions, in 'L-CA - L-CV' space, where the homogeneity tests of Hosking and Wallis (HW) and Anderson-Darling (AD) are preferable.

 \code{Lspace.limits} displays limits for regional L-moments in the 'L-CA - L-CV'.
}

%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{EXP}}, \code{\link{GENLOGIS}}, \code{\link{GENPAR}}, \code{\link{LOGNORM}}, \code{\link{GUMBEL}}, \code{\link{GEV}}, \code{\link{P3}}}

\examples{
Lmoment.ratio.diagram()
Lspace.HWvsAD()
Lspace.limits()

data(hydroSIMN)
annualflows[c(1:10),]
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]
rlm <- regionalLmoments(x,cod)
Lmoment.ratio.diagram()
points(rlm["lcaR"],rlm["lkurR"],col="red",pch=19)

Lspace.HWvsAD()
points(rlm["lcaR"],rlm["lcvR"],col="red",pch=19)
}
\keyword{hplot}
