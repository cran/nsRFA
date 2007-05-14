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
 Lmoment.ratio.diagram ()
 Lspace.HWvsAD ()
 Lspace.limits ()
}
%- maybe also `usage' for other objects documented here.
%\arguments{\item{}{}}
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
\references{
Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.


Viglione A., Laio F., Claps P. (2006) A comparison of homogeneity tests in Regional Frequency Analysis. Water Resourches Rese
arch, In press.
}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

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
