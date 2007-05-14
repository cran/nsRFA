% 2007-03-01, Alberto Viglione
%
\name{AD.dist}
\alias{AD.dist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Distance matrix for growth curves}
\description{
 Distance matrix for growth curves.
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
%\details{}
\value{
 \code{AD.dist} returns the distance matrix between growth-curves built with the Anderson-Darling statistic.
}
\references{
Viglione A., Laio F., Claps P. (2006) A comparison of homogeneity tests for regional frequency analysis, Water Resourches Research, In press.

Viglione A., Claps P., Laio F. (2006) Utilizzo di criteri di prossimit\`a nell'analisi regionale del deflusso annuo, XXX Convegno di Idraulica e Costruzioni Idrauliche - IDRA 2006, Roma, 10-15 Settembre 2006.

Viglione A. (2007) Metodi statistici non-supervised per la stima di grandezze idrologiche in siti non strumentati, PhD thesis, In press.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

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
