% 2007-03-01, Alberto Viglione
%
\name{roi}
\alias{roi}
\alias{roi.hom}
\alias{roi.st.year}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Region of influence}
\description{
 Formation of clusters for Regional Frequency Analysis: region of influence (Burn, 1990).
}
\usage{
 roi (p.ungauged, p.gauged, cod.p, x=NULL, cod=NULL)
 roi.hom (p.ungauged, p.gauged, cod.p, x, cod,
   test="HW", limit=2, Nsim=500, index=2)
 roi.st.year (p.ungauged, p.gauged, cod.p, x, cod,
   test="HW", station.year=500, Nsim=500, index=2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing data from many samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{index}{if \code{index}=1 samples are divided by their average value;
              if \code{index}=2 (default) samples are divided by their median value}
 \item{p.ungauged}{parameters of the ungauged site (1 row)}
 \item{p.gauged}{parameters of gauged sites}
 \item{cod.p}{code of gauged sites}
 \item{test}{homogeneity test to apply: \code{"HW"} (default) or \code{"AD"} (in \code{roi.st.year} you can choose \code{"HW and AD"} too}
 \item{limit}{limit over which regions must be considered heterogeneous: for example 2 for \code{"HW"} or .95 for \code{"AD"}}
 \item{Nsim}{number of simulations in \code{"HW"} or \code{"AD"} tests}
 \item{station.year}{number of station years to form the region}
}
\details{The Euclidean distance is used.
Given \eqn{p} different classification variables, the distance between two elements \eqn{i} and \eqn{j} is:
\deqn{d_{i j} = \sqrt{\frac{1}{p} \sum_{h=1}^{p} (x_{h i} - x_{h j})^2}}
where \eqn{x_{h i}} is the value of the \eqn{h}-th variable of the \eqn{i}-th element.
}
\value{
 \code{roi} returns the \sQuote{region of influence} for the site defined with \code{p.ungauged}.
 It the gauged sites ordered according to the euclidean distance against the site of interest (the distance is evaluated in the space defined by parameters \code{p.ungauged} and \code{p.gauged}).
 If \code{x=NULL} and \code{cod=NULL} (default), a data.frame with the ordered sites and the distances against the site of interest is returned.
 If \code{x} and \code{cod} are provided, the data.frame will contain also statistics of samples (number of data \code{n} and L-moments).

 \code{roi.hom} returns the \sQuote{region of influence} for the site defined with \code{p.ungauged}.
 It returns codes of gauged sites that form an homogeneous region according to the Hosking and Wallis \code{"HW"} or Anderson-Darling \code{"AD"} tests.
 The region is formed using distances in the space defined by parameters \code{p.ungauged} and \code{p.gauged}.

 \code{roi.st.year} returns the \sQuote{region of influence} for the site defined with \code{p.ungauged}.
 It returns codes of gauged sites that form a region and the risult of homogeneity tests, according to the station-year criterion.
 It also return the similarity ranking factor \eqn{S_i}, the weights \eqn{w_i} and the regional L-moments as evaluated in the Flood Estimation Handbook (Robson and Reed, 1999).
 The region is formed using distances in the space defined by parameters \code{p.ungauged} and \code{p.gauged}.
}
\references{
Burn D. (1990) Evaluation of regional flood frequency analysis with a region of influence approach. Water Resources Research, 26, pp. 2257-2265.

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Pre
ss, Cambridge, UK.

Viglione A. (2007) Metodi statistici non-supervised per la stima di grandezze idrologiche in siti non strumentati, PhD thesis, In press.

Robson, A. and Reed, D. (1999) Statistical procedures for flood frequency estimation. In Flood Estimation Handbook, volume 3. Institute of Hydrology Crowmarsh Gifford, Wallingford, Oxfordshire. 
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{traceWminim}}, \code{\link{AD.dist}}, \code{\link{HOMTESTS}} for the definition of the Hosking and Wallis \code{"HW"} or Anderson-Darling \code{"AD"} tests.}

\examples{
data(hydroSIMN)
parameters
summary(parameters)

annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]

roi(parameters[5,3:5],parameters[-5,3:5],parameters[-5,1])
roi(parameters[5,3:5],parameters[-5,3:5],parameters[-5,1],x,cod)

# roi.hom
#roi.hom(parameters[5,3:5],parameters[-5,3:5],parameters[-5,1],x,cod)
                            # it takes some time
#roi.hom(parameters[5,3:5],parameters[-5,3:5],parameters[-5,1],x,cod,
#        test="AD",limit=.95)      # it takes some time

#roi.hom(parameters[8,3:5],parameters[-8,3:5],
#         parameters[-8,1],x,cod)    # it takes some time


# roi.st.year
roi.st.year(parameters[5,3:5],parameters[-5,3:5],
            parameters[-5,1],x,cod)
roi.st.year(parameters[5,3:5],parameters[-5,3:5],parameters[-5,1],
            x,cod,test="AD",station.year=100)

}
\keyword{cluster}
