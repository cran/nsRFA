% 2006-06-16, Alberto Viglione
%
\name{CLUSTERS}
\alias{CLUSTERS}
\alias{traceWminim}
\alias{sumtraceW}
\alias{nearest}
\alias{AD.dist}
\alias{roi.hom}
\alias{roi.st.year}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Cluster analysis}
\description{
 Formation of clusters for Regional Frequency Analysis: disjoint regions (\code{traceWminim}) or region of influence (\code{roi.hom} and \code{roi.st.year}).
}
\usage{
 traceWminim (X, centers)
 sumtraceW (clusters, X)
 nearest (clusters, X)
 AD.dist (x, cod, index=2)
 roi.hom (p.ungauged, p.gauged, cod.p, x, cod,
   test="HW", limit=2, Nsim=500, index=2)
 roi.st.year (p.ungauged, p.gauged, cod.p, x, cod,
   test="HW", station.year=500, Nsim=500, index=2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{X}{a numeric matrix of characteristics, or an object that can be coerced to
          such a matrix (such as a numeric vector or a data frame with
          all numeric columns)}
 \item{centers}{the number of clusters}
 \item{clusters}{a numeric vector containing the subdivision of \code{X} in clusters}
 \item{x}{vector representing data from many samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{index}{if \code{index}=1 samples are divided by their average value;
              if \code{index}=2 (default) samples are divided by their median value}
 \item{p.ungauged}{parameters of the ungauged site (1 row)}
 \item{p.gauged}{parameters of gauged sites}
 \item{cod.p}{code of gauged sites}
 \item{test}{homogeneity test to apply: \code{"HW"} (default) or \code{"AD"}}
 \item{limit}{limit over which regions must be considered heterogeneous: for example 2 for \code{"HW"} or .95 for \code{"AD"}}
 \item{Nsim}{number of simulations in \code{"HW"} or \code{"AD"} tests}
 \item{station.year}{number of station years to form the region}
}
\details{In all cases, disjoint regions (Hosking and Wallis, 1997) or region of influence (Burn, 1990), the Euclidean distance is used.
Given \eqn{p} different classification variables, the distance between two elements \eqn{i} and \eqn{j} is:
\deqn{d_{i j} = \sqrt{\frac{1}{p} \sum_{h=1}^{p} (x_{h i} - x_{h j})^2}}
where \eqn{x_{h i}} is the value of the \eqn{h}-th variable of the \eqn{i}-th element.

The function \code{traceWminim} is a composition of a jerarchical algorithm, the Ward (1963) one, and an optimisation procedure consisting in the minimisation of:
\deqn{W = \sum_{i=1}^k \left( \sum_{j=1}^{n_i} \delta_{i j}^2 \right)}
where
\eqn{k} is the number of clusters (obtained initially with Ward's algorithm), \eqn{n_i} is the number of sites in the \eqn{i}-th cluster and \eqn{\delta_{i j}} is the Euclidean distance between the \eqn{j}-th element of the \eqn{i}-th group and the center of mass of the \eqn{i}-th cluster.
\eqn{W} is calculated with \code{sumtraceW}.
The algorithm consist in moving a site from one cluster to another if this makes \eqn{W} decrease.
}
\value{
 \code{traceWminim} gives a vector defining the subdivision of elements characterized by \code{X} in n=\code{centers} clusters.

 \code{sumtraceW} gives \eqn{W} (it is used by \code{traceWminim}).

 \code{nearest} gives the nearest site to the centers of mass of clusters (it is used by \code{traceWminim}).

 \code{AD.dist} returns the distance matrix between growth-curves built with the Anderson-Darling statistic.

 \code{roi.hom} returns the \sQuote{region of influence} for the site defined with \code{p.ungauged}.
 It returns codes of gauged sites that form an homogeneous region according to the Hosking and Wallis \code{"HW"} or Anderson-Darling \code{"AD"} tests.
 The region is formed using distances in the space defined by parameters \code{p.ungauged} and \code{p.gauged}.

 \code{roi.st.year} returns the \sQuote{region of influence} for the site defined with \code{p.ungauged}.
 It returns codes of gauged sites that form a region and the risult of homogeneity tests, according to the station-year criterion.
 The region is formed using distances in the space defined by parameters \code{p.ungauged} and \code{p.gauged}.

}
\references{
Burn D. (1990) Evaluation of regional flood frequency analysis with a region of influence approach. Water Resources Research, 26, pp. 2257-2265.

Everitt, B. (1974) Cluster Analysis. Social Science Research Council. Halsted Press, New York.

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Pre
ss, Cambridge, UK.

Viglione A., Laio F., Claps P. (2006) A comparison of homogeneity tests for regional frequency analysis, Water Resourches Research, In press.

Viglione A., Claps P., Laio F. (2006) Utilizzo di criteri di prossimit\`a nell'analisi regionale del deflusso annuo, XXX Conv
egno di Idraulica e Costruzioni Idrauliche - IDRA 2006, Roma, 10-15 Settembre 2006.

Viglione A. (2007) Metodi statistici non-supervised per la stima di grandezze idrologiche in siti non strumentati, PhD thesis, In press.

Ward J. (1963) Hierarchical grouping to optimize an objective function, Journal of the American Statistical Association, 58, pp. 236-244.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{HOMTESTS}} for the definition of the Hosking and Wallis \code{"HW"} or Anderson-Darling \code{"AD"} tests.}

\examples{
data(hydroSIMN)
parameters
summary(parameters)

# traceWminim
param <- parameters[c("Hm","Ybar")]
n <- dim(param)[1]; k <- dim(param)[2]
param.norm <- (param - matrix(mean(param),nrow=n,ncol=k,byrow=TRUE))/matrix(sd(param),
              nrow=n,ncol=k,byrow=TRUE)
clusters <- traceWminim(param.norm,4); 
names(clusters) <- parameters["cod"][,]
clusters

annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]

fac <- factor(annualflows["cod"][,],levels=names(clusters[clusters==1]))
x1 <- annualflows[!is.na(fac),"dato"]
cod1 <- annualflows[!is.na(fac),"cod"]
#HW.tests(x1,cod1)          # it takes some time

fac <- factor(annualflows["cod"][,],levels=names(clusters[clusters==3]))
x3 <- annualflows[!is.na(fac),"dato"]
cod3 <- annualflows[!is.na(fac),"cod"]
#HW.tests(x3,cod3)          # it takes some time

# Ad.dist
#AD.dist(x,cod)             # it takes some time

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
