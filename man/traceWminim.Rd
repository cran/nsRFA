% 2007-03-01, Alberto Viglione
%
\name{traceWminim}
\alias{traceWminim}
\alias{sumtraceW}
\alias{nearest}
\title{Cluster analysis: disjoint regions}
\description{
 Formation of disjoint regions for Regional Frequency Analysis.
}
\usage{
 traceWminim (X, centers)
 sumtraceW (clusters, X)
 nearest (clusters, X)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{X}{a numeric matrix of characteristics, or an object that can be coerced to
          such a matrix (such as a numeric vector or a data frame with
          all numeric columns)}
 \item{centers}{the number of clusters}
 \item{clusters}{a numeric vector containing the subdivision of \code{X} in clusters}
}
\details{The Euclidean distance is used.
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
}
\references{
Everitt, B. (1974) Cluster Analysis. Social Science Research Council. Halsted Press, New York.

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Pre
ss, Cambridge, UK.

Viglione A., Claps P., Laio F. (2006) Utilizzo di criteri di prossimit\`a nell'analisi regionale del deflusso annuo, XXX Convegno di Idraulica e Costruzioni Idrauliche - IDRA 2006, Roma, 10-15 Settembre 2006.

Viglione A. (2007) Metodi statistici non-supervised per la stima di grandezze idrologiche in siti non strumentati, PhD thesis, In press.

Ward J. (1963) Hierarchical grouping to optimize an objective function, Journal of the American Statistical Association, 58, pp. 236-244.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{roi}}, \code{\link{AD.dist}}.}

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
}
\keyword{cluster}
