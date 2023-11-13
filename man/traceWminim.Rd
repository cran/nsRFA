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
\deqn{d_{i j} = \sqrt{\frac{1}{p} \sum_{h=1}^{p} (x_{h i} - x_{h j})^2}}{d_ij = \sqrt{1/p sum[h from 1 to p](x_hi - x_hj)^2}}
where \eqn{x_{h i}}{x_hi} is the value of the \eqn{h}-th variable of the \eqn{i}-th element.

The function \code{traceWminim} is a composition of a jerarchical algorithm, the Ward (1963) one, and an optimisation procedure consisting in the minimisation of:
\deqn{W = \sum_{i=1}^k \left( \sum_{j=1}^{n_i} \delta_{i j}^2 \right)}{W = sum[i from 1 to k](sum[j from 1 to ni]delta_ij^2)}
where
\eqn{k} is the number of clusters (obtained initially with Ward's algorithm), \eqn{n_i}{ni} is the number of sites in the \eqn{i}-th cluster and \eqn{\delta_{i j}}{delta_ij} is the Euclidean distance between the \eqn{j}-th element of the \eqn{i}-th group and the center of mass of the \eqn{i}-th cluster.
\eqn{W} is calculated with \code{sumtraceW}.
The algorithm consist in moving a site from one cluster to another if this makes \eqn{W} decrease.
}
\value{
 \code{traceWminim} gives a vector defining the subdivision of elements characterized by \code{X} in n=\code{centers} clusters.

 \code{sumtraceW} gives \eqn{W} (it is used by \code{traceWminim}).

 \code{nearest} gives the nearest site to the centers of mass of clusters (it is used by \code{traceWminim}).
}

%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{roi}}, \code{\link{AD.dist}}.}

\examples{
data(hydroSIMN)
parameters
summary(parameters)

# traceWminim
param <- parameters[c("Hm","Ybar")]
n <- dim(param)[1]; k <- dim(param)[2]
param.norm <- (param - matrix(apply(param,2,mean),nrow=n,ncol=k,
               byrow=TRUE))/matrix(apply(param,2,sd),
               nrow=n,ncol=k,byrow=TRUE)
clusters <- traceWminim(param.norm,4); 
names(clusters) <- parameters["cod"][,]
clusters

annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]

fac <- factor(annualflows["cod"][,],
              levels=names(clusters[clusters==1]))
x1 <- annualflows[!is.na(fac),"dato"]
cod1 <- annualflows[!is.na(fac),"cod"]
#HW.tests(x1,cod1)          # it takes some time

fac <- factor(annualflows["cod"][,],
              levels=names(clusters[clusters==3]))
x3 <- annualflows[!is.na(fac),"dato"]
cod3 <- annualflows[!is.na(fac),"cod"]
#HW.tests(x3,cod3)          # it takes some time
}
\keyword{cluster}
