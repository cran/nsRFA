% 2007-03-12, Alberto Viglione
%
\name{varLmoments}
\alias{varLmoments}
\alias{varLCV}
\alias{varLCA}
\alias{varLkur}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Exact variance structure of sample L-moments}
\description{
 \code{varLmoments} provides distribution-free unbiased estimators of the variances and covariances of sample L-moments.
}
\usage{
 varLmoments (x, matrix=TRUE)
 varLCV (x)
 varLCA (x)
 varLkur (x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing a data-sample}
 \item{matrix}{if \code{TRUE} (default), the matrix of estimates of the variance structure (variance and covariance) i
               of sample L-moments is returned; if \code{FALSE}, a vector containing \eqn{var(l_1)}{var(l1)}, 
               \eqn{var(l_2)}{var(l2)}, \eqn{var(l_3)}{var(l3)}, \eqn{var(l_4)}{var(l4)}, \eqn{var(t)}, 
               \eqn{var(t_3)}{var(t3)} and \eqn{var(t_4)}{var(t4)} is returned.}
}
\details{
The estimation of the exact variance structure of sample L-moments is based on Elamir et Seheult (2004).
}
\value{
 \code{varLmoments} gives the matrix of unbiased estimates of the variance structure of sample L-moments: 
  this is a 4x4 matrix containg \eqn{var(l_1)}{var(l1)}, \eqn{var(l_2)}{var(l2)}, \eqn{var(l_3)}{var(l3)}, 
  \eqn{var(l_4)}{var(l4)} on the main diagonal, 
  and the correspondant covariances elsewhere (\eqn{cov(l_1,l_2)}{cov(l1,l2)}, \eqn{cov(l_1,l_3)}{cov(l1,l3)}, etc.);

 \code{varLCV} gives the unbiased estimate of the variance of sample coefficient of L-variation of \code{x};

 \code{varLCA} gives the unbiased estimate of the variance of sample L-skewness of \code{x};

 \code{varLkur} gives the unbiased estimate of the variance of sample L-kurtosis of \code{x}.
}
\references{Elamir, E.A.H. and Seheult, A.H. (2004) Exact variance structure of sample L-moments, Journal of Statistical Planning and Inference 124, 337-359.}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{var}}, \code{\link{Lmoments}}.}

\examples{
x <- rnorm(30,10,2)
varLmoments(x)
varLmoments(x, FALSE)

varLCV(x)
varLCA(x)
varLkur(x)

data(hydroSIMN)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]
dvarLmom <- function(x) {diag(varLmoments(x))}
sapply(split(x,cod),dvarLmom)

}
\keyword{univar}
