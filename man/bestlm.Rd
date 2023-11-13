% 2008-08-28, Alberto Viglione
%
\name{bestlm}
\alias{bestlm}
\alias{print.bestlm}
\alias{summary.bestlm}
\alias{.leaps.setup}
\alias{.leaps.exhaustive}
\alias{.leaps.backward}
\alias{.leaps.forward}
\alias{.leaps.seqrep}
\alias{.summary.regsubsets}
\alias{.leaps}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Subsets regression}
\description{
 \code{bestlm} performs an exhaustive search for the best subsets of the variables in \code{ind} 
 for predicting \code{dip} in linear regression.
}
\usage{
 bestlm (dip, ind, kmax=4, nbest=3)
 \method{print}{bestlm} (x, ...)
 \method{summary}{bestlm} (object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{object of class \code{bestlm}, output of function \code{bestlm}}
 \item{object}{object of class \code{bestlm}, output of function \code{bestlm}}
 \item{dip}{vector n x 1 of dependent variable to be predicted}
 \item{ind}{matrix n x K of the K independent variables (candidate predictors)}
 \item{kmax}{maximum size (number of regressors) to report}
 \item{nbest}{number of subsets of each size to report}
 \item{...}{other arguments}
}
\details{
This function has been obtained using the function \code{leaps} of the R package \code{leaps}.
It is based on the Alan Miller's FORTRAN routines.

Warning:
the function will stop with an error if \code{ind} is not of full rank or if it has more than 31 columns.
}
\value{
  \code{bestlm} returns the following values:

   \code{subselect} matrix (kmax*nbest)x(ncol(ind)) with the sets of chosen linear models, ordered
                    in function of the adjusted coefficient of determination (R2adj);

   \code{R2adj} the ordered adjusted coefficient of determination;
}


%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}


\seealso{\code{\link{lm}}, \code{\link{REGRDIAGNOSTICS}}.}

\examples{
data(hydroSIMN)

bestlm(parameters[,"Dm"], parameters[,-c(1:2)])

regr <- lm(Dm ~ Am + S2000 + NORD, parameters); regr
summary(regr)
}
\keyword{regression}
