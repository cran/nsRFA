% 2005-09-23, Alberto Viglione
%
\name{MLlaio2004}
\alias{MLlaio2004}
\alias{ML_estimation}
\alias{moment_estimation}
\alias{logLgumb}
\alias{logLgev}
\alias{logLgam}
\alias{sample_generator}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Maximum likelihood parameters estimation}
\description{
 Maximum Likelihood estimation of parameters for extreme-value distributions, from Laio (2004).
}
\usage{
 ML_estimation (x, dist="NORM")
 moment_estimation (x, dist="NORM")
 logLgumb (T, x)
 logLgev (T, x)
 logLgam (T1, x)
 sample_generator (n, T, dist="NORM")
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{data sample}
 \item{dist}{distribution: normal \code{"NORM"}, Gumbel \code{"GUMBEL"}, Generalized Extreme Value \code{"GEV"}, Pearson type III \code{"P3"} and, only for \code{sample_generator}, Exponential \code{"EXP"}}
 \item{T}{parameters (position, scale, shape, ...)}
 %\item{w}{Transformed test statistic (Laio, 2004)}
 \item{T1}{position parameter of the Gamma (Pearson type III) distribution}
 \item{n}{sample length}
}
%\details{}
\value{
  \code{ML_estimation} estimate the parameters of the distribution \code{dist} from a sample \code{x} using the maximum likelihood approach. 

  \code{moment_estimation} estimate the parameters of the distribution \code{dist} from a sample \code{x} using the moment method.

  \code{logLgumb} is the negative log-likelihood function for the Gumbel distribution.

  \code{logLgev} is the negative log-likelihood function for the Generalized Extreme Value distribution.

  \code{logLgam} is the negative log-likelihood function for the Pearson type III distribution.

  \code{sample_generator} generates a sample of length \code{n} from the distribution \code{dist} with parameters \code{T}.

}
\references{
Laio, F., Cramer-von Mises and Anderson-Darling goodness of fit tests for extreme value distributions with unknown parameters, Water Resour. Res., 40, W09308, doi:10.1029/2004WR003204.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{GOFlaio2004}}.}

\examples{
sm <- sample_generator(100, c(0,1), dist="EV1")
moment_estimation (sm, dist="GEV")
ML_estimation (sm, dist="GEV")
Fx (sm, c(-0.051,0.97,-0.024), dist="GEV")
sample_generator (100, c(-0.051,0.97,-0.024), dist="GEV")

moment_estimation (sm, dist="P3")
ML_estimation (sm, dist="P3")
Fx (sm, c(-1.6,0.78,2.8), dist="P3")

}
\keyword{htest}
