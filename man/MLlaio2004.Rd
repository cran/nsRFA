% 2005-09-23, Alberto Viglione
%
\name{MLlaio2004}
\alias{MLlaio2004}
\alias{ML_estimation}
\alias{moment_estimation}
%\alias{logLgumb}
%\alias{logLgev}
%\alias{logLgam}
%\alias{sample_generator}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Maximum likelihood parameters estimation}
\description{
 Maximum Likelihood estimation of parameters for extreme-value distributions, from Laio (2004).
}
\usage{
 ML_estimation (x, dist="NORM")
 moment_estimation (x, dist="NORM")
% logLgumb (T, x)
% logLgev (T, x)
% logLgam (T1, x)
% sample_generator (n, T, dist="NORM")
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{data sample}
 \item{dist}{distribution: normal \code{"NORM"}, Gumbel \code{"GUMBEL"}, Generalized Extreme Value \code{"GEV"}, Pearson type III \code{"P3"} and, only for \code{sample_generator}, Exponential \code{"EXP"}}
% \item{T}{parameters (position, scale, shape, ...)}
 %\item{w}{Transformed test statistic (Laio, 2004)}
% \item{T1}{position parameter of the Gamma (Pearson type III) distribution}
% \item{n}{sample length}
}
%\details{}
\value{
  \code{ML_estimation} estimate the parameters of the distribution \code{dist} from a sample \code{x} using the maximum likelihood approach. 

  \code{moment_estimation} estimate the parameters of the distribution \code{dist} from a sample \code{x} using the moment method.

%  \code{logLgumb} is the negative log-likelihood function for the Gumbel distribution.

%  \code{logLgev} is the negative log-likelihood function for the Generalized Extreme Value distribution.

%  \code{logLgam} is the negative log-likelihood function for the Pearson type III distribution.

%  \code{sample_generator} generates a sample of length \code{n} from the distribution \code{dist} with parameters \code{T}.

}
%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{GOFlaio2004}}.}

\examples{
# sample from an EV1 distribution
sm <- rand.gumb(100, 0, 1)
moment_estimation (sm, dist="GEV")
ML_estimation (sm, dist="GEV")

F.GEV(sm, -0.051, 0.97, -0.024)
rand.GEV (100, -0.051, 0.97, -0.024)
moment_estimation (sm, dist="P3")
ML_estimation (sm, dist="P3")
}
\keyword{htest}
