% 2005-09-23, Alberto Viglione
%
\name{DIAGNOSTICS}
\alias{DIAGNOSTICS}
\alias{R2}
\alias{RMSE}
\alias{MAE}
\alias{RMSEP}
\alias{MAEP}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Diagnostics of models}
\description{
 Diagnostics of model results, it compares estimated values \code{y} with observed values \code{x}. 
}
\usage{
 R2 (x, y, na.rm=FALSE)
 RMSE (x, y, na.rm=FALSE) 
 MAE (x, y, na.rm=FALSE)
 RMSEP (x, y, na.rm=FALSE)
 MAEP (x, y, na.rm=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{observed values}
 \item{y}{estimated values}
 \item{na.rm}{logical. Should missing values be removed?}
}
\details{
 If \eqn{x_i}{xi} are the observed values, \eqn{y_i}{yi} the estimated values, with \eqn{i=1,...,n}, and \eqn{\bar{x}}{<x>} the sample mean of \eqn{x_i}{xi}, then:
  \deqn{R^2 = 1 - \frac{\sum_1^n (x_i-y_i)^2}{\sum_1^n x_i^2 - n \bar{x}^2}}{R2 = 1 - (sum(xi-yi)^2)/(sum(xi^2-n<x>^2)}
  \deqn{RMSE = \sqrt{\frac{1}{n} \sum_1^n (x_i-y_i)^2}}{RMSE = sqrt(1/n sum(xi-yi)^2)}
  \deqn{MAE = \frac{1}{n} \sum_1^n |x_i-y_i|}{MAE = 1/n sum(|xi-yi|)}
  \deqn{RMSEP = \sqrt{\frac{1}{n} \sum_1^n ((x_i-y_i)/x_i)^2}}{RMSEP = sqrt(1/n sum((xi-yi)/xi)^2)}
  \deqn{MAEP = \frac{1}{n} \sum_1^n |(x_i-y_i)/x_i|}{MAEP = 1/n sum(|(xi-yi)/xi|}
 See \url{http://en.wikipedia.org/wiki/Coefficient_of_determination}, \url{http://en.wikipedia.org/wiki/Mean_squared_error} and \url{http://en.wikipedia.org/wiki/Mean_absolute_error} for other details.
}
\value{
 \code{R2} returns the coefficient of determination \eqn{R^2}{R2} of a model.

 \code{RMSE} returns the root mean squared error of a model.

 \code{MAE} returns the mean absolute error of a model.

 \code{RMSE} returns the percentual root mean squared error of a model.

 \code{MAE} returns the percentual mean absolute error of a model.
}
%\references{}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{lm}}, \code{\link{summary.lm}}, \code{\link{predict.lm}}, \code{\link{REGRDIAGNOSTICS}}}

\examples{
data(hydroSIMN)

datregr <- parameters
regr0 <- lm(Dm ~ .,datregr); summary(regr0)
regr1 <- lm(Dm ~ Am + Hm + Ybar,datregr); summary(regr1)

obs <- parameters[,"Dm"]
est0 <- regr0$fitted.values
est1 <- regr1$fitted.values

R2(obs, est0)
R2(obs, est1)

RMSE(obs, est0)
RMSE(obs, est1)

MAE(obs, est0)
MAE(obs, est1)

RMSEP(obs, est0)
RMSEP(obs, est1)

MAEP(obs, est0)
MAEP(obs, est1)
}
\keyword{models}
