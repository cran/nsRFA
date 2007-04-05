% 2005-09-23, Alberto Viglione
%
\name{Lmoments}
\alias{Lmoments}
\alias{regionalLmoments}
\alias{LCV}
\alias{LCA}
\alias{Lkur}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Hosking and Wallis sample L-moments}
\description{
 \code{Lmoments} provides the estimate of L-moments of a sample or regional L-moments of a region.
}
\usage{
 Lmoments (x)
 regionalLmoments (x,cod)
 LCV (x)
 LCA (x)
 Lkur (x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing a data-sample (or data from many samples defined with \code{cod} in the case of \code{regionalLmoments})}
 \item{cod}{array that defines the data subdivision among sites}
}
\details{
The estimation of L-moments is based on a sample of size \eqn{n}, arranged in ascending order.
Let \eqn{x_{1:n} \le x_{2:n} \le \dots \le x_{n:n}}{x(1:n) <= x(2:n) <= ... <= x(n:n)} be the ordered sample.
An unbiased estimator of the probability weighted moments \eqn{\beta_r}{\betar} is:
\deqn{b_r = n^{-1} \sum_{j=r+1}^n \frac{(j-1)(j-2)\dots(j-r)}{(n-1)(n-2)\dots(n-r)} x_{j:n}}{%
      br = 1/n sum[j from r+1 to n](x(j:n) (j-1)(j-2)...(j-r)/(n-1)/(n-2)/.../(n-r))}

The sample L-moments are defined by:
\deqn{l_1 = b_0}{l1 = b0}
\deqn{l_2 = 2b_1 - b_0}{l2 = 2b1 - b0}
\deqn{l_3 = 6b_2 - 6b_1 + b_0}{l3 = 6b2 - 6b1 + b0}
\deqn{l_4 = 20b_3-30b_2+12b_1-b_0}{l4 = 20b3 - 30b2 + 12b1 - b0}
and in general
\deqn{l_{r+1} = \sum_{k=0}^r \frac{(-1)^{r-k}(r+k)!}{(k!)^2(r-k)!} b_k}{%
      l(r+1) = sum[k from 0 to r](b_k (-1)^(r-k) (r+k)! / (k!)^2 / (r-k)!)}
where \eqn{r=0,1,\dots,n-1}{r=0, 1, ..., n-1}.

The sample L-moment ratios are defined by
\deqn{t_r=l_r/l_2}{tr = lr / l2}
and the sample L-CV by
\deqn{t=l_2/l_1}{t = l2 / l1}

Sample regional L-CV, L-skewness and L-kurtosis coefficients are defined as
\deqn{t^R = \frac{\sum_{i=1}^k n_i t^{(i)}}{ \sum_{i=1}^k n_i}}{%
      t^R = sum[i from 1 to k](ni t^(i)) / sum[i from 1 to k](ni)}
\deqn{t_3^R =\frac{ \sum_{i=1}^k n_i t_3^{(i)}}{ \sum_{i=1}^k n_i}}{%
      t3^R = sum[i from 1 to k](ni t3^(i)) / sum[i from 1 to k](ni)}
\deqn{t_4^R =\frac{ \sum_{i=1}^k n_i t_4^{(i)}}{\sum_{i=1}^k n_i}}{%
      t4^R = sum[i from 1 to k](ni t4^(i)) / sum[i from 1 to k](ni)}
}
\value{
 \code{Lmoments} gives the L-moments (\eqn{l_1}, \eqn{l_2}, \eqn{t}, \eqn{t_3}, \eqn{t_4}), \code{regionalLmoments} gives the regional weighted L-moments (\eqn{l_1^R}, \eqn{l_2^R}, \eqn{t^R}, \eqn{t_3^R}, \eqn{t_4^R}), \code{LCV} gives the coefficient of L-variation, \code{LCA} gives the L-skewness and \code{Lkur} gives the L-kurtosis of \code{x}.
}
\references{Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{mean}}, \code{\link{var}}, \code{\link{sd}}, \code{\link{HOMTESTS}}.}

\examples{
x <- rnorm(30,10,2)
Lmoments(x)

data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]
split(x,cod)
camp <- split(x,cod)$"45"
Lmoments(camp)
sapply(split(x,cod),Lmoments)

regionalLmoments(x,cod)
}
\keyword{univar}
