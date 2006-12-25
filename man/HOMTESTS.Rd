% 2005-09-23, Alberto Viglione
%
\name{HOMTESTS}
\alias{HOMTESTS}
\alias{ADbootstrap.test}
\alias{HW.tests}
\alias{DK.test}
\alias{discordancy}
\alias{criticalD}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Homogeneity tests}
\description{
 Homogeneity tests for Regional Frequency Analysis.
}
\usage{
 ADbootstrap.test (x, cod, Nsim=500, index=2)
 HW.tests (x, cod, Nsim=500)
 DK.test (x, cod)
 discordancy (x, cod)
 criticalD ()
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing data from many samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{Nsim}{number of regions simulated with the bootstrap of the original region}
 \item{index}{if \code{index}=1 samples are divided by their average value;
              if \code{index}=2 (default) samples are divided by their median value}
}
\details{
\bold{The Hosking and Wallis heterogeneity measures}

The idea underlying Hosking and Wallis (1993) heterogeneity
statistics is to measure the sample variability of the L-moment
ratios and compare it to the variation that would be expected in a
homogeneous region. The latter is estimated through repeated
simulations of homogeneous regions with samples drawn from a four
parameter kappa distribution (see e.g., Hosking and Wallis,
1997, pp. 202-204). 
More in detail, the steps are the following:
with regards to the \eqn{k} samples belonging to the region under analysis, find
the sample L-moment ratios (see, Hosking and Wallis, 1997)
pertaining to the \eqn{i}-th site: these are the
L-coefficient of variation (L-CV),
\deqn{t^{(i)}=\frac{\frac{1}{n_i}\sum_{j=1}^{n_i}\left(
\frac{2(j-1)}{(n_i-1)}-1\right)Y_{i,j}}{\frac{1}{n_i}\sum_{j=1}^{n_i}Y_{i,j}}}
the coefficient of L-skewness,
\deqn{t_3^{(i)}=\frac{\frac{1}{n_i}\sum_{j=1}^{n_i}\left(
\frac{6(j-1)(j-2)}{(n_i-1)(n_i-2)}-\frac{6(j-1)}{(n_i-1)}+1\right)Y_{i,j}}{\frac{1}{n_i}\sum_{j=1}^{n_i}\left(
\frac{2(j-1)}{(n_i-1)}-1\right)Y_{i,j}}}
and the coefficient of L-kurtosis
\deqn{t_4^{(i)}=\frac{\frac{1}{n_i}\sum_{j=1}^{n_i}\left(\frac{20(j-1)(j-2)(j-3)}{(n_i-1)(n_i-2)(n_i-3)}-
\frac{30(j-1)(j-2)}{(n_i-1)(n_i-2)}+\frac{12(j-1)}{(n_i-1)}-1\right)Y_{i,j}}{\frac{1}{n_i}\sum_{j=1}^{n_i}\left(
\frac{2(j-1)}{(n_i-1)}-1\right)Y_{i,j}}}
Note that the L-moment ratios are not affected by the
normalization by the index value, i.e. it is the same to use
\eqn{X_{i,j}} or \eqn{Y_{i,j}} in Equations.

Define the regional averaged L-CV, L-skewness and
L-kurtosis coefficients,
\deqn{t^R = \frac{\sum_{i=1}^k n_i t^{(i)}}{ \sum_{i=1}^k n_i}}
\deqn{t_3^R =\frac{ \sum_{i=1}^k n_i t_3^{(i)}}{ \sum_{i=1}^k n_i}}
\deqn{t_4^R =\frac{ \sum_{i=1}^k n_i t_4^{(i)}}{\sum_{i=1}^k n_i}}
and compute the statistic
\deqn{V = \left\{ \sum_{i=1}^k n_i (t^{(i)} - t^R )^2 / \sum_{i=1}^k n_i\right\} ^{1/2}}

Fit the parameters of a
four-parameters kappa distribution to the regional averaged L-moment ratios
\eqn{t^R}, \eqn{t_3^R} and \eqn{t_4^R}, and then generate a large number
\eqn{N_{sim}} of realizations of sets of \eqn{k} samples. The \eqn{i}-th site sample in each set 
has a kappa distribution as its parent and
record length equal to \eqn{n_i}. For each simulated
homogeneous set, calculate the statistic \eqn{V}, obtaining \eqn{N_{sim}} values. 
On this vector of \eqn{V} values determine the mean \eqn{\mu_V} and standard
deviation \eqn{\sigma_V} that relate to the hypothesis of homogeneity
(actually, under the composite hypothesis of homogeneity and kappa
parent distribution).

An heterogeneity measure, which is called here
\eqn{HW_1}, is finally found as
\deqn{\theta_{HW_1} = \frac{V - \mu_V}{\sigma_V}}
\eqn{\theta_{HW_1}} can be approximated by a normal distributed with zero
mean and unit variance: following Hosking and Wallis (1997),
the region under analysis can therefore be regarded as
\sQuote{acceptably homogeneous} if \eqn{\theta_{HW_1}<1}, \sQuote{possibly
heterogeneous} if \eqn{1 \leq \theta_{HW_1} < 2}, and \sQuote{definitely
heterogeneous} if \eqn{\theta_{HW_1}\geq2}. Hosking and Wallis
(1997) suggest that these limits should be treated as useful
guidelines. Even if the \eqn{\theta_{HW_1}} statistic is constructed
like a significance test, significance levels obtained from such a
test would in fact be accurate only under special assumptions: to have
independent data both serially and between sites, and the true
regional distribution being kappa.

Hosking and Wallis (1993) also give an alternative heterogeneity measure
(that we call \eqn{HW_2}), in which \eqn{V} is
replaced by:
\deqn{V_2 = \sum_{i=1}^k n_i \left\{ (t^{(i)} - t^R)^2 + (t_3^{(i)} -
t_3^R)^2\right\}^{1/2} / \sum_{i=1}^k n_i}
The test statistic in this case becomes
\deqn{\theta_{HW_2} = \frac{V_2 - \mu_{V_2}}{\sigma_{V_2}}}
with similar acceptability limits as the \eqn{HW_1} statistic. 
Hosking and Wallis (1997) judge \eqn{\theta_{HW_2}} to be inferior to
\eqn{\theta_{HW_1}} and say that it rarely yields values larger than 2 even for grossly heterogeneous regions.




\bold{The bootstrap Anderson-Darling test}

A test that does not make any assumption on the parent distribution is the
Anderson-Darling (\eqn{AD}) rank test (Scholz and Stephens, 1987).
The \eqn{AD} test is the generalization of the classical
Anderson-Darling goodness of fit test (e.g., D'Agostino and
Stephens, 1986), and it is used to test the hypothesis that \eqn{k}
independent samples belong to the same population without
specifying their common distribution function.

The test is based on the comparison between local and regional
empirical distribution functions. The empirical distribution
function, or sample distribution function, is defined by
\eqn{F(x)=\frac{j}{\eta}, x_{(j)}\leq x < x_{(j+1)}}, where \eqn{\eta} is
the size of the sample and \eqn{x_{(j)}} are the order statistics,
i.e. the observations arranged in ascending order. Denote the
empirical distribution function of the $i$-th sample (local) by \eqn{\hat{F}_i
(x)}, and that of the pooled sample of all \eqn{N = n_1 + ... + n_k}
observations (regional) by \eqn{H_N (x)}. The \eqn{k}-sample Anderson-Darling test
statistic is then defined as
\deqn{\theta_{AD} = \sum_{i=1}^k n_i \int _{{\rm all}\ x} \frac{[
\hat{F}_i (x) - H_N (x) ]^2}{H_N (x) [ 1 - H_N (x) ] } dH_N (x)}

If the pooled ordered sample is \eqn{Z_1 < ... < Z_N}, the
computational formula to evaluate \eqn{\theta_{AD}} is:
\deqn{\theta_{AD} = \frac{1}{N} \sum_{i=1}^k \frac{1}{n_i}
\sum_{j=1}^{N-1} \frac{(N M_{ij} - j n_i)^2 }{j (N-j)}}
where \eqn{M_{ij}} is the number of observations in the \eqn{i}-th sample
that are not greater than \eqn{Z_j}. The homogeneity test can be
carried out by comparing the obtained \eqn{\theta_{AD}} value to the
tabulated percentage points reported by Scholz and Stephens
(1987) for different significance levels.

The statistic \eqn{\theta_{AD}} depends on the sample values only
through their ranks. This guarantees that the test statistic
remains unchanged when the samples undergo monotonic
transformations, an important stability property not possessed by
\eqn{HW} heterogeneity measures. However, problems arise in applying this test in a
common index value procedure. In fact, the index
value procedure corresponds to dividing each site sample by a different
value, thus modifying the ranks in the pooled sample. In
particular, this has the effect of making the
local empirical distribution functions much more similar to the
other, providing an impression of homogeneity even when the
samples are highly heterogeneous. The effect is analogous to that
encountered when applying goodness-of-fit tests to distributions
whose parameters are estimated from the same sample used for the
test (e.g., D'Agostino and Stephens, 1986; Laio,
2004). In both cases, the percentage points for the test should be
opportunely redetermined. This can be done with a nonparametric bootstrap approach
presenting the following steps:
build up the pooled sample \eqn{\cal S} of the observed
non-dimensional data.
Sample with replacement from \eqn{\cal S} and generate \eqn{k}
artificial local samples, of size \eqn{n_1, \dots ,n_k}.
Divide each sample for its index value, and calculate
\eqn{\theta^{(1)}_{AD}}.
Repeat the procedure for \eqn{N_{sim}} times and obtain a sample
of \eqn{\theta^{(j)}_{AD}}, \eqn{j=1,\dots , N_{sim}} values, whose
empirical distribution function can be used as an approximation of
\eqn{G_{H_0}(\theta_{AD})}, the distribution of \eqn{\theta_{AD}} under
the null hypothesis of homogeneity.
The acceptance limits for the test, corresponding to any
significance level \eqn{\alpha}, are then easily determined as the
quantiles of \eqn{G_{H_0}(\theta_{AD})} corresponding to a probability
\eqn{(1-\alpha)}.

We will call the test obtained with the above procedure the bootstrap Anderson-Darling test, hereafter referred to as \eqn{AD}.



\bold{Durbin and Knott test}

The last considered homogeneity test derives from a
goodness-of-fit statistic originally proposed by Durbin and
Knott (1971). The test is formulated to measure discrepancies in
the dispersion of the samples, without accounting for the possible
presence of discrepancies in the mean or skewness of the data.
Under this aspect, the test is similar to the \eqn{HW_1} test, while it
is analogous to the \eqn{AD} test for the fact that it is a rank test.
The original goodness-of-fit test is very simple: suppose to have
a sample \eqn{X_i}, \eqn{i = 1, ..., n}, with hypothetical
distribution \eqn{F(x)}; under the null hypothesis the random variable
\eqn{F(X_i)} has a uniform distribution in the \eqn{(0,1)} interval, and
the statistic \eqn{D = \sum_{i=1}^n \cos[2 \pi F(X_i)]} is
approximately normally distributed with mean 0 and variance 1
(Durbin and Knott, 1971). \eqn{D} serves the purpose of
detecting discrepancy in data dispersion: if the variance of \eqn{X_i}
is greater than that of the hypothetical distribution \eqn{F(x)}, \eqn{D} is significantly greater than
0, while \eqn{D} is significantly below 0 in the reverse case.
Differences between the mean (or the median) of \eqn{X_i} and \eqn{F(x)}
are instead not detected by \eqn{D}, which guarantees that the
normalization by the index value does not affect the test.

The extension to homogeneity testing of the  Durbin and
Knott (\eqn{DK}) statistic is straightforward: we substitute the empirical
distribution function obtained with the pooled observed data,
\eqn{H_N(x)}, for \eqn{F(x)} in \eqn{D}, obtaining at each site a statistic
\deqn{D_i = \sum_{j=1}^{n_i} \cos[2 \pi H_N(X_j)]}
which is normal under the hypothesis of homogeneity. The statistic
\eqn{\theta_{DK} = \sum_{i=1}^k D_i^2} has then a chi-squared
distribution with \eqn{k-1} degrees of freedom, which allows one to
determine the acceptability limits for the test, corresponding to
any significance level \eqn{\alpha}. 




\bold{Comparison among tests}

The comparison (Viglione et al, 2006) shows that the Hosking and Wallis heterogeneity measure \eqn{HW_1} (only based on L-CV) is preferable when skewness is low, while the bootstrap Anderson-Darling test should be used for more skewed regions.
As for \eqn{HW_2}, the Hosking and Wallis heterogeneity measure based on L-CV and L-CA, it is shown once more how much it lacks power.

Our suggestion is to guide the choice of the test according to a compromise between power and Type I error of the \eqn{HW_1} and \eqn{AD} tests.
The L-moment space is divided into two regions: 
if the \eqn{t_3^R} coefficient for the region under analysis is lower than 0.23, we propose to use the Hosking and Wallis heterogeneity measure \eqn{HW_1};
if \eqn{t_3^R > 0.23}, the bootstrap Anderson-Darling test is preferable.




}
\value{
  \code{ADbootstrap.test} and \code{DK.test} test gives its test statistic and its distribution value \eqn{P}.
  If \eqn{P} is, for example, 0.92, samples shouldn't be considered heterogeneous with significance level minor of 8%.

  \code{HW.tests} gives the two Hosking and Wallis heterogeneity measures \eqn{H_1} and \eqn{H_2}; following Hosking and Wallis (1997), the region under analysis can therefore be regarded as \sQuote{acceptably homogeneous} if \eqn{H_1<1}, \sQuote{possibly heterogeneous} if \eqn{1 \leq H_1 < 2}, and \sQuote{definitely heterogeneous} if \eqn{H \geq 2}.

  \code{discordancy} returns the discordancy measure \eqn{D} of Hosking and Wallis for all sites. 
  Hosking and Wallis suggest to consider a site discordant if \eqn{D \geq 3} if \eqn{N \geq 15} (where \eqn{N} is the number of sites considered in the region). For \eqn{N<15} the critical values of \eqn{D} can be listed with \code{criticalD}.
}
\references{
D'Agostino R., Stephens M. (1986) Goodness-of-Fit Techniques, chapter Tests based on EDF statistics. Marcel Dekker, New York.

Durbin J., Knott M. (1971) Components of Cramer-von Mises statistics. London School of Economics and Political Science, pp. 290-307.

Hosking J., Wallis J. (1993) Some statistics useful in regional frequency analysis. Water Resources Research, 29 (2), pp. 271-281.

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.

Laio, F., Cramer-von Mises and Anderson-Darling goodness of fit tests for extreme value distributions with unknown parameters
, Water Resour. Res., 40, W09308, doi:10.1029/2004WR003204.

Scholz F., Stephens M. (1987) K-sample Anderson-Darling tests. Journal of American Statistical Association, 82 (399), pp. 918
-924.

Viglione A., Laio F., Claps P. (2006) A comparison of homogeneity tests in Regional Frequency Analysis. Water Resourches Rese
arch, In press.

Viglione A. (2007) Metodi statistici non-supervised per la stima di grandezze idrologiche in siti non strumentati, PhD thesis
, In press.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{CLUSTERS}}, \code{\link{KAPPA}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]
split(x,cod)

#ADbootstrap.test(x,cod,Nsim=100)   # it takes some time
#HW.tests(x,cod)                    # it takes some time
DK.test(x,cod)

fac <- factor(annualflows["cod"][,],levels=c(34:38))
x2 <- annualflows[!is.na(fac),"dato"]
cod2 <- annualflows[!is.na(fac),"cod"]

ADbootstrap.test(x2,cod2,Nsim=100)
ADbootstrap.test(x2,cod2,index=1,Nsim=200)
HW.tests(x2,cod2,Nsim=100)
DK.test(x2,cod2)

discordancy(x,cod)

criticalD()
}
\keyword{htest}
