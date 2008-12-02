% 2008-07-09, Alberto Viglione
%
\name{BayesianMCMC}
\alias{BayesianMCMC}
\alias{print.BayesianMCMC}
\alias{plot.BayesianMCMC}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Bayesian MCMC frequency analysis}
\description{
 Bayesian Markov Chain Monte Carlo algorithm for flood frequency analysis with historical and other information.
}
\usage{
 BayesianMCMC (xcont, xhist=NA, infhist=NA, suphist=NA, 
               nbans=NA, seuil=NA, nbpas=1000, nbchaines=3, 
               confint=c(0.05, 0.95), dist="GEV",
               apriori=function(...){1}, 
               parameters0=NA, varparameters0=NA)
 \method{plot}{BayesianMCMC} (x, which=1, ask=FALSE, ...)
 \method{print}{BayesianMCMC} (x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{object of class \code{BayesianMCMC}, output of function \code{BayesianMCMC}}
 \item{xcont}{vector of systematic data}
 \item{xhist}{vector of historical data}
 \item{infhist}{inferior limit for historical data}
 \item{suphist}{superior limit for historical data}
 \item{nbans}{period (in years) over which the threshold has been exceeded by historical data}
 \item{seuil}{threshold exceeded by historical data}
 \item{nbpas}{number of iterations for the MCMC algorithm}
 \item{nbchaines}{number of chains for the MCMC algorithm}
 \item{confint}{confidence limits for the flood quantiles}
 \item{dist}{distribution: normal \code{"NORM"}, log-normal with 2 parameters \code{"LN"}, Exponential \code{"EXP"}, Gumbel \code{"GUMBEL"}, Generalized Extreme Value \code{"GEV"}, Generalized Logistic \code{"GENLOGIS"}, Generalized Pareto \code{"GENPAR"}, log-normal with 3 parameters \code{"LN3"}, Pearson type III \code{"P3"}, log-Pearson type III \code{"LP3"}}
 \item{apriori}{function of the parameters of the model `proportional to' their a-priori guessed distribution.
                The default fuction returns always 1, i.e. there is no a-priori information}
 \item{parameters0}{initial values of the parameters for the MCMC algorithm}
 \item{varparameters0}{initial values of the parameter variances for the MCMC algorithm}
 \item{which}{a number of a vector of numbers that defines the graph to plot (see details)}
 \item{ask}{if TRUE, the interactive mode is run}
 \item{...}{other arguments}
}
\details{
\bold{Bayesian inference}

Bayesian inference uses a numerical estimate of the degree of belief in a hypothesis before evidence has been observed and calculates a numerical estimate of the degree of belief in the hypothesis after evidence has been observed.
The name `Bayesian' comes from the frequent use of Bayes' theorem in the inference process.
In our case the problem is: which is the probability that a frequency function \eqn{P} (of type defined in \code{dist}) has parameters \eqn{\theta}{T}, given that we have observed the realizations \eqn{D} (defined in \code{xcont}, \code{xhist}, \code{infhist}, \code{suphist}, \code{nbans}, \code{seuil}).
The Bayes' theorem writes
\deqn{P(\theta|D) = \frac{P(D|\theta) \cdot P(\theta)}{P(D)}}{P(T|D) = P(D|T)*P(T)/P(D)}
where 
\eqn{P(\theta|D)}{P(T|D)} is the conditional probability of \eqn{\theta}{T}, given \eqn{D} (it is also called the posterior probability because it is derived from or depends upon the specified value of \eqn{D}) and is the result we are interested in;
\eqn{P(\theta)}{P(T)} is the prior probability or marginal probability of \eqn{\theta}{T} (`prior' in the sense that it does not take into account any information about \eqn{D}), and can be given using the input \code{apriori} (it can be used to account for causal information); 
\eqn{P(D|\theta)}{P(D|T)} is the conditional probability of \eqn{D} given \eqn{\theta}{T} and it is defined choosing \code{dist} and depending on the availability of historical data;
\eqn{P(D)} is the prior or marginal probability of \eqn{D}, and acts as a normalizing constant.
Intuitively, Bayes' theorem in this form describes the way in which one's beliefs about observing \eqn{\theta}{T} are updated by having observed \eqn{D}.

Since complex models cannot be processed in closed form by a Bayesian analysis, namely because of the extreme difficulty in computing the normalization factor \eqn{P(D)}, simulation-based Monte Carlo techniques as the MCMC approaches are used.





\bold{MCMC Metropolis algorithm}

Markov chain Monte Carlo (MCMC) methods (which include random walk Monte Carlo methods), are a class of algorithms for sampling from probability distributions based on constructing a Markov chain that has the desired distribution as its equilibrium distribution. The state of the chain after a large number of steps is then used as a sample from the desired distribution.
The quality of the sample improves as a function of the number of steps.

The MCMC is performed here through a simple Metropolis algorithm, i.e. a Metropolis-Hastings algorithm with symmetric proposal density.
The Metropolis-Hastings algorithm can draw samples from any probability distribution \eqn{P(x)}, requiring only that a function proportional to the density can be calculated at \eqn{x}.
In Bayesian applications, the normalization factor is often extremely difficult to compute, so the ability to generate a sample without knowing this constant of proportionality is a major virtue of the algorithm.
The algorithm generates a Markov chain in which each state \eqn{x_t + 1}{xt+1} depends only on the previous state \eqn{x_t}{xt}.
The algorithm uses a Gaussian proposal density \eqn{N(x_t, \sigma_x)}{N(xt, Sx)}, which depends on the current state \eqn{x_t}{xt}, to generate a new proposed sample \eqn{x'}.
This proposal is accepted as the next value \eqn{x_t + 1 = x'}{xt + 1 = x'} if \eqn{\alpha}{a} drawn from \eqn{U(0,1)} satisfies
\deqn{\alpha < \frac{P(x')}{P(x_t)}}{a < P(x')/P(xt)}
If the proposal is not accepted, then the current value of \eqn{x} is retained (\eqn{x_t + 1 = x_t}{xt+1=xt}).

The Markov chain is started from a random initial value \eqn{x_0}{x0} and the algorithm is run for many iterations until this initial state is forgotten. 
These samples, which are discarded, are known as burn-in. 
The remaining set of accepted values of \eqn{x} represent a sample from the distribution \eqn{P(x)}.
As a Gaussian proposal density (or a lognormal one for definite-positive parameters) is used, the variance parameter \eqn{\sigma_x^2}{Sx^2} has to be tuned during the burn-in period.
This is done by calculating the acceptance rate, which is the fraction of proposed samples that is accepted in a window of the last \eqn{N} samples. 
The desired acceptance rate depends on the target distribution, however it has been shown theoretically that the ideal acceptance rate for a one dimensional Gaussian distribution is approx 50\%, decreasing to approx 23\% for an N-dimensional Gaussian target distribution.
If \eqn{\sigma_x^2}{Sx^2} is too small the chain will mix slowly (i.e., the acceptance rate will be too high, so the sampling will move around the space slowly and converge slowly to \eqn{P(x)}).
If \eqn{\sigma_x^2}{Sx^2} is too large the acceptance rate will be very low because the proposals are likely to land in regions of much lower probability density.
The desired acceptance rate is fixed here to 34\%.


The MCMC algorithm is based on a code developed by Eric Gaume on Scilab.
It is still unstable and not all the distributions have been tested.
}
\value{
  \code{BayesianMCMC} returns the following values:

   \code{parameters} matrix (nbpas)x(nbchaines) with the simulated sets of parameters with the MCMC algorithm;

   \code{parametersML} set of parameters correspondent to the maximum likelihood;

   \code{returnperiods} return periods for which \code{quantilesML} and \code{intervals} are calculated;

   \code{quantilesML} quantiles correspondent to \code{returnperiods} for the distribution whose parameters are \code{parametersML};

   \code{intervals} confidence intervals for the quantiles \code{quantilesML} for limits \code{confint};

   \code{varparameters} matrix (nbpas)x(nbchaines)x(number of parameters) with the simulated variances for the MCMC algorithm;

   \code{vraisdist} likelihoods for the sets \code{parameters};




  \code{plot.BayesianMCMC} plots the following figures:

   \code{1} data as plotting position, fitted distribution (maximum likelihood) and confidence intervals;

   \code{2} diagnostic plot of the MCMC simulation (parameters);

   \code{3} diagnostic plot of the MCMC simulation (likelyhood and MCMC acceptance rate);

   \code{4} posterior distribution of parameters obtained with the MCMC simulation (cloud plots);

   \code{5} a-priori distribution of parameters (contour plots);
}


%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}


\seealso{.}

\examples{
set.seed(2988)
serie <- rand.GEV(120, xi=40, alfa=20, k=-0.4)
serie100 <- serie[1:100]
serie100[serie100 < 250] <- NA
serie20 <- serie[101:120]
serie <- c(serie100, serie20)


plot(serie, type="h", ylim=c(0, 600), xlab="", 
     ylab="Annual flood peaks [m3/s]", lwd=3)
abline(h=0)
points(serie100, col=2)

# Using only sistematic data
only_sist <- BayesianMCMC (xcont=serie20, xhist=NA, infhist=NA, suphist=NA, 
                           nbans=NA, seuil=NA,
                           nbpas=5000, nbchaines=3, 
                           confint=c(0.05, 0.95), dist="GEV")
plot(only_sist, which=c(1:3), ask=TRUE, ylim=c(1,600))



# Adding the information that the threshold 250 m3/s was exceeded 
#   3 times in the past 100 years
with_hist_thresh <- BayesianMCMC (xcont=serie20, xhist=NA, infhist=rep(250,3), 
                                  suphist=NA, nbans=100, seuil=250,
                                  nbpas=5000, nbchaines=3, 
                                  confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_thresh, which=c(1:3), ask=TRUE, ylim=c(1,600))



# Assuming that the 3 historical events are known with high uncertainty
with_hist_limits <- BayesianMCMC (xcont=serie20, xhist=NA, 
                                  infhist=c(320,320,250), 
                                  suphist=c(360,400,270), 
                                  nbans=100, seuil=250,
                                  nbpas=5000, nbchaines=3, 
                                  confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_limits, which=c(1:3), ask=TRUE, ylim=c(1,600))



# Assuming that the 3 historical events are perfectly known
with_hist_known <- BayesianMCMC (xcont=serie20, xhist=serie100[!is.na(serie100)], 
                                 infhist=NA, suphist=NA, 
                                 nbans=100, seuil=250,
                                 nbpas=5000, nbchaines=3, 
                                 confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_known, which=c(1:3), ask=TRUE, ylim=c(1,600))




# Using one reasonable a-priori distribution
fNORM3 <- function (x) {
 # x = vector of values
 # mu = vector of means
 mu = c(44, 26, -0.40)
 # CM = covariance matrix
 CM = matrix(c(13, 7.8, -0.055,
               7.8, 15, -0.42,
               -0.055, -0.42, 0.056), nrow=3, ncol=3)
 CMm1 <- solve(CM)
 term2 <- exp(-((x - mu) \%*\% CMm1 \%*\% (x - mu))/2)
 term1 <- 1/(2*pi)^(3/2)/sqrt(det(CM))
 term1*term2
}

with_hist_known2 <- BayesianMCMC (xcont=serie20, xhist=serie100[!is.na(serie100)], 
                                  infhist=NA, suphist=NA,
                                  nbans=100, seuil=250,
                                  nbpas=5000, nbchaines=3, apriori=fNORM3,
                                  confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_known2, 5)
plot(with_hist_known2, 4)
plot(with_hist_known, 4)
plot(with_hist_known)
plot(with_hist_known2)


\dontrun{
# Using one non-reasonable a-priori distribution
fNORM3 <- function (x) {
 # x = vector of values
 # mu = vector of means
 mu = c(30, 50, -0.10)
 # CM = covariance matrix
 CM = matrix(c(13, 7.8, -0.055,
               7.8, 15, -0.42,
               -0.055, -0.42, 0.056), nrow=3, ncol=3)
 CMm1 <- solve(CM)
 term2 <- exp(-((x - mu) \%*\% CMm1 \%*\% (x - mu))/2)
 term2
}

with_hist_known3 <- BayesianMCMC (xcont=serie20, xhist=serie100[!is.na(serie100)], 
                                  infhist=NA, suphist=NA,
                                  nbans=100, seuil=250,
                                  nbpas=5000, nbchaines=3, apriori=fNORM3,
                                  confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_known3, 5)
plot(with_hist_known3, 4)
plot(with_hist_known, 4)
plot(with_hist_known)
plot(with_hist_known3)
}
}
\keyword{htest}
