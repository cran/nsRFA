% 2008-07-09, Alberto Viglione
%
\name{BayesianMCMC}
\alias{BayesianMCMC}
\alias{print.BayesianMCMC}
\alias{plot.BayesianMCMC}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Bayesian MCMC frequency analysis}
\description{
 Bayesian Markov Chain Monte Carlo algorithm for flood frequency analysis with historical information.
}
\usage{
 BayesianMCMC (xcont, xhist=NA, infhist=NA, suphist=NA, nbans=NA, seuil=NA,
               nbpas=1000, nbchaines=3, confint=c(0.05, 0.95), dist="GEV",
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
 \item{dist}{distribution: normal \code{"NORM"}, log-normal with 2 parameters \code{"LN"}, Exponential \code{"EXP"}, Gumbel \code{"GUMBEL"}, Generalized Extreme Value \code{"GEV"}, Generalized Logistic \code{"GENLOGIS"}, Generalized Pareto \code{"GENPAR"}, log-normal with 3 parameters \code{"LOGNORM"}, Pearson type III \code{"P3"}, log-Pearson type III \code{"LP3"}}
 \item{parameters0}{initial values of the parameters for the MCMC algorithm}
 \item{varparameters0}{initial values of the parameter variances for the MCMC algorithm}
 \item{which}{a number of a vector of numbers that defines the graph to plot (see details)}
 \item{ask}{if TRUE, the interactive mode is run}
 \item{...}{other arguments}
}
\details{The algorithm is based on a code developed by Eric Gaume on Scilab.
         It is still unstable and not all the distributions have been tested.}
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

   \code{3} posterior distribution of parameters obtained with the MCMC simulation (cloud plots);
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


plot(serie, type="h", ylim=c(0, 600), xlab="", ylab="Annual flood peaks [m3/s]", lwd=3)
abline(h=0)
points(serie100, col=2)

# Using only sistematic data
only_sist <- BayesianMCMC (xcont=serie20, xhist=NA, infhist=NA, suphist=NA, 
                           nbans=NA, seuil=NA,
                           nbpas=5000, nbchaines=3, 
                           confint=c(0.05, 0.95), dist="GEV")
plot(only_sist, which=c(1:3), ask=TRUE, ylim=c(1,600))



# Adding the information that the threshold 250 m3/s was exceeded 3 times in the past 100 years
with_hist_thresh <- BayesianMCMC (xcont=serie20, xhist=NA, infhist=rep(250,3), suphist=NA, 
                                  nbans=100, seuil=250,
                                  nbpas=5000, nbchaines=3, 
                                  confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_thresh, which=c(1:3), ask=TRUE, ylim=c(1,600))



# Assuming that the 3 historical events are known with high uncertainty
with_hist_limits <- BayesianMCMC (xcont=serie20, xhist=NA, infhist=c(320,320,250), suphist=c(360,400,270), 
                                  nbans=100, seuil=250,
                                  nbpas=5000, nbchaines=3, 
                                  confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_limits, which=c(1:3), ask=TRUE, ylim=c(1,600))



# Assuming that the 3 historical events are perfectly known
with_hist_known <- BayesianMCMC (xcont=serie20, xhist=serie100[!is.na(serie100)], infhist=NA, suphist=NA, 
                                 nbans=100, seuil=250,
                                 nbpas=5000, nbchaines=3, 
                                 confint=c(0.05, 0.95), dist="GEV")
plot(with_hist_known, which=c(1:3), ask=TRUE, ylim=c(1,600))


}
\keyword{htest}
