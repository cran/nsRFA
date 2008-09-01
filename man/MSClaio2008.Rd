% 2008-08-19, Alberto Viglione
%
\name{MSClaio2008}
\alias{MSClaio2008}
\alias{print.MSClaio2008}
\alias{plot.MSClaio2008}
\alias{summary.MSClaio2008}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Model Selection Criteria}
\description{
 Model selection criteria for the frequency analysis of hydrological extremes, from Laio et al (2008).
}
\usage{
 MSClaio2008 (sample, dist=c("NORM","LN","GUMBEL","EV2","GEV","P3","LP3"), 
              crit=c("AIC", "AICc", "BIC", "ADC"))
 \method{print}{MSClaio2008} (x, digits=max(3, getOption("digits") - 3), ...)
 \method{summary}{MSClaio2008} (object, ...)
 \method{plot}{MSClaio2008} (x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{sample}{data sample}
 \item{dist}{distributions: normal \code{"NORM"}, 2 parameter log-normal \code{"LN"}, Gumbel \code{"GUMBEL"}, Frechet \code{"EV2"}, Generalized Extreme Value \code{"GEV"}, Pearson type III \code{"P3"}, log-Pearson type III \code{"LP3"}}
 \item{crit}{Model-selection criteria: Akaike Information Criterion \code{"AIC"}, Akaike Information Criterion corrected \code{"AICc"}, Bayesian Information Criterion \code{"BIC"}, Anderson-Darling Criterion \code{"ADC"}}
 \item{x}{object of class \code{MSClaio2008}, output of \code{MSClaio2008()}}
 \item{object}{object of class \code{MSClaio2008}, output of \code{MSClaio2008()}}
 \item{digits}{minimal number of "significant" digits, see 'print.default'}
 \item{...}{other arguments}
}
\details{
The following lines are extracted from Laio et al. (2008). See the paper for more details and references.

\bold{Model selection criteria}

The problem of model selection can be formalized as follows: a sample of \eqn{n} data, \eqn{D=(x_1, \dots, x_n)}{D=(x1, ..., xn)}, arranged in ascending order is available, sampled from an unknown parent distribution \eqn{f(x)}; 
\eqn{N_m}{Nm} operating models, \eqn{M_j}{Mj}, \eqn{j=1,\dots, N_m}{j=1, ..., Nm}, are used to represent the data. 
The operating models are in the form of probability distributions, \eqn{M_j = g_j(x,\hat{\theta})}{Mj = gj(x, O*)}, with parameters \eqn{\hat{\theta}}{O*} estimated from the available data sample \eqn{D}. 
The scope of model selection is to identify the model \eqn{M_{opt}}{Mopt} which is better suited to represent the data, i.e. the model which is closer in some sense to the parent distribution \eqn{f(x)}. 

Three different model selection criteria are considered here, namely, the Akaike Information Criterion (AIC), the Bayesian Information Criterion (BIC), and the Anderson-Darling Criterion (ADC).
Of the three methods, the first two belong to the category of classical literature approaches, while the third derives from a heuristic interpretation of the results of a standard goodness-of-fit test (see Laio, 2004).


\bold{Akalike Information Criterion}

The Akaike information Criterion (AIC) for the j-th operational model can be computed as
\deqn{AIC_j = -2 ln (L_j(\hat{\theta})) + 2 p_j}{%
      AICj = -2 ln (Lj(O*)) + 2 pj}
where
\deqn{L_j(\hat{\theta}) = \prod_{i=1}^n g_j(x_i, \hat{\theta})}{%
      Lj(O*) = prod (gj(xi, O*))}
is the likelihood function, evaluated at the point \eqn{\theta=\hat{\theta}}{O = O*} corresponding to the maximum likelihood estimator of the parameter vector \eqn{\theta}{O} and \eqn{p_j}{pj} is the number of estimated parameter of the j-th operational model.
In practice, after the computation of the \eqn{AIC_j}{AICj}, for all of the operating models, one selects the model with the minimum AIC value, \eqn{AIC_{min}}{AICmin}.

When the sample size, \eqn{n}, is small, with respect to the number of estimated parameters, \eqn{p}, the AIC may perform inadequately. In those cases a second-order variant of AIC, called AICc, should be used:
\deqn{AICc_j = -2 ln (L_j(\hat{\theta})) + 2 p_j (n/(n - p_j - 1))}{%
      AICcj = -2 ln (Lj(O*)) + 2 pj (n/(n - pj - 1))}
Indicatively, AICc should be used when \eqn{n/p < 40}.


\bold{Bayesian Information Criterion}

The Bayesian Information Criterion (BIC) for the j-th operational model reads
\deqn{BIC_j = -2 ln (L_j(\hat{\theta})) + ln(n) p_j}{%
      BICj = -2 ln (Lj(O*)) + ln(n) pj}
In practical application, after the computation of the \eqn{BIC_j}{BICj}, for all of the operating models, one selects the model with the minimum BIC value, \eqn{BIC_{min}}{BICmin}.


\bold{Anderson-Darling Criterion}

The Anderson-Darling criterion has the form:
\deqn{ADC_j = 0.0403 + 0.116 ((\Delta_{AD,j} - \epsilon_j)/\beta_j)^{(\eta_j/0.851)}}{%
      ADCj = 0.0403 + 0.116 ((ADj - epsj)/betaj)^(etaj/0.851)}
if \eqn{1.2 \epsilon_j < \Delta_{AD,j}}{1.2 epsj < ADj},
\deqn{ADC_j = [0.0403 + 0.116 ((0.2 \epsilon_j)/\beta_j)^{(\eta_j/0.851)}] (\Delta_{AD,j} - 0.2 \epsilon_j / \epsilon_j)}{%
      ADCj = [0.0403 + 0.116 ((0.2 epsj)/betaj)^(etaj/0.851) (ADj - 0.2 epsj / epsj}
if \eqn{1.2 \epsilon_j \ge \Delta_{AD,j}}{1.2 epsj >= ADj},
where \eqn{\Delta_{AD,j}}{ADj} is the discrepancy measure characterizing the criterion, the Anderson-Darling statistic \code{A2} in \code{\link{GOFlaio2004}}, and \eqn{\epsilon_j}{epsj}, \eqn{\beta_j}{betaj} and \eqn{\eta_j}{etaj} are distribution-dependent coefficients that are tabled by Laio [2004, Tables 3 and 5] for a set of seven distributions commonly employed for the frequency analysis of extreme events. 
In practice, after the computation of the \eqn{ADC_j}{ADCj}, for all of the operating models, one selects the model with the minimum ADC value, \eqn{ADC_{min}}{ADCmin}.
}
\value{
 \code{MSClaio2008} returns the value of the criteria \code{crit} (see Details) chosen applied to the \code{sample}, for every distribution \code{dist}.

 \code{plot.MSClaio2008} plots the empirical distribution function of \code{sample} (Weibull plotting position) on a log-normal probability plot, plots the candidate distributions \code{dist} (whose parameters are evaluated with the maximum likelihood technique, see \code{\link{MLlaio2004}}, and highlights the ones chosen by the criteria \code{crit}.)
}

%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{\code{\link{GOFlaio2004}}, \code{\link{MLlaio2004}}.}

\examples{
data(FEH1000)

sitedata <- am[am[,1]==53004, ] # data of site 53004
serieplot(sitedata[,4], sitedata[,3])
MSC <- MSClaio2008(sitedata[,4])
MSC
summary(MSC)
plot(MSC)

sitedata <- am[am[,1]==69023, ]	# data of site 69023
serieplot(sitedata[,4], sitedata[,3])
MSC <- MSClaio2008(sitedata[,4], crit=c("AIC", "ADC"))
MSC
summary(MSC)
plot(MSC)

sitedata <- am[am[,1]==83802, ] # data of site 83802
serieplot(sitedata[,4], sitedata[,3])
MSC <- MSClaio2008(sitedata[,4], dist=c("GEV", "P3", "LP3"))
MSC
summary(MSC)
plot(MSC)

# short sample, high positive L-CA
sitedata <- am[am[,1]==40012, ] # data of site 40012
serieplot(sitedata[,4], sitedata[,3])
MSC <- MSClaio2008(sitedata[,4])
MSC
summary(MSC)
plot(MSC)

# negative L-CA
sitedata <- am[am[,1]==68002, ] # data of site 68002
serieplot(sitedata[,4], sitedata[,3])
MSC <- MSClaio2008(sitedata[,4])
MSC
summary(MSC)
plot(MSC)

}
\keyword{htest}
