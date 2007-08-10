% 2005-09-23, Alberto Viglione
%
\name{HW.original}
\alias{HW.original}
\alias{print.HWorig}
\alias{plot.HWorig}
\alias{LMR}
\alias{PEL}
\alias{SAMLMR}
\alias{SAMLMU}
\alias{SAMPWM}
\alias{REGLMR}
\alias{REGTST}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Original Hosking and Wallis Fortran routine}
\description{
 The original Fortran routine by Hosking is here used to analyse a region.
}
\usage{
 HW.original (data, cod, Nsim=500)
 print.HWorig (x, ...)
 plot.HWorig (x, interactive=TRUE, ...)
 LMR (PARA, distr="EXP")
 PEL (XMOM, distr="EXP")
 SAMLMR (X, A=0, B=0)
 SAMLMU (X)
 SAMPWM (X, A=0, B=0)
 REGLMR (data, cod)
 REGTST (data, cod, A=0, B=0, Nsim=500)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{object of class \code{HWorig}}
 \item{data}{vector representing data from many samples defined with \code{cod}}
 \item{cod}{array that defines the data subdivision among sites}
 \item{Nsim}{number of regions simulated with the bootstrap of the original region}
 \item{interactive}{logical: if TRUE the graphic showing is interactive}
 \item{...}{additional parameter for \code{print}}
 \item{PARA}{parameters of the distribution (vector)}
 \item{distr}{distribution:
              \code{EXP} = Exponential (2 parameters: xi, alfa);
              \code{GAM} = Gamma (2 parameters: alfa, beta);
              \code{GEV} = Generalized extreme value (3 parameters: xi, alfa, k);
              \code{GLO} = Generalized logistic (3 parameters: xi, alfa, k);
              \code{GNO} = Generalized Normal (3 parameters: xi, alfa, k);
              \code{GPA} = Generalized Pareto (3 parameters: xi, alfa, k);
              \code{GUM} = Gumbel (2 parameters: xi, alfa);
              \code{KAP} = Kappa (4 parameters: xi, alfa, k, h);
              \code{NOR} = Normal (2 parameters: mu, sigma);
              \code{PE3} = Pearson type III (3 parameters: mu, sigma, gamm);
              \code{WAK} = Wakeby (5 parameters: xi, alfa, beta, gamm, delta).}
 \item{XMOM}{the L-moment ratios of the distribution, in order \eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}, \eqn{\tau_5}...}
 \item{X}{a data vector}
 \item{A, B}{Parameters of plotting position: for unbiased estimates (of the \eqn{\lambda}'s) set A=B=zero. 
             Otherwise, plotting-position estimators are used, based on the plotting position \eqn{(j+a)/(n+b)} 
             for the \eqn{j}'th smallest of \eqn{n} observations. 
             For example,  A=-0.35 and B=0.0 yields the estimators recommended by 
             Hosking et al. (1985, technometrics) for the GEV distribution.}
}
\details{
Documentation of the original Fortran routines by Hosking available at
\url{http://www.research.ibm.com/people/h/hosking/lmoments.html}.

Differences among \code{HW.original} and \code{HW.tests} should depend on differences among \code{PEL} and \code{par.kappa} for the kappa distribution.
A numerical algorithm is used to resolve the implicit Equations (A.99) and (A.100) in Hosking and Wallis (1997, pag. 203-204).
The algorithms in \code{PEL} and \code{par.kappa} are different.
Anyway the risults of the tests should converge asymptotically.
}
\value{
 \code{HW.original} returns an object of class \code{HWorig} (what the Fortran subroutine REGTST return).

 \code{LMR} calculates the L-moment ratios of a distribution given its parameters.

 \code{PEL} calculates the parameters of a distribution given its L-moments.

 \code{SAMLMR} calculates the sample L-moments ratios of a data-set.

 \code{SAMLMU} calculates the `unbiased' sample L-moments ratios of a data-set.

 \code{SAMPWM} calculates the sample probability weighted moments of a data-set.

 \code{REGLMR} calculates regional weighted averages of the sample L-moments ratios.

 \code{REGTST} calculates statistics useful in regional frequency analysis.
 1) Discordancy measure, d(i), for individual sites in a region.
    Large values might be used as a flag to indicate potential errors
    in the data at the site.  "large" might be 3 for regions with 15
    or more sites, but less (exact values in array dc1) for smaller
    regions.
 2) Heterogeneity measures, H(j), for a region based upon either:-
    j=1: the weighted s.d. of the l-cvs or
    j=2: the average distance from the site to the regional average
         on a graph of l-cv vs. l-skewness
    j=3: the average distance from the site to the regional average
         on a graph of l-skewness vs. l-kurtosis.
    In practice H(1) is probably sufficient.  a value greater than
    (say) 1.0 suggests that further subdivision of the region should
    be considered as it might improve quantile estimates.
 3) Goodness-of-fit measures, Z(k), for 5 candidate distributions:
    k=1: generalized logistic
    k=2: generalized extreme value
    k=3: generalized normal (lognormal)
    k=4: pearson type iii (3-parameter gamma)
    k=5: generalized pareto.
    Provided that the region is acceptably close to homogeneous,
    the fit may be judged acceptable at 10% significance level
    if Z(k) is less than 1.645 in absolute value.

 For further details see Hosking and Wallis (1997),
 "Regional frequency analysis: an approach based on L-moments",
 cambridge university press, chapters 3-5.  
}
\references{
Hosking J., Wallis J. (1993) Some statistics useful in regional frequency analysis. Water Resources Research, 29 (2), pp. 271-281.

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.
}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

\note{
\bold{IBM software disclaimer}

LMOMENTS: Fortran routines for use with the method of L-moments

Permission to use, copy, modify and distribute this software for any purpose and without fee is hereby granted, 
provided that this copyright and permission notice appear on all copies of the software. 
The name of the IBM Corporation may not be used in any advertising or publicity pertaining to the use 
of the software. IBM makes no warranty or representations about the suitability of the software for any purpose. 
It is provided "AS IS" without any express or implied warranty, including the implied warranties of merchantability, 
fitness for a particular purpose and non-infringement. IBM shall not be liable for any direct, indirect, 
special or consequential damages resulting from the loss of use, data or projects, whether in an action of contract 
or tort, arising out of or in connection with the use or performance of this software.
}

\seealso{\code{\link{HW.tests}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,]
split(x,cod)

HW.original(x,cod)

fac <- factor(annualflows["cod"][,],levels=c(34:38))
x2 <- annualflows[!is.na(fac),"dato"]
cod2 <- annualflows[!is.na(fac),"cod"]

HW.original(x2,cod2)

plot(HW.original(x2,cod2))
}
\keyword{htest}
