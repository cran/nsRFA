% 2005-09-23, Alberto Viglione
%
\name{DISTPLOTS}
\alias{DISTPLOTS}
\alias{plotpos}
\alias{plotposRP}
\alias{unifplot}
\alias{normplot}
\alias{lognormplot}
\alias{gumbelplot}
\alias{pointspos}
\alias{pointsposRP}
\alias{unifpoints}
\alias{normpoints}
\alias{gumbelpoints}
\alias{regionalplotpos}
\alias{regionalnormplot}
\alias{regionallognormplot}
\alias{regionalgumbelplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Empirical distribution plots}
\description{
 Sample values are plotted against their empirical distribution in graphs where points belonging to a particular distribution should lie on a straight line.
}
\usage{
 plotpos (x, a=0, orient="xF", ...)
 plotposRP (x, a=0, orient="xF", ...)
 unifplot (x, a=0, orient="xF", line=TRUE, ...)
 normplot (x, a=0, orient="xF", line=TRUE, ...)
 lognormplot (x, a=0, orient="xF", line=TRUE, ...)
 gumbelplot (x, a=0, orient="xF", line=TRUE, ...)
 pointspos (x, a=0, orient="xF", ...)
 pointsposRP (x, a=0, orient="xF", ...)
 unifpoints (x, a=0, orient="xF", ...)
 normpoints (x, a=0, orient="xF", ...)
 gumbelpoints (x, a=0, orient="xF", ...)
 regionalplotpos (x, cod, a=0, orient="xF", ...)
 regionalnormplot (x, cod, a=0, orient="xF", ...)
 regionallognormplot (x, cod, a=0, orient="xF", ...)
 regionalgumbelplot (x, cod, a=0, orient="xF", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing a data-sample}
 \item{a}{plotting position parameter, normally between 0 and 0.5 (the default value here, corresponding to the Hazen plotting position, see details)}
 \item{orient}{if \code{orient="xF"} the abscissa will be x and the ordinate F}
 \item{line}{if TRUE (default) a straight line indicating the normal, lognormal or Gumbel distribution with parameters estimated from \code{x} is plotted}
 \item{cod}{array that defines the data subdivision among sites}
 \item{...}{graphical parameters as \code{xlab}, \code{ylab}, \code{main}, ...}
}
\details{
 A brief introduction on Probability Plots (or Quantile-Quantile plots) is available on \url{http://en.wikipedia.org/wiki/Q-Q_plot}.
 For plotting positions see \url{http://en.wikipedia.org/wiki/Plotting_position}.

For the quantiles of the comparison distribution typically the Weibull formula \eqn{k/(n + 1)} is used (default here). 
Several different formulas have been used or proposed as symmetrical plotting positions. 
Such formulas have the form 
\deqn{(k - a)/(n + 1 - 2a)} 
for some value of \eqn{a} in the range from 0 to 1/2. 
The above expression \eqn{k/(n+1)} is one example of these, for \eqn{a=0}. 
The Filliben plotting position has \eqn{a = 0.3175} and the Cunanne plotting position has \eqn{a = 0.4} should be nearly quantile-unbiased for a range of distributions.
The Hazen plotting position, widely used by engineers, has \eqn{a = 0.5}.
The Blom's plotting position, \eqn{a = 3/8}, gives nearly unbiased quantiles for the normal distribution, while the Gringeton plotting position, \eqn{a = 0.44}, is optimized for the largest observations from a Gumbel distribution.
For the generalized Pareto, the GEV and related distributions of the Type I (Gumbel) and Weibull, \eqn{a = 0.35} is suggested.

For large sample size, \eqn{n}, there is little difference between these various expressions.
}
\value{
 Representation of the values of \code{x} vs their empirical probability function \eqn{F} in a cartesian, uniform, normal, lognormal or Gumbel plot. 
 \code{plotpos} and \code{unifplot} are analogous except for the axis notation, \code{unifplot} has the same notation as \code{normplot}, \code{lognormplot}, ...
 \code{plotposRP} is analogous to \code{plotpos} but the frequencies \eqn{F} are expressed as Return Periods \eqn{T=1/(1-F)}.
 With the default settings, \eqn{F} is defined with the Weibull plotting position \eqn{F=k/(n+1)}.
 The straight line (if \code{line}=TRUE) indicate the uniform, normal, lognormal or Gumbel distribution with parameters estimated from \code{x}.
 The regional plots draw samples of a region on the same plot. 

 \code{pointspos}, \code{normpoints}, ... are the analogous of \code{points}, they can be used to add points or lines to \code{plotpos}, \code{normplot}, ...
 \code{normpoints} can be used either in \code{normplot} or \code{lognormplot}.

}
%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\seealso{These functons are analogous to \code{\link{qqnorm}}; for the distributions, see \code{\link{Normal}}, \code{\link{Lognormal}}, \code{\link{LOGNORM}}, \code{\link{GUMBEL}}.}

\examples{
x <- rnorm(30,10,2)
plotpos(x)
normplot(x)
normplot(x,xlab=expression(D[m]),ylab=expression(hat(F)),
         main="Normal plot",cex.main=1,font.main=1)
normplot(x,line=FALSE)

x <- rlnorm(30,log(100),log(10))
normplot(x)
lognormplot(x)

x <- rand.gumb(30,1000,100)
normplot(x)
gumbelplot(x)

x <- rnorm(30,10,2)
y <- rnorm(50,10,3)
z <- c(x,y)
codz <- c(rep(1,30),rep(2,50))
regionalplotpos(z,codz)
regionalnormplot(z,codz,xlab="z")
regionallognormplot(z,codz)
regionalgumbelplot(z,codz)

plotpos(x)
pointspos(y,pch=2,col=2)

x <- rnorm(50,10,2)
F <- seq(0.01,0.99,by=0.01)
qq <- qnorm(F,10,2)
plotpos(x)
pointspos(qq,type="l")

normplot(x,line=FALSE)
normpoints(x,type="l",lty=2,col=3)

lognormplot(x)
normpoints(x,type="l",lty=2,col=3)

gumbelplot(x)
gumbelpoints(x,type="l",lty=2,col=3)

# distributions comparison in probabilistic graphs
x <- rnorm(50,10,2)
F <- seq(0.001,0.999,by=0.001)
loglikelhood <- function(param) {-sum(dgamma(x, shape=param[1], 
                scale=param[2], log=TRUE))}
parameters <- optim(c(1,1),loglikelhood)$par
qq <- qgamma(F,shape=parameters[1],scale=parameters[2])
plotpos(x)
pointspos(qq,type="l")

normplot(x,line=FALSE)
normpoints(qq,type="l")

lognormplot(x,line=FALSE)
normpoints(qq,type="l")

}
\keyword{hplot}
