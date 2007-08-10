% 2005-09-23, Alberto Viglione
%
\name{DISTPLOTS}
\alias{DISTPLOTS}
\alias{plotpos}
\alias{unifplot}
\alias{normplot}
\alias{lognormplot}
\alias{gumbelplot}
\alias{pointspos}
\alias{unifpoints}
\alias{normpoints}
\alias{gumbelpoints}
\alias{regionalplotpos}
\alias{regionalnormplot}
\alias{regionallognormplot}
\alias{regionalgumbelplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Distribution plots}
\description{
 Sample values are plotted against their empirical distribution in graphs where points belonging to a particular distribution should lie on a straight line.
}
\usage{
 plotpos (x, ...)
 unifplot (x, line=TRUE, ...)
 normplot (x, line=TRUE, ...)
 lognormplot (x, line=TRUE, ...)
 gumbelplot (x, line=TRUE, ...)
 pointspos (x, ...)
 unifpoints (x, ...)
 normpoints (x, ...)
 gumbelpoints (x, ...)
 regionalplotpos (x, cod, ...)
 regionalnormplot (x, cod, ...)
 regionallognormplot (x, cod, ...)
 regionalgumbelplot (x, cod, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{vector representing a data-sample}
 \item{line}{if TRUE (default) a straigth line indicating the normal, lognormal or Gumbel distribution with parameters estimated from \code{x} is plotted}
 \item{cod}{array that defines the data subdivision among sites}
 \item{...}{graphical parameters as \code{xlab}, \code{ylab}, \code{main}, ...}
}
\details{
 A brief introduction on Probability Plots (or Quantile-Quantile plots) is available on \url{http://en.wikipedia.org/wiki/Q-Q_plot}.
}
\value{
 Representation of the values of \code{x} vs their empirical probability function \eqn{F} in a cartesian, uniform, normal, lognormal or Gumbel plot. 
 \code{plotpos} and \code{unifplot} are analogous except for the axis notation, \code{unifplot} has the same notation as \code{normplot}, \code{lognormplot}, ...
 \eqn{F} is defined with the plotting position \eqn{F=(n-0.5)/n}.
 The straigth line (if \code{line}=TRUE) indicate the uniform, normal, lognormal or Gumbel distribution with parameters estimated from \code{x}.
 The regional plots draw samples of a region on the same plot. 

 \code{pointspos}, \code{normpoints}, ... are the analogous of \code{points}, they can be used to add points or lines to \code{plotpos}, \code{normplot}, ...
 \code{normpoints} can be used either in \code{normplot} or \code{lognormplot}.

}
%\references{}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{}

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
