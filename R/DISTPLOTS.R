# ---------------------------------------------------------------- #

plotpos <- function(x,...) {

  # INPUT
  # x = colonna
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta cartesiana della serie

  ordinato <- sort(x)
  n <- length(ordinato)
  plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=plotpos
  plot(x,F,...)
  grid()

}


# ---------------------------------------------------------------- #

pointspos <- function(x,...) {

  # INPUT
  # x = colonna
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta cartesiana della serie

  ordinato <- sort(x)
  n <- length(ordinato)
  plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)

  x=ordinato
  F=plotpos
  points(x,F,...)
}


# ---------------------------------------------------------------- #

normplot <- function(x,line=TRUE,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica normale della serie

  ordinato <- sort(x)
  m <- mean(ordinato)
  s <- sd(ordinato)
  u <- (ordinato-m)/s
  n <- length(ordinato)
  plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
  uteorica <- qnorm(plotpos)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  plot(x,F,axes=FALSE,...)
  axis(1)
  axis(2,at=qq,labels=pp)
  abline(h=qq, col="gray", lty="dotted")
  grid(nx=NULL, ny=NA)
  box()
  if(line==TRUE) lines(ordinato,u,lty=2)

}


# ---------------------------------------------------------------- #

normpoints <- function(x,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica normale della serie

  ordinato <- sort(x)
  m <- mean(ordinato)
  s <- sd(ordinato)
  u <- (ordinato-m)/s
  n <- length(ordinato)
  plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
  uteorica <- qnorm(plotpos)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)

  x=ordinato
  F=uteorica
  points(x,F,...)
}


# ---------------------------------------------------------------- #

lognormplot <- function(x,line=TRUE,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica lognormale della serie

  ordinato <- sort(x)
  logordinato <- log(ordinato)
  m <- mean(logordinato)
  s <- sd(logordinato)
  u <- (logordinato-m)/s
  n <- length(ordinato)
  plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
  uteorica <- qnorm(plotpos)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  plot(x,F,log="x",axes=FALSE,...)
  axis(1)
  axis(2,at=qq,labels=pp)
  abline(h=qq, col="gray", lty="dotted")
  grid(nx=NULL, ny=NA)
  box()
  if (line==TRUE) lines(ordinato,u,lty=2)

}


# ---------------------------------------------------------------- #

gumbelplot <- function(x,line=TRUE,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica di Gumbel della serie

  ordinato <- sort(x)
  m <- mean(ordinato)
  s <- sd(ordinato)
  beta <- (s*sqrt(6))/pi
  mu <- m - 0.5772*beta
  adim <- (ordinato-mu)/beta
  n <- length(ordinato)
  plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
  uteorica <- -log(-log(plotpos))
  ordinatoteorico <- s*uteorica + m
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- -log(-log(pp))

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  plot(x,F,axes=FALSE,...)
  axis(1)
  axis(2,at=qq,labels=pp)
  abline(h=qq, col="gray", lty="dotted")
  #verticalgrid <- invF.gumb (c(.001,.01,.05,.1,.25,.50,.75,.9,.95,.99,.999),mu,beta)
  #abline(v=verticalgrid, col="gray", lty="dotted")
  grid(nx=NULL, ny=NA)
  box()
  if (line==TRUE) lines(ordinato,adim,lty=2)

}


# ---------------------------------------------------------------- #

gumbelpoints <- function(x,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica di Gumbel della serie

  ordinato <- sort(x)
  m <- mean(ordinato)
  s <- sd(ordinato)
  beta <- (s*sqrt(6))/pi
  mu <- m - 0.5772*beta
  adim <- (ordinato-mu)/beta
  n <- length(ordinato)
  plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
  uteorica <- -log(-log(plotpos))
  ordinatoteorico <- s*uteorica + m
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- -log(-log(pp))

  x=ordinato
  F=uteorica
  points(x,F,...)
}


# ---------------------------------------------------------------- #

unifplot <- function (x, line=TRUE, ...) {
    ordinato <- sort(x)
    n <- length(ordinato)
    plotpos <- seq(1,n)/n
    uteorica <- qunif(plotpos)
    pp <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
    qq <- qunif(pp)
    x = ordinato
    F = uteorica
    plot(x, F, axes = FALSE, ...)
    axis(1)
    axis(2, at = qq, labels = pp)
    abline(h = qq, col = "gray", lty = "dotted")
    #grid(nx = NULL, ny = NA)
    box()
    if (line == TRUE)
        lines(ordinato, qunif(ordinato), lty = 2)
}


# ---------------------------------------------------------------- #

unifpoints <- function (x, ...) {
    ordinato <- sort(x)
    n <- length(ordinato)
    plotpos <- seq(1,n)/n
    uteorica <- qunif(plotpos)
    x = ordinato
    F = uteorica
    points(x, F, ...)
}


# ---------------------------------------------------------------- #

regionalplotpos <- function(x,cod,...) {

  # INPUT
  # x = colonna
  # cod = .
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta cartesiana della serie

  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)

  x=c(min(x),max(x))
  F=c(0,1)
  plot(x,F,type="n",...)
  for(i in 1:k) {
   ordinato <- sort(X[cod==liv[i]])
   n <- length(ordinato)
   plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
   points(ordinato,plotpos,pch=i,col=i)
  }
  grid()

}


# --------------------------------------------------------------------------- #

regionalnormplot <- function(x,cod,...) {

  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(qnorm(0.5/(max(ni))),qnorm((max(ni)-0.5)/max(ni)))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)
  plot(x,F,type="n",axes=FALSE,...)
  axis(1)
  axis(2,at=qq,labels=pp)
  for(i in 1:k) {
   ordinato <- sort(X[cod==liv[i]])
   n <- length(ordinato)
   plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
   uteorica <- qnorm(plotpos)
   points(ordinato,uteorica,pch=i,col=i)
  }
  box()
  abline(h=qq, col="gray", lty="dotted")
  grid(nx=NULL, ny=NA)
}


# --------------------------------------------------------------------------- #

regionallognormplot <- function(x,cod,...) {

  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(qnorm(0.5/(max(ni))),qnorm((max(ni)-0.5)/max(ni)))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)
  plot(x,F,log="x",type="n",axes=FALSE,...)
  axis(1)
  axis(2,at=qq,labels=pp)
  for(i in 1:k) {
   ordinato <- sort(X[cod==liv[i]])
   n <- length(ordinato)
   plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
   uteorica <- qnorm(plotpos)
   points(ordinato,uteorica,pch=i,col=i)
  }
  box()
  abline(h=qq, col="gray", lty="dotted")
  grid(nx=NULL, ny=NA)
}


# ---------------------------------------------------------------------- #

regionalgumbelplot <- function(x,cod,...) {

  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(-log(-log(0.5/(max(ni)))),-log(-log((max(ni)-0.5)/max(ni))))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- -log(-log(pp))
  plot(x,F,type="n",axes=FALSE,...)
  axis(1)
  axis(2,at=qq,labels=pp)
  for(i in 1:k) {
   ordinato <- sort(X[cod==liv[i]])
   n <- length(ordinato)
   plotpos <- seq(0.5/n, (n-0.5)/n, by=1/n)
   uteorica <- -log(-log(plotpos))
   points(ordinato,uteorica,pch=i,col=i)
  }
  box()
  abline(h=qq, col="gray", lty="dotted")
  grid(nx=NULL, ny=NA)
}

