# the .thresML is a minimum threshold below which the loglikelihood is considered is not considered 
# (i.e., if the calculated loglikelihood is lower to .thresML, it is set to .thresML for numerical reasons).
.thresML = -6000


BayesianMCMC <- function (xcont, xhist=NA, infhist=NA, suphist=NA, nbans=NA, seuil=NA,
                          nbpas=1000, nbchaines=3, confint=c(0.05, 0.95), dist="GEV",
                          apriori=function(...){1}, parameters0=NA, varparameters0=NA) {

 # This is the main function! other functions are called in this one and are written below.
 # To know the meaning of the inputs and the outputs you can type help(BayesianMCMC)

 reject <- round(nbpas/10)   # number of initial steps that will be excluded from the posterior distribution
 nbpas <- nbpas + reject
 returnperiods <- 10^(seq(0.1, 4, by=.1))
 nonexceedF <- 1 - 1/returnperiods

 # here, if initial parameters and parameter variances are not provided, a guess is made using the function .chooseparameters0 (see below)
 if (any(is.na(parameters0))) {
  parameters0_est <- parameters0
  parameters0 <- .chooseparameters0(xcont, xhist, infhist, suphist, dist)
  parameters0[!is.na(parameters0_est)] <- parameters0_est[!is.na(parameters0_est)]
 }
 if (any(is.na(varparameters0))) {
  varparameters0_est <- varparameters0
  varparameters0 <- (0.1*parameters0)^2
  varparameters0[!is.na(varparameters0_est)] <- varparameters0_est[!is.na(varparameters0_est)]
 }
 lpar <- length(parameters0)
 
 
 # here I initialize the arrays which will contain the results of the MCMC
 parameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))	# array 3D
 varparameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))    # array 3D
 vraisdist <- array(data=NA, dim=c(nbpas, nbchaines))
 #vraistest <- rep(NA, nbchaines)
 #nbsaut <- rep(NA, nbchaines)
 #propsaut <- rep(0, nbchaines)
 qq <- array(data=NA, dim=c(nbpas, length(nonexceedF), nbchaines))    # array 3D 

 # the following if-else selects the type of likelihood depending on the input provided
 if(all(is.na(c(xhist, infhist, suphist, seuil)))) {
  # Calcul sur les seules données systèmatiques
  funzionetest <- call(".lnvrais5", quote(parameters[1,,j]), quote(xcont), quote(dist))
  funzionecand <- call(".lnvrais5", quote(parameterscand), quote(xcont), quote(dist))
 }
 else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
  # Calcul avec info censurée mais débits historiques connus (Stedinger et Cohn, Naulet cas b)
  funzionetest <- call(".lnvrais1", quote(parameters[1,,j]), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais1", quote(parameterscand), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist)) 
 }
 else if (all(is.na(c(xhist, suphist))) & all(!is.na(c(infhist, seuil, nbans)))) {
  # Calcul avec info censurée mais débits historiques non connus (Stedinger et Cohn, Naulet cas a)
  funzionetest <- call(".lnvrais2", quote(parameters[1,,j]), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais2", quote(parameterscand), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
 }
 else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
  # Calcul avec prise en compte des seuls intervalles d'estimation de débit 
  funzionetest <- call(".lnvrais4", quote(parameters[1,,j]), quote(xcont), quote(infhist), quote(suphist), 
                                    quote(nbans), quote(seuil), quote(dist))
  funzionecand <- call(".lnvrais4", quote(parameterscand), quote(xcont), quote(infhist), quote(suphist), 
                                    quote(nbans), quote(seuil), quote(dist))
 }
 else stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): inconsistency in input data")

 # Algorithm with no repetitions
 #   # initialisation
 #   for (j in 1:nbchaines) {
 #    parameters[1,,j] <- rnorm(rep(1, lpar), mean=parameters0, sd=sqrt(varparameters0))
 #    varparameters[1,,j] <- varparameters0
 #    #vraistest[j] <- .lnvrais5(parameters[1,,j], xcont, dist)
 #    vraistest <- eval(funzionetest)
 #    vraisdist[1,j] <- vraistest
 #    qq[1,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[1,,j], dist=dist)
 #    nbsaut <- 0
 #    i=1
 #    cont=0
 #    while (i < nbpas) {
 #     cont <- cont+1
 #     parameterscand <- rnorm(rep(1, lpar), mean=parameters[i,,j], sd=sqrt(varparameters[i,,j]))
 #     #vraiscand <- .lnvrais5(parameterscand, xcont, dist)
 #     vraiscand <- eval(funzionecand)
 #     valtest <- min((exp(vraiscand - vraistest)), 1)
 #     test <- runif(1)
 #     if ((valtest > test) & (vraiscand > .thresML)) {
 #      i <- i+1
 #      nbsaut <- nbsaut + 1
 #      parameters[i,,j] <- parameterscand
 #      varparameters[i,,j] <- varparameters[i-1,,j]
 #      vraistest <- vraiscand
 #      vraisdist[i,j] <- vraistest
 #      qq[i,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[i,,j], dist=dist)
 #     }
 #     propsaut <- nbsaut/cont
 #     if (propsaut < 0.33) {
 #      varparameters[i,,j] <- varparameters[i,,j]*(1+(propsaut-0.34)/0.34/1000)
 #     }
 #     else if (propsaut > 0.35) {
 #      varparameters[i,,j] <- varparameters[i,,j]*(1+(propsaut-0.34)/0.34/1000)
 #     }
 #    }
 #   }


 # Algorithm with repetitions
 # initialisation
 propsaut <- array(data=NA, dim=c(nbpas, nbchaines))
 for (j in 1:nbchaines) {
  parameters[1,,j] <- .parameterscandMOD(parameters0, varparameters0, dist)   # forst step (see .parameterscandMOD below)
  varparameters[1,,j] <- varparameters0
  #vraistest[j] <- .lnvrais5(parameters[1,,j], xcont, dist)
  vraistest <- eval(funzionetest) + log(apriori(parameters[1,,j]))   # it is a log-likelihhod
  vraisdist[1,j] <- vraistest
  qq[1,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[1,,j], dist=dist)   # first quantiles (see .quantilesMOD below)
  nbsaut <- 0
  for (i in 2:nbpas) {
   parameterscand <- .parameterscandMOD(parameters[i-1,,j], varparameters[i-1,,j], dist)   # other steps
   vraiscand <- eval(funzionecand) + log(apriori(parameterscand))   # it is a log-likelihhod
   valtest <- min((exp(vraiscand - vraistest)), 1)
   test <- runif(1)
   #if ((valtest > test) & (vraiscand > .thresML)) {
   if (valtest > test) {   # I move to the new set of parameters
    nbsaut <- nbsaut + 1
    parameters[i,,j] <- parameterscand
    vraistest <- vraiscand
   }
   else {   # I remain where I am
    parameters[i,,j] <- parameters[i-1,,j]
   }
   vraisdist[i,j] <- vraistest
   qq[i,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[i,,j], dist=dist)   # other quantiles
   
   # acceptance rate
   propsaut[i,j] <- nbsaut/i
   if (propsaut[i,j] < 0.33) {
    varparameters[i,,j] <- varparameters[i-1,,j]*(1+(propsaut[i,j]-0.34)/0.34/1000)
   }
   else if (propsaut[i,j] > 0.35) {
    varparameters[i,,j] <- varparameters[i-1,,j]*(1+(propsaut[i,j]-0.34)/0.34/1000)
   }
   else {
    varparameters[i,,j] <- varparameters[i-1,,j]
   }
  }
 }
 nbpas <- nbpas - reject
 qq <- qq[-c(1:reject),,]
 parameters <- parameters[-c(1:reject),,]
  dimnames(parameters) <- list(c(1:nbpas), paste("par", seq(1,lpar)), paste("chain", seq(1,nbchaines)))
 varparameters <- varparameters[-c(1:reject),,]
  dimnames(varparameters) <- list(c(1:nbpas), paste("varp", seq(1,lpar)), paste("chain", seq(1,nbchaines)))
 vraisdist <- vraisdist[-c(1:reject),]
  dimnames(vraisdist) <- list(c(1:nbpas), paste("chain", seq(1,nbchaines)))
 propsaut <- propsaut[-c(1:reject),]
  dimnames(propsaut) <- list(c(1:nbpas), paste("chain", seq(1,nbchaines)))

 # The maximum likelihood is here performed simply taking the maximum of vraisdist
 dummy1 <- which.max(apply(vraisdist, 2, max))
 dummy2 <- apply(vraisdist, 2, which.max)[dummy1]
 parametersML <- parameters[dummy2,,dummy1]
 intervals <- apply(qq, 2, quantile, probs=confint, na.rm=TRUE)
 qqML <- qq[dummy2,,dummy1]
 output <- list(xcont=xcont, xhist=xhist, infhist=infhist, suphist=suphist, nbans=nbans, seuil=seuil, 
                nbpas=nbpas, nbchaines=nbchaines, dist=dist, confint=confint, apriori=apriori,
                parameters=parameters, quantiles=qq, varparameters=varparameters, 
                parameters0=parameters0, varparameters0=varparameters0,
                vraisdist=vraisdist, propsaut=propsaut,
                returnperiods=returnperiods, intervals=intervals,
                parametersML=parametersML, quantilesML=qqML)
 class(output) <- "BayesianMCMC"
 return(output)
}


# ----------------------------- #

print.BayesianMCMC <- function (x, ...) {
 dummy <- data.frame(cbind(x$quantilesML, t(x$intervals)), row.names=signif(x$returnperiods, 4))
 names(dummy)[1] <- "ML"
 print(dummy)
}


# ----------------------------- #

plot.BayesianMCMC <- function (x, which=1, ask=FALSE, ...) {
 if (ask) {
  op <- par(ask = TRUE)
  on.exit(par(op))
 }
 show <- rep(FALSE, 100)
 show[which] <- TRUE
 if (show[1]) .plotdiagnMCMC02(x, ...)
 if (show[2]) .plotdiagnMCMC01(x, ...)
 if (show[3]) .plotdiagnMCMC04(x, ...)
 if (show[4]) .plotdiagnMCMC05(x, ...)
 if (show[5]) .plotdiagnMCMC06(x, ...)   # a-priori distribution
}


# --------------------- #

.plotdiagnMCMC01 <- function(x, ...) {
 # Diagnostic plot of the parameters
 lpar <- dim(x$parameters)[2]
 #graphics.off()
 #x11()
 op <- par(mfrow=c(lpar, 3))
  for (j in 1:lpar) {
   limiti <- range(x$parameters[,j,])
   plot(x$parameters[,j,1], type="l", col=2, ylim=limiti, ylab=paste("par",j), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$parameters[,j,i], col=1+i)
   }
   abline(h=x$parametersML[j])
   #hist(x$parameters[,j,1], border=2, breaks=11, #seq(limiti[1], limiti[2], length=11), 
   #     xlab=paste("par",j), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   #for (i in 2:x$nbchaines) {
   # hist(x$parameters[,j,i], border=1+i, breaks=11, add=TRUE)
   #}
   #abline(v=x$parametersML[j])
   ht <- hist(x$parameters[,j,], plot=FALSE)
   br <- ht$breaks
   hist(x$parameters[,j,1], border=2, breaks=br, xlab=paste("par",j), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   for (i in 2:x$nbchaines) {
    hist(x$parameters[,j,i], border=1+i, breaks=br, add=TRUE)
   }
   abline(v=x$parametersML[j])
   plot(x$varparameters[,j,1], type="l", col=2, ylim=range(x$varparameters[,j,]), ylab=paste("var par",j), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$varparameters[,j,i], col=1+i)
   }
  }
 par(op)
}

.plotdiagnMCMC02 <- function(x, ...) {
 # Plot of the frequency curve
 #if(is.na(x$nbans)) nbans <- 0 else nbans <- x$nbans
 #if(is.na(x$seuil)) seuil <- Inf else seuil <- x$seuil
 T <- c(1,1000)
 X <- c(0, 1.3*max(c(x$xcont, x$xhist, x$infhist, x$suphist), na.rm=TRUE))
 plot(T, X, type="n", log="x", ...)
 grid(equilogs=FALSE)
 ret <- .pointspos3 (x$xcont, x$xhist, x$infhist, x$suphist, x$nbans, x$seuil)
 lines(x$returnperiods, x$quantilesML)
 lines(x$returnperiods, x$intervals[1,], lty=2)
 lines(x$returnperiods, x$intervals[2,], lty=2)
 invisible(ret)
}

.plotdiagnMCMC03 <- function(x, ...) {
 Nsim=10000
 #if(all(is.na(c(x$xhist, x$infhist, x$suphist, x$seuil)))) {
 # # Calcul sur les seules données systèmatiques
 # funzione <- call(".lnvrais5", quote(x$parametersML), quote(xcont), quote(dist))
 #}
 #else if (all(is.na(c(x$infhist, x$suphist))) & all(!is.na(c(x$xhist, x$seuil, x$nbans)))) {
 # # Calcul avec info censurée mais débits historiques connus (Stedinger et Cohn, Naulet cas b)
 # funzione <- call(".lnvrais1", quote(x$parametersML), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist))
 #}
 #else if (all(is.na(c(x$xhist, x$suphist))) & all(!is.na(c(x$infhist, x$seuil, x$nbans)))) {
 # # Calcul avec info censurée mais débits historiques non connus (Stedinger et Cohn, Naulet cas a)
 # funzione <- call(".lnvrais2", quote(x$parametersML), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
 #}
 #else if (all(is.na(c(x$xhist))) & all(!is.na(c(x$infhist, x$suphist, x$seuil, x$nbans)))) {
 # # Calcul avec prise en compte des seuls intervalles d'estimation de débit
 # funzione <- call(".lnvrais4", quote(x$parametersML), quote(xcont), quote(infhist), quote(suphist),
 #                                   quote(nbans), quote(seuil), quote(dist))
 #}
}

.plotdiagnMCMC04 <- function(x, ...) {
 # Plot of the acceptance rate and the likelihood
 op <- par(mfrow=c(2, 2))
  limiti <- range(x$vraisdist)
  plot(x$vraisdist[,1], type="l", col=2, ylim=limiti, ylab=paste("ln likelyhood"), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$vraisdist[,i], col=1+i)
   }
  ht <- hist(x$vraisdist, plot=FALSE)
  br <- ht$breaks
  hist(x$vraisdist[,1], border=2, breaks=br, xlab=paste("ln likelyhood"), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   for (i in 2:x$nbchaines) {
    hist(x$vraisdist[,i], border=1+i, breaks=br, add=TRUE)
   }

  limiti <- range(x$propsaut)
  plot(x$propsaut[,1], type="l", col=2, ylim=limiti, ylab=paste("acceptance rate"), xlab="")
   for (i in 2:x$nbchaines) {
    lines(x$propsaut[,i], col=1+i)
   }
  ht <- hist(x$propsaut, plot=FALSE)
  br <- ht$breaks
  hist(x$propsaut[,1], border=2, breaks=br, xlab=paste("acceptance rate"), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   for (i in 2:x$nbchaines) {
    hist(x$propsaut[,i], border=1+i, breaks=br, add=TRUE)
   }
 par(op)
}



.plotdiagnMCMC05 <- function(x, ...) {
 # A posteriori distribution of the parameters
 if (length(x$parametersML) == 2) {
  plot(x$parameters[,1:2,1], col=2, pch=".")
  for (i in 2:x$nbchaines) {
   points(x$parameters[,1:2,i], col=1+i, pch=".")
  }
 } 
 else if (length(x$parametersML) == 3) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  #layout(matrix(c(1,1,1,2,2,2,3,3,3,4,5,6), 6, 2, byrow=FALSE))
  layout(matrix(c(1,2,3,0), 2, 2, byrow=FALSE))
  plot(x$parameters[,c(1,2),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(1,2),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[1], x$parametersML[2], pch=19, cex=1.5, col=7)

  plot(x$parameters[,c(1,3),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(1,3),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[1], x$parametersML[3], pch=19, cex=1.5, col=7)
  
  plot(x$parameters[,c(3,2),1], col=2, pch=".", cex=2)
  for (i in 2:x$nbchaines) {
   points(x$parameters[,c(3,2),i], col=1+i, pch=".", cex=2)
  }
  points(x$parametersML[3], x$parametersML[2], pch=19, cex=1.5, col=7)
 
  #  par(mar=c(1,1,1,1))
  #  plot(density(x$parameters[,1,1]), col=1, main="", axes=FALSE)
  #  for (i in 2:x$nbchaines) {
  #   lines(density(x$parameters[,1,i]), col=i)
  #  }
  #  mtext("par 1", 3, -1.5, adj=0.03, cex=.8)
  #  box()

  #  plot(density(x$parameters[,2,1]), col=1, main="", axes=FALSE)
  #  for (i in 2:x$nbchaines) {
  #   lines(density(x$parameters[,2,i]), col=i)
  #  }
  #  mtext("par 2", 3, -1.5, adj=0.03, cex=.8)
  #  box()

  #  plot(density(x$parameters[,3,1]), col=1, main="", axes=FALSE)
  #  for (i in 2:x$nbchaines) {
  #   lines(density(x$parameters[,3,i]), col=i)
  #  }
  #  mtext("par 3", 3, -1.5, adj=0.03, cex=.8)
  #  box()
  par(def.par)
 }
}



.plotdiagnMCMC06 <- function(x, ...) {
 # plot the a-priori distribution of the parameters
 if (length(x$parametersML) == 2) {
  plot(x$parameters[,1:2,1], type="n", ...)
  xr <- range(x$parameters[,1,1])
  yr <- range(x$parameters[,2,1])
  xx <- seq(xr[1], xr[2], length=50)
  yy <- seq(yr[1], yr[2], length=50)
  xy <- expand.grid(xx, yy)
  zz <- rep(NA, length=dim(xy)[1])
  for (i in 1:dim(xy)[1]) {
   zz[i] <- x$apriori(as.numeric(xy[i,]))
  }
  contour(xx, yy, zz, add=TRUE, col="darkgray")
 } 
 else if (length(x$parametersML) == 3) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  x1r <- range(x$parameters[,1,1])
  x2r <- range(x$parameters[,2,1])
  x3r <- range(x$parameters[,3,1])
  xx1 <- seq(x1r[1], x1r[2], length=20)
  xx2 <- seq(x2r[1], x2r[2], length=20)
  xx3 <- seq(x3r[1], x3r[2], length=20)
  xx123 <- expand.grid(xx1, xx2, xx3)
  N <- dim(xx123)[1]
  zz <- rep(NA, length=N)
  for (i in 1:N) {
   zz[i] <- x$apriori(as.numeric(xx123[i,]))
  }
  layout(matrix(c(1,2,3,0), 2, 2, byrow=FALSE))
  plot(x$parameters[,c(1,2),1], type="n")
  contour(xx1, xx2, tapply(zz, xx123[,c(1,2)], sum), add=TRUE, col="darkgray")

  plot(x$parameters[,c(1,3),1], type="n")
  contour(xx1, xx3, tapply(zz, xx123[,c(1,3)], sum), add=TRUE, col="darkgray")
  
  plot(x$parameters[,c(3,2),1], type="n")
  contour(xx3, xx2, tapply(zz, xx123[,c(3,2)], sum), add=TRUE, col="darkgray")
  par(def.par)
 }
}





# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

.pointspos3 <- function (xcont, xhist, infhist, suphist, nbans, seuil, ...)
{
 if(all(is.na(c(xhist, infhist, suphist, seuil)))) {
  # Calcul sur les seules données systèmatiques (Cunnane plotting position)
  xcont <- sort(xcont, decreasing=TRUE)
  F <- c(1:length(xcont))
  F <- 1 - ((F - 0.4)/(length(xcont) + 1 - 2*0.4))
  T <- 1/(1 - F)
  x <- xcont
  #plot(T, x, log="x", type="n", ...)
  #grid(equilogs = FALSE)
  points(T, x)
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
  # Calcul avec info censurée mais débits historiques connus (Stedinger et Cohn, Naulet cas b)
  xseuil <- 1-((sum(xcont > seuil) + length(xhist))/(nbans + length(xcont)))
  xcont2 <- xcont[xcont <= seuil]
  soprasoglia <- sum(xcont > seuil)
  xhist2 <- c(xcont[xcont > seuil], xhist)

  xx1 <- 1 + length(xcont2) - rank(xcont2, ties.method="first")
  xx1 <- xseuil*(1 - ((xx1 - 0.4)/(length(xcont2) + 1 - 2*0.4)))
  T1 <- 1/(1 - xx1)

  xx2 <- 1 + length(xhist2) - rank(xhist2, ties.method="first")
  xx2 <- 1 - ((xx2 - 0.4)/(length(xhist2) + 1 - 2*0.4)*(1 - xseuil))
  T2 <- 1/(1 - xx2)

  x <- c(xcont2, xhist2)
  T <- c(T1, T2)
  #plot(T, x, log="x", type="n", ...)
  #grid(equilogs = FALSE)
  points(T1, xcont2)
  if (soprasoglia < 1) {
   points(T2, xhist2, pch=19)
  }
  else
  {
   points(T2[1:soprasoglia], xhist2[1:soprasoglia])
   points(T2[-c(1:soprasoglia)], xhist2[-c(1:soprasoglia)], pch=19)
  }
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(xhist, suphist))) & all(!is.na(c(infhist, seuil, nbans)))) {
  # Calcul avec info censurée mais débits historiques non connus (Stedinger et Cohn, Naulet cas a)
  xseuil <- 1-((sum(xcont > seuil) + length(infhist))/(nbans + length(xcont)))
  xcont2 <- xcont[xcont <= seuil]
  soprasoglia <- sum(xcont > seuil)
  infhist2 <- c(xcont[xcont > seuil], infhist)

  xx1 <- 1 + length(xcont2) - rank(xcont2, ties.method="first")
  xx1 <- xseuil*(1 - ((xx1 - 0.4)/(length(xcont2) + 1 - 2*0.4)))
  T1 <- 1/(1 - xx1)

  xx2 <- 1 + length(infhist2) - rank(infhist2, ties.method="first")
  xx2 <- 1 - ((xx2 - 0.4)/(length(infhist2) + 1 - 2*0.4)*(1 - xseuil))
  T2 <- 1/(1 - xx2)

  x <- c(xcont2, infhist2)
  T <- c(T1, T2)
  #plot(T, x, log="x", type="n", ...)
  #grid(equilogs = FALSE)

  points(T1, xcont2)
  if (soprasoglia < 1) {
   points(T2, infhist2, pch=24, bg=1)
   segments(T2, infhist2, T2, 10*max(x), lty=3)
  }
  else
  {
   points(T2[1:soprasoglia], infhist2[1:soprasoglia])
   points(T2[-c(1:soprasoglia)], infhist2[-c(1:soprasoglia)], pch=24, bg=1)
   segments(T2[-c(1:soprasoglia)], infhist2[-c(1:soprasoglia)], T2[-c(1:soprasoglia)], 10*max(x), lty=3)
  }
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
  # Calcul avec prise en compte des seuls intervalles d'estimation de débit
  xseuil <- 1-((sum(xcont > seuil) + length(infhist))/(nbans + length(xcont)))
  xcont2 <- xcont[xcont <= seuil]
  soprasoglia <- sum(xcont > seuil)
  infhist2 <- c(xcont[xcont > seuil], infhist)
  suphist2 <- c(xcont[xcont > seuil], suphist)
  mezzohist2 <- (suphist2 + infhist2)/2

  xx1 <- 1 + length(xcont2) - rank(xcont2, ties.method="first")
  xx1 <- xseuil*(1 - ((xx1 - 0.4)/(length(xcont2) + 1 - 2*0.4)))
  T1 <- 1/(1 - xx1)

  xx2 <- 1 + length(mezzohist2) - rank(mezzohist2, ties.method="first")
  xx2 <- 1 - ((xx2 - 0.4)/(length(mezzohist2) + 1 - 2*0.4)*(1 - xseuil))
  T2 <- 1/(1 - xx2)

  x <- c(xcont2, infhist2)
  T <- c(T1, T2)
  #plot(T, x, log="x", type="n", ...)
  #grid(equilogs = FALSE)
  points(T1, xcont2)
  if (soprasoglia < 1) {
   points(T2, infhist2, pch=24, bg=1)
   points(T2, suphist2, pch=25, bg=1)
   segments(T2, infhist2, T2, suphist2, lty=3)
  }
  else
  {
   points(T2[1:soprasoglia], infhist2[1:soprasoglia])
   points(T2[-c(1:soprasoglia)], infhist2[-c(1:soprasoglia)], pch=24, bg=1)
   points(T2[-c(1:soprasoglia)], suphist2[-c(1:soprasoglia)], pch=25, bg=1)
   segments(T2[-c(1:soprasoglia)], infhist2[-c(1:soprasoglia)], T2[-c(1:soprasoglia)], suphist2[-c(1:soprasoglia)], lty=3)
  }
  invisible(cbind(T,x))
 }
}

# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.chooseparameters0 <- function (xcont, xhist, infhist, suphist, dist) {
  campionissimo <- sort(c(xcont, xhist, infhist, suphist))   # valore sballato ma non cannatissimo
  ll <- as.numeric(Lmoments(campionissimo))
  if (dist=="GEV") {
   parameters0 <- unlist(par.GEV(ll[1], ll[2], ll[4]))
  }
  else if (dist=="NORM") {
   parameters0 <- c(mean(campionissimo), sd(campionissimo))
  }
  else if (dist=="EXP") {
   parameters0 <- unlist(par.exp(ll[1], ll[2]))
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont)
  }
  else if (dist=="GENLOGIS") {
   parameters0 <- unlist(par.genlogis(ll[1], ll[2], ll[4]))
  }
  else if (dist=="GENPAR") {
   parameters0 <- unlist(par.genpar(ll[1], ll[2], ll[4]))
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont)
  }
  else if ((dist=="GUMBEL")||(dist=="EV1")) {
   parameters0 <- unlist(par.gumb(ll[1], ll[2]))
  }
  else if (dist=="KAPPA") {
   parameters0 <- unlist(par.kappa(ll[1], ll[2], ll[4], ll[5]))
  }
  else if ((dist=="LOGNORM")||(dist=="LN3")) {
   parameters0 <- unlist(par.lognorm(ll[1], ll[2], ll[4]))
  }
  else if ((dist=="LN")||(dist=="LN2")) {
   parameters0 <- c(mean(log(campionissimo)), sd(log(campionissimo))) 
  }
  else if ((dist=="P3")||(dist=="GAM")) {
   parameters0 <- unlist(par.gamma(ll[1], ll[2], ll[4])[1:3])
   if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont) # asimmetria positiva!!!
  }
  else stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): distribution unknown")
  return(parameters0)
}


# ---------------------------- #

.parameterscandMOD <- function (parameters, varparameters, dist) {
  # Perform a step using different distributions (normal or lognormal) depending on the parameter
  # (essentially if it must be positive or can be also negative)
  if (dist=="GEV") {
   #parameterscand <- rnorm(rep(1, 3), mean=parameters, sd=sqrt(varparameters))
   # I know that the scale parameter is positive
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if (dist=="NORM") {
   #parameterscand <- rnorm(rep(1, 2), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if (dist=="EXP") {
   #parameterscand <- rnorm(rep(1, 2), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if (dist=="GENLOGIS") {
   #parameterscand <- rnorm(rep(1, 3), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if (dist=="GENPAR") {
   #parameterscand <- rnorm(rep(1, 3), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if ((dist=="GUMBEL")||(dist=="EV1")) {
   #parameterscand <- rnorm(rep(1, 2), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if (dist=="KAPPA") {
   parameterscand <- rnorm(rep(1, 4), mean=parameters, sd=sqrt(varparameters))
  }
  else if ((dist=="LOGNORM")||(dist=="LN3")) {
   #parameterscand <- rnorm(rep(1, 3), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else if ((dist=="LN")||(dist=="LN2")) {
   #parameterscand <- rnorm(rep(1, 2), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter1 <- rnorm(1, mean=parameters, sd=sqrt(varparameters))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter1, parameter2)
  }
  else if ((dist=="P3")||(dist=="GAM")) {
   #parameterscand <- rnorm(rep(1, 3), mean=parameters, sd=sqrt(varparameters))
   LNvarparameters <- log(1 + varparameters[2]/parameters[2]^2)
   LNparameters <- log(parameters[2]) - LNvarparameters/2
   parameter13 <- rnorm(rep(1,2), mean=parameters[c(1,3)], sd=sqrt(varparameters[c(1,3)]))
   parameter2 <- rlnorm(1, meanlog=LNparameters, sdlog=sqrt(LNvarparameters))
   parameterscand <- c(parameter13[1], parameter2, parameter13[2])
  }
  else stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): distribution unknown")
  return(parameterscand)
}



# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.quantilesMOD <- function (F, parameters, dist="GEV") {
 # calculates the quantiles for a given distribution, a given parameter set and a given non-exceedance probability
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]  
  qq <- invF.GEV(F, xi, alfa, k)
 }
 else if (dist=="NORM") {
  qq <- qnorm(F, mean=parameters[1], sd=parameters[2])
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  qq <- invF.exp(F, xi, alfa)
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.genlogis(F, xi, alfa, k)
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.genpar(F, xi, alfa, k)
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  qq <- invF.gumb(F, xi, alfa)
 }
 else if (dist=="KAPPA") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  h <- parameters[4]
  qq <- invF.kappa(F, xi, alfa, k, h)
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.lognorm(F, xi, alfa, k)
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  #varlogLN <- log(1 + parameters[2]/parameters[1]^2)
  #meanlogML <- log(parameters[1]) - varlogLN/2
  #qq <- qlnorm(F, meanlog=meanlogML, sdlog=sqrt(varlogLN))
  qq <- qlnorm(F, meanlog=parameters[1], sdlog=parameters[2])
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  qq <- invF.gamma(F, xi, beta, alfa)
 }
 else stop(".quantilesMOD(F, parameters, dist): distribution unknown")

 return(qq)
}

# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

.lnvrais5 <- function (parameters, xcont, dist="GEV") {
 # Calcul sur les seules données systèmatiques
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]  
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(f.GEV(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if (dist=="NORM") {
  lnvrais <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvrais <- sum(log(f.exp(xcont, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0 & all(xcont > xi)) {
   lnvrais <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvrais <- sum(log(f.gumb(xcont, xi, alfa)))
 }
 #else if (dist=="KAPPA") {
 # xi <- parameters[1]
 # alfa <- parameters[2]
 # k <- parameters[3]
 # h <- parameters[4]
 # if (sum((k*(xcont-xi)/alfa) > 1)==0) {
 #  lnvrais <- sum(log(f.kappa(xcont, xi, alfa, k, h)))
 # }
 # else lnvrais <- .thresML
 #}
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvrais <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  if (alfa > 0) {
   lnvrais <- sum(log(f.gamma(xcont, xi, beta, alfa)))
  }
  else lnvrais <- .thresML
 }
 else stop(".lnvrais5(parameters, xcont, dist): distribution unknown")
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}




# ---------------- #

.lnvrais1 <- function (parameters, xcont, xhist, nbans, seuil, dist="GEV") {
 longhist <- length(xhist)
 # Calcul avec info censurée mais débits historiques connus (Stedinger et Cohn, Naulet cas b)
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]  
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.GEV(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.GEV(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(xhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(f.GEV(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) *
                     log(pnorm(seuil, mean=parameters[1], sd=parameters[2]))
  lnvraishist <- lnvraishist + sum(log(dnorm(xhist, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.exp(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * 
                     log(F.exp(seuil, xi, alfa)) + sum(log(f.exp(xhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.genlogis(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(xhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(f.genlogis(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0 & all(xcont > xi)) {
  #if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (((k*(seuil-xi)/alfa) < 1) & (seuil > xi)) {
  #if (((k*(seuil-xi)/alfa) < 1)) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.genpar(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(xhist-xi)/alfa) > 1)==0 & all(xhist > xi)) {
  #if (sum((k*(xhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(f.genpar(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- (nbans - longhist) * log(F.gumb(seuil, xi, alfa)) + sum(log(f.gumb(xhist, xi, alfa)))
 }
 #else if (dist=="KAPPA") {
 # xi <- parameters[1]
 # alfa <- parameters[2]
 # k <- parameters[3]
 # h <- parameters[4]
 # lnvraiscont <- sum(log(f.kappa(xcont, xi, alfa, k, h)))
 # lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.kappa(seuil, xi, alfa, k, h)) + sum(log(1 - F.kappa(xhist, xi, alfa, k, h)))
 #}
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.lognorm(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(xhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(f.lognorm(xhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * 
                     log(plnorm(seuil, meanlog=parameters[1], sdlog=parameters[2]))
  lnvraishist <- lnvraishist + sum(log(dlnorm(xhist, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  if (alfa > 0) {
   lnvraiscont <- sum(log(f.gamma(xcont, xi, beta, alfa)))
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.gamma(seuil, xi, beta, alfa)) + 
                  sum(log(f.gamma(xhist, xi, beta, alfa)))
  }
  else {
   lnvraiscont <- .thresML
   lnvraishist <- .thresML
  }
 }
 else stop(".lnvrais1(parameters, xcont, xhist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}





# ---------------- #

.lnvrais2 <- function (parameters, xcont, infhist, nbans, seuil, dist="GEV") {
 longhist <- length(infhist)
 # Calcul avec info censurée mais débits historiques non connus (Stedinger et Cohn, Naulet cas a)
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]  
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.GEV(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.GEV(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(1 - F.GEV(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) *
                     log(pnorm(seuil, mean=parameters[1], sd=parameters[2]))
  lnvraishist <- lnvraishist + sum(log(1 - pnorm(infhist, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.exp(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.exp(seuil, xi, alfa)) + sum(log(1 - F.exp(infhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.genlogis(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(1 - F.genlogis(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0 & all(xcont > xi)) {
   lnvraiscont <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (((k*(seuil-xi)/alfa) < 1) & (seuil > xi)) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.genpar(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0 & all(infhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(1 - F.genpar(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.gumb(seuil, xi, alfa)) + sum(log(1 - F.gumb(infhist, xi, alfa)))
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.lognorm(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(1 - F.lognorm(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) *
                     log(plnorm(seuil, meanlog=parameters[1], sdlog=parameters[2]))
  lnvraishist <- lnvraishist + sum(log(1 - plnorm(infhist, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(f.gamma(xcont, xi, beta, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.gamma(seuil, xi, beta, alfa)) + 
                 sum(log(1 - F.gamma(infhist, xi, beta, alfa)))
 }
 else stop(".lnvrais2(parameters, xcont, infhist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}





# ---------------- #

.lnvrais4 <- function (parameters, xcont, infhist, suphist, nbans, seuil, dist="GEV") {
 longhist <- length(infhist)
 # Calcul avec prise en compte des seuls intervalles d'estimation de débit
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]  
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.GEV(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.GEV(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.GEV(suphist, xi, alfa, k) - F.GEV(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="NORM") {
  lnvraiscont <- sum(log(dnorm(xcont, mean=parameters[1], sd=parameters[2])))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) *
                     log(pnorm(seuil, mean=parameters[1], sd=parameters[2]))
  lnvraishist <- lnvraishist + sum(log(pnorm(suphist, mean=parameters[1], sd=parameters[2]) -
                                       pnorm(infhist, mean=parameters[1], sd=parameters[2])))
 }
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.exp(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.exp(seuil, xi, alfa)) + 
                 sum(log(F.exp(suphist, xi, alfa) - F.exp(infhist, xi, alfa)))
 }
 else if (dist=="GENLOGIS") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.genlogis(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.genlogis(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist<- lnvraishist + sum(log(F.genlogis(suphist, xi, alfa, k) - F.genlogis(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if (dist=="GENPAR") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0 & all(xcont > xi)) {
   lnvraiscont <- sum(log(f.genpar(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if (((k*(seuil-xi)/alfa) < 1) & (seuil > xi)) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.genpar(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0 & all(infhist > xi)) {
   lnvraishist <- lnvraishist + sum(log(F.genpar(suphist, xi, alfa, k) - F.genpar(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="GUMBEL")||(dist=="EV1")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.gumb(seuil, xi, alfa)) + 
                 sum(log(F.gumb(suphist, xi, alfa) - F.gumb(infhist, xi, alfa)))
 }
 else if ((dist=="LOGNORM")||(dist=="LN3")) {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvraiscont <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvraiscont <- .thresML
  if ((k*(seuil-xi)/alfa) < 1) {
   lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.lognorm(seuil, xi, alfa, k))
  }
  else lnvraishist <- .thresML
  if (sum((k*(infhist-xi)/alfa) > 1)==0) {
   lnvraishist <- lnvraishist + sum(log(F.lognorm(suphist, xi, alfa, k) - F.lognorm(infhist, xi, alfa, k)))
  }
  else lnvraishist <- .thresML
 }
 else if ((dist=="LN")||(dist=="LN2")) {
  lnvraiscont <- sum(log(dlnorm(xcont, meanlog=parameters[1], sdlog=parameters[2])))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) *
                     log(plnorm(seuil, meanlog=parameters[1], sdlog=parameters[2]))
  lnvraishist <- lnvraishist + sum(log(plnorm(suphist, meanlog=parameters[1], sdlog=parameters[2]) - 
                                       plnorm(infhist, meanlog=parameters[1], sdlog=parameters[2])))
 }
 else if ((dist=="P3")||(dist=="GAM")) {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  lnvraiscont <- sum(log(f.gamma(xcont, xi, beta, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.gamma(seuil, xi, beta, alfa)) + 
                 sum(log(F.gamma(suphist, xi, beta, alfa) - F.gamma(infhist, xi, beta, alfa)))
 }
 else stop(".lnvrais4(parameters, xcont, infhist, suphist, nbans, seuil, dist): distribution unknown")
 lnvrais <- lnvraiscont + lnvraishist
 if (is.nan(lnvrais)) lnvrais <- .thresML
 if (lnvrais < .thresML) lnvrais <- .thresML

 return(lnvrais)
}



#   # ------------------------------------------------------------------------------------------------------ #
#   
#   .regBayesianMCMC <- function (xcont, xhist=NA, infhist=NA, suphist=NA, nbans=NA, seuil=NA,
#                                 index=2,
#                                 nbpas=1000, nbchaines=3, confint=c(0.05, 0.95), dist="GEV",
#                                 apriori=function(...){1}, parameters0=NA, varparameters0=NA) {
#    # index = 1=mean, 2=median
#    reject <- round(nbpas/10)
#    nbpas <- nbpas + reject
#    returnperiods <- 10^(seq(0.1, 4, by=.1))
#    #returnperiods[3] <- 2
#    nonexceedF <- 1 - 1/returnperiods
#    if (any(is.na(parameters0))) {
#     parameters0_est <- parameters0
#     parameters0 <- .chooseparameters0(xcont, xhist, infhist, suphist, dist)
#     parameters0[!is.na(parameters0_est)] <- parameters0_est[!is.na(parameters0_est)]
#    }
#    if (any(is.na(varparameters0))) {
#     varparameters0_est <- varparameters0
#     varparameters0 <- (0.1*parameters0)^2
#     varparameters0[!is.na(varparameters0_est)] <- varparameters0_est[!is.na(varparameters0_est)]
#    }
#    lpar <- length(parameters0)
#    
#   
#    parameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))	# array 3D
#    varparameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))    # array 3D
#    vraisdist <- array(data=NA, dim=c(nbpas, nbchaines))
#    qq <- array(data=NA, dim=c(nbpas, length(nonexceedF), nbchaines))    # array 3D 
#   
#    if(all(is.na(c(xhist, infhist, suphist, seuil)))) {
#     # Calcul sur les seules données systèmatiques
#     funzionetest <- call(".lnvrais5", quote(parameters[1,,j]), quote(xcont), quote(dist))
#     funzionecand <- call(".lnvrais5", quote(parameterscand), quote(xcont), quote(dist))
#    }
#    else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
#     # Calcul avec info censurée mais débits historiques connus (Stedinger et Cohn, Naulet cas b)
#     funzionetest <- call(".lnvrais1", quote(parameters[1,,j]), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist))
#     funzionecand <- call(".lnvrais1", quote(parameterscand), quote(xcont), quote(xhist), quote(nbans), quote(seuil), quote(dist)) 
#    }
#    else if (all(is.na(c(xhist, suphist))) & all(!is.na(c(infhist, seuil, nbans)))) {
#     # Calcul avec info censurée mais débits historiques non connus (Stedinger et Cohn, Naulet cas a)
#     funzionetest <- call(".lnvrais2", quote(parameters[1,,j]), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
#     funzionecand <- call(".lnvrais2", quote(parameterscand), quote(xcont), quote(infhist), quote(nbans), quote(seuil), quote(dist))
#    }
#    else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
#     # Calcul avec prise en compte des seuls intervalles d'estimation de débit 
#     funzionetest <- call(".lnvrais4", quote(parameters[1,,j]), quote(xcont), quote(infhist), quote(suphist), 
#                                       quote(nbans), quote(seuil), quote(dist))
#     funzionecand <- call(".lnvrais4", quote(parameterscand), quote(xcont), quote(infhist), quote(suphist), 
#                                       quote(nbans), quote(seuil), quote(dist))
#    }
#    else stop("regBayesianMCMC(xcont, xhist, infhist, suphist, index, nbpas, nbchaines, dist): inconsistency in input data")
#   
#    # Algorithm with repetitions
#    # initialisation
#    propsaut <- array(data=NA, dim=c(nbpas, nbchaines))
#    for (j in 1:nbchaines) {
#     parameters[1,,j] <- rnorm(rep(1, lpar), mean=parameters0, sd=sqrt(varparameters0))
#     # constrain the scale parameter
#     parameters[1,,j] <- .constrainscale(parameters[1,,j], index, dist)
#     varparameters[1,,j] <- varparameters0
#     #vraistest[j] <- .lnvrais5(parameters[1,,j], xcont, dist)
#     vraistest <- eval(funzionetest) + log(apriori(parameters[1,,j]))   # it is a log-likelihhod
#     vraisdist[1,j] <- vraistest
#     qq[1,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[1,,j], dist=dist)
#     nbsaut <- 0
#     for (i in 2:nbpas) {
#      parameterscand <- rnorm(rep(1, lpar), mean=parameters[i-1,,j], sd=sqrt(varparameters[i-1,,j]))
#      # constrain the scale parameter
#      parameterscand <- .constrainscale(parameterscand, index, dist)
#      vraiscand <- eval(funzionecand) + log(apriori(parameterscand))   # it is a log-likelihhod
#      valtest <- min((exp(vraiscand - vraistest)), 1)
#      test <- runif(1)
#      #if ((valtest > test) & (vraiscand > .thresML)) {
#      if (valtest > test) {
#       nbsaut <- nbsaut + 1
#       parameters[i,,j] <- parameterscand
#       vraistest <- vraiscand
#      }
#      else {
#       parameters[i,,j] <- parameters[i-1,,j]
#      }
#      vraisdist[i,j] <- vraistest
#      qq[i,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[i,,j], dist=dist)
#      propsaut[i,j] <- nbsaut/i
#      if (propsaut[i,j] < 0.33) {
#       varparameters[i,,j] <- varparameters[i-1,,j]*(1+(propsaut[i,j]-0.34)/0.34/1000)
#      }
#      else if (propsaut[i,j] > 0.35) {
#       varparameters[i,,j] <- varparameters[i-1,,j]*(1+(propsaut[i,j]-0.34)/0.34/1000)
#      }
#      else {
#       varparameters[i,,j] <- varparameters[i-1,,j]
#      }
#     }
#    }
#    nbpas <- nbpas - reject
#    qq <- qq[-c(1:reject),,]
#    parameters <- parameters[-c(1:reject),,]
#     dimnames(parameters) <- list(c(1:nbpas), paste("par", seq(1,lpar)), paste("chain", seq(1,nbchaines)))
#    varparameters <- varparameters[-c(1:reject),,]
#     dimnames(varparameters) <- list(c(1:nbpas), paste("varp", seq(1,lpar)), paste("chain", seq(1,nbchaines)))
#    vraisdist <- vraisdist[-c(1:reject),]
#     dimnames(vraisdist) <- list(c(1:nbpas), paste("chain", seq(1,nbchaines)))
#    propsaut <- propsaut[-c(1:reject),]
#     dimnames(propsaut) <- list(c(1:nbpas), paste("chain", seq(1,nbchaines)))
#   
#    dummy1 <- which.max(apply(vraisdist, 2, max))
#    dummy2 <- apply(vraisdist, 2, which.max)[dummy1]
#    parametersML <- parameters[dummy2,,dummy1]
#    intervals <- apply(qq, 2, quantile, probs=confint, na.rm=TRUE)
#    qqML <- qq[dummy2,,dummy1]
#    output <- list(xcont=xcont, xhist=xhist, infhist=infhist, suphist=suphist, nbans=nbans, seuil=seuil, 
#                   index=index, nbpas=nbpas, nbchaines=nbchaines, confint=confint, dist=dist, apriori=apriori,
#                   parameters=parameters, quantiles=qq, varparameters=varparameters, 
#                   parameters0=parameters0, varparameters0=varparameters0,
#                   vraisdist=vraisdist, propsaut=propsaut,
#                   returnperiods=returnperiods, intervals=intervals,
#                   parametersML=parametersML, quantilesML=qqML)
#    class(output) <- "BayesianMCMC"
#    return(output)
#   }
#   
#   
#   # ------------------------------------------------ #
#   .constrainscale <- function (parameters, index=2, dist="GEV") {
#    # Calcul sur les seules données systèmatiques
#    if (dist=="GEV") {
#     xi <- parameters[1]
#     alfa <- parameters[2]
#     k <- parameters[3]  
#     if (index==1) {
#      xi <- 1 - alfa*(1 - gamma(1+k))/k
#     }
#     else if (index==2) {
#      xi <- 1 - alfa*(1 - (log(2))^k)/k
#     }
#     parameters <- c(xi, alfa, k)
#    }
#    else if (dist=="EXP") {
#    }
#    else if (dist=="GENLOGIS") {
#    }
#    else if (dist=="GENPAR") {
#    }
#    else if ((dist=="GUMBEL")||(dist=="EV1")) {
#    }
#    #else if (dist=="KAPPA") {
#    #}
#    else if ((dist=="LOGNORM")||(dist=="LN3")) {
#    }
#    else if ((dist=="P3")||(dist=="GAM")) {
#    }
#    else stop(".constrainscale (parameters, index, dist): distribution unknown")
#   
#    return(parameters)
#   }
#   
#   
