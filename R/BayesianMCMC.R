.thresML = -6000

BayesianMCMC <- function (xcont, xhist=NA, infhist=NA, suphist=NA, nbans=NA, seuil=NA,
                          nbpas=1000, nbchaines=3, confint=c(0.05, 0.95), dist="GEV") {
 reject <- round(nbpas/10)
 nbpas <- nbpas + reject
 returnperiods <- 10^(seq(0.1, 4, by=.1))
 nonexceedF <- 1 - 1/returnperiods
 ll <- as.numeric(Lmoments(c(xcont, xhist, infhist, suphist)))	# valore sballato ma non cannatissimo
 if (dist=="GEV") {
  parameters0 <- unlist(par.GEV(ll[1], ll[2], ll[4]))
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
 else if (dist=="GUMBEL") {
  parameters0 <- unlist(par.gumb(ll[1], ll[2]))
 }
 else if (dist=="KAPPA") {
  parameters0 <- unlist(par.kappa(ll[1], ll[2], ll[4], ll[5]))
 }
 else if (dist=="LOGNORM") {
  parameters0 <- unlist(par.lognorm(ll[1], ll[2], ll[4]))
 }
 else if (dist=="P3") {
  parameters0 <- unlist(par.gamma(ll[1], ll[2], ll[4])[1:3])
  if(min(xcont) < parameters0[1]) parameters0[1] <- min(xcont)	# asimmetria positiva!!!
 }
 else stop("BayesianMCMC(xcont, xhist, infhist, suphist, nbpas, nbchaines, dist): distribution unknown")
 varparameters0 <- (0.1*parameters0)^2
 lpar <- length(parameters0)

 parameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))	# array 3D
 varparameters <- array(data=NA, dim=c(nbpas, lpar, nbchaines))    # array 3D
 vraisdist <- array(data=NA, dim=c(nbpas, nbchaines))
 #vraistest <- rep(NA, nbchaines)
 #nbsaut <- rep(NA, nbchaines)
 #propsaut <- rep(0, nbchaines)
 qq <- array(data=NA, dim=c(nbpas, length(nonexceedF), nbchaines))    # array 3D 

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

 # initialisation
 for (j in 1:nbchaines) {
  parameters[1,,j] <- rnorm(rep(1, lpar), mean=parameters0, sd=sqrt(varparameters0))
  varparameters[1,,j] <- varparameters0
  #vraistest[j] <- .lnvrais5(parameters[1,,j], xcont, dist)
  vraistest <- eval(funzionetest)
  vraisdist[1,j] <- vraistest
  qq[1,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[1,,j], dist=dist)
  nbsaut <- 0
  i=1
  cont=0
  while (i < nbpas) {
   cont <- cont+1
   parameterscand <- rnorm(rep(1, lpar), mean=parameters[i,,j], sd=sqrt(varparameters[i,,j]))
   #vraiscand <- .lnvrais5(parameterscand, xcont, dist)
   vraiscand <- eval(funzionecand)
   valtest <- min((exp(vraiscand - vraistest)), 1)
   test <- runif(1)
   if ((valtest > test) & (vraiscand > .thresML)) {
    i <- i+1
    nbsaut <- nbsaut + 1
    parameters[i,,j] <- parameterscand
    varparameters[i,,j] <- varparameters[i-1,,j]
    vraistest <- vraiscand
    vraisdist[i,j] <- vraistest
    qq[i,,j] <- .quantilesMOD(F=nonexceedF, parameters=parameters[i,,j], dist=dist)
   }
   propsaut <- nbsaut/cont
   if (propsaut < 0.33) {
    varparameters[i,,j] <- varparameters[i,,j]*(1+(propsaut-0.34)/0.34/1000)
   }
   else if (propsaut > 0.35) {
    varparameters[i,,j] <- varparameters[i,,j]*(1+(propsaut-0.34)/0.34/1000)
   }
  }
 }
 qq <- qq[-c(1:reject),,]
 parameters <- parameters[-c(1:reject),,]
 varparameters <- varparameters[-c(1:reject),,]
 vraisdist <- vraisdist[-c(1:reject),]
 nbpas <- nbpas - reject

 dummy1 <- which.max(apply(vraisdist, 2, max))
 dummy2 <- apply(vraisdist, 2, which.max)[dummy1]
 parametersML <- parameters[dummy2,,dummy1]
 intervals <- apply(qq, 2, quantile, probs=confint, na.rm=TRUE)
 qqML <- qq[dummy2,,dummy1]
 output <- list(xcont=xcont, xhist=xhist, infhist=infhist, suphist=suphist, nbans=nbans, seuil=seuil, 
                nbpas=nbpas, nbchaines=nbchaines, dist=dist,
                parameters=parameters, varparameters=varparameters, vraisdist=vraisdist, returnperiods=returnperiods, intervals=intervals,
                parametersML=parametersML, quantilesML=qqML)
 class(output) <- "BayesianMCMC"
 return(output)
}


# ----------------------------- #

print.BayesianMCMC <- function (x, ...) {
 dummy <- data.frame(cbind(x$quantilesML, t(x$intervals)), row.names=signif(x$returnperiods, 3))
 names(dummy)[1] <- "ML"
 print(dummy)
}


# ----------------------------- #

plot.BayesianMCMC <- function (x, which=c(1:2), ask=TRUE, ...) {
 if (ask) {
  op <- par(ask = TRUE)
  on.exit(par(op))
 }
 show <- rep(FALSE, 100)
 show[which] <- TRUE
 if (show[1]) .plotdiagnMCMC01(x, ...)
 if (show[2]) .plotdiagnMCMC02(x, ...)
}


# --------------------- #

.plotdiagnMCMC01 <- function(x, ...) {
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
   hist(x$parameters[,j,1], border=2, breaks=11, #seq(limiti[1], limiti[2], length=11), 
        xlab=paste("par",j), main="", xlim=limiti, ylim=c(0,x$nbpas/3))
   for (i in 2:x$nbchaines) {
    hist(x$parameters[,j,i], border=1+i, breaks=11, add=TRUE)
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
 #if(is.na(x$nbans)) nbans <- 0 else nbans <- x$nbans
 #if(is.na(x$seuil)) seuil <- Inf else seuil <- x$seuil
 plot(c(1,1000), c(0, 1.3*max(c(x$xcont, x$xhist, x$infhist, x$suphist), na.rm=TRUE)), type="n", log="x", xlab="T [yrs]", ylab="x", ...)
 grid(equilogs=FALSE)
 .pointspos3 (x$xcont, x$xhist, x$infhist, x$suphist, x$nbans, x$seuil)
 lines(x$returnperiods, x$quantilesML)
 lines(x$returnperiods, x$intervals[1,], lty=2)
 lines(x$returnperiods, x$intervals[2,], lty=2)
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
# --------------------------------------------------------------- #

.pointspos3 <- function (xcont, xhist, infhist, suphist, nbans, seuil, ...)
{
 if(all(is.na(c(xhist, infhist, suphist, seuil)))) {
  # Calcul sur les seules données systèmatiques
  xcont <- sort(xcont, decreasing=TRUE)
  F <- c(1:length(xcont))
  F <- 1 - ((F - 0.4)/(length(xcont) + 1 - 2*0.4))
  T <- 1/(1 - F)
  x <- xcont
  #plot(T, x, log="x", type="n", ...)
  #grid(equilogs = FALSE)
  points(T, x)
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

  x <- c(xcont2, mezzohist2)
  T <- c(T1, T2)
  #plot(T, x, log="x", type="n", ...)
  #grid(equilogs = FALSE)
  points(T1, xcont2)
  if (soprasoglia < 1) {
   points(T2, infhist2, pch=24, bg=1)
   points(T2, suphist2, pch=25, bg=1)
   segments(T2, infhist2, T2, suphist2, lty=2)
  }
  else
  {
   points(T2[1:soprasoglia], infhist2[1:soprasoglia])
   points(T2[-c(1:soprasoglia)], infhist2[-c(1:soprasoglia)], pch=24, bg=1)
   points(T2[-c(1:soprasoglia)], suphist2[-c(1:soprasoglia)], pch=25, bg=1)
   segments(T2[-c(1:soprasoglia)], infhist2[-c(1:soprasoglia)], T2[-c(1:soprasoglia)], suphist2[-c(1:soprasoglia)], lty=3)
  }
 }
}

# ------------------------------------------------------------------------------------------ #

.quantilesMOD <- function (F, parameters, dist="GEV") {
 if (dist=="GEV") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]  
  qq <- invF.GEV(F, xi, alfa, k)
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
 else if (dist=="GUMBEL") {
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
 else if (dist=="LOGNORM") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  qq <- invF.lognorm(F, xi, alfa, k)
 }
 else if (dist=="P3") {
  xi <- parameters[1]
  beta <- parameters[2]
  alfa <- parameters[3]
  qq <- invF.gamma(F, xi, beta, alfa)
 }
 else stop(".quantilesMOD(F, parameters, dist): distribution unknown")

 return(qq)
}

# ---------------- #

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
 else if (dist=="GUMBEL") {
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
 else if (dist=="LOGNORM") {
  xi <- parameters[1]
  alfa <- parameters[2]
  k <- parameters[3]
  if (sum((k*(xcont-xi)/alfa) > 1)==0) {
   lnvrais <- sum(log(f.lognorm(xcont, xi, alfa, k)))
  }
  else lnvrais <- .thresML
 }
 else if (dist=="P3") {
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
 else if (dist=="EXP") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.exp(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.exp(seuil, xi, alfa)) + sum(log(f.exp(xhist, xi, alfa)))
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
 else if (dist=="GUMBEL") {
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
 else if (dist=="LOGNORM") {
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
 else if (dist=="P3") {
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
 else if (dist=="GUMBEL") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.gumb(seuil, xi, alfa)) + sum(log(1 - F.gumb(infhist, xi, alfa)))
 }
 else if (dist=="LOGNORM") {
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
 else if (dist=="P3") {
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
 else if (dist=="GUMBEL") {
  xi <- parameters[1]
  alfa <- parameters[2]
  lnvraiscont <- sum(log(f.gumb(xcont, xi, alfa)))
  lnvraishist <- log(choose(nbans, longhist)) + (nbans - longhist) * log(F.gumb(seuil, xi, alfa)) + 
                 sum(log(F.gumb(suphist, xi, alfa) - F.gumb(infhist, xi, alfa)))
 }
 else if (dist=="LOGNORM") {
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
 else if (dist=="P3") {
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



# ------------------------------------------------------------------------------------------------------ #

