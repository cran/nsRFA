

# ----------------------------------------------------------------------------- #

traceWminim <- function (X,centers) {

  # B.Everitt [1974] Cluster Analysis, pages 42-43

  # INPUT
  # X        A numeric matrix of data, or an object that can be coerced to
  #          such a matrix (such as a numeric vector or a data frame with
  #          all numeric columns).
  # centers  The number of clusters.

  n <- dim(X)[1]
  k <- dim(X)[2]

  #X.norm <- (X - matrix(mean(X),nrow=n,ncol=k,byrow=TRUE))/matrix(sd(X),nrow=n,ncol=k,byrow=TRUE)

  d <- dist(X, method = "euclidean")
  tree <- hclust(d, method = "ward")

  clusters <- cutree(tree,centers)

  fine=FALSE
  cont=1
  while(fine == FALSE) {
    traceW.1 <- sumtraceW(clusters,X)
    scambi <- nearest(clusters,X)
    traceW.2 <- rep(NA,centers)
    for (i in 1:centers) {
      clusters.mod <- clusters; clusters.mod[scambi[i]] <- i
      traceW.2[i] <- sumtraceW(clusters.mod,X)
    }
    min.traceW.2 <- min(traceW.2); pos.min.traceW.2 <- which(traceW.2==min(traceW.2))
    if(min.traceW.2 > traceW.1) {
      fine=TRUE
    }
    else {
      clusters[scambi[pos.min.traceW.2]] <- pos.min.traceW.2
    }
    cont <- cont+1
  }


  return(clusters)

}


# ------------------------------------------------------------------------------- #

sumtraceW <- function (clusters,X) {

  # INPUT
  # X         A numeric matrix of data, or an object that can be coerced to
  #           such a matrix (such as a numeric vector or a data frame with
  #           all numeric columns).
  # clusters  A numeric vector containing the subdivision of X in clusters

  #clusters <- round(as.numeric(clusters))
  k <- max(clusters)
  traceW <- rep(NA,k)
  for (i in 1:k) {
    gruppo <- X[clusters==i,]
    mean.vector <- apply(gruppo,2,mean)
    distanze <- dist(rbind(mean.vector,gruppo), method = "euclidean")
    distanze <- as.matrix(distanze)[,1][-1]
    W <- crossprod(t(distanze))
    traceW[i] <- sum(diag(W))
  }

  sumtraceW <- sum(traceW)

  return(sumtraceW)
}


# ------------------------------------------------------------------------------- #

nearest <- function (clusters,X) {

  # INPUT
  # X         A numeric matrix of data, or an object that can be coerced to
  #           such a matrix (such as a numeric vector or a data frame with
  #           all numeric columns).
  # clusters  A numeric vector containing the subdivision of X in clusters

  k <- max(clusters)
  near <- rep(NA,k)
  for (i in 1:k) {
    gruppo <- X[which(clusters==i),]
    altri <- X[-which(clusters==i),]
    mean.vector <- apply(gruppo,2,mean)
    distanze <- dist(rbind(mean.vector,altri), method = "euclidean")
    distanze <- as.matrix(distanze)[,1][-1]
    near[i] <- names(which(distanze==min(distanze)))
  }

  return(near)
}


# --------------------------------------------------------------------------------------- #

AD.dist <- function (x,cod,index=2) {

  if (length(x)!=length(cod)) {stop('x and cod must have the same length')}

  ni <- tapply(x,cod,length)
  k <- nlevels(as.factor(cod))
  livelli <- levels(as.factor(cod))
  N <- sum(ni)

  if(index==1) {
    indexflood <- function(x) {m <- mean(x); return(m)}
  }
  else if(index==2) {
    indexflood <- function(x) {m <- median(x); return(m)}
  }
  med <- tapply(x,cod,indexflood)
  x.adim <- x/unsplit(med,cod)

  matrice <- matrix(NA,ncol=k,nrow=k)
  diag(matrice) <- 0

  for (i in 1:(k-1)) {
   for (j in (i+1):k) {
    fac <- factor(cod,levels=livelli[c(i,j)])
    vettore <- x.adim[!is.na(fac)]
    codij <- cod[!is.na(fac)]
    dist <- ksampleA2(vettore,codij)
    matrice[i,j] <- dist
    matrice[j,i] <- dist
   }
  }
  matrice.d <- as.dist(matrice)
  return(matrice.d)
}


# -------------------------------------------------------------------- #

roi.hom <- function (p.ungauged,p.gauged,cod.p,x,cod,test="HW",limit=2,Nsim=500,index=2) {

 parametri <- rbind(p.ungauged,p.gauged)
 n <- dim(parametri)[1]
 k <- dim(parametri)[2]
 matrice.distanze <- dist(parametri)
 distanze <- as.matrix(matrice.distanze)[-1,1]
 names(distanze) <- cod.p
 distanze.ordinate <- sort(distanze)

 i=20
 t=FALSE
 if (test=="HW") {
  while((t==FALSE)&&(i>1)) {
   bacini <- names(distanze.ordinate)[1:i]
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   H1 <- HW.tests(x.reg,cod.reg,Nsim)[1]
   if (H1<=limit) t=TRUE
   i=i-1
  }
 }
 else if (test=="AD") {
  while((t==FALSE)&&(i>1)) {
   bacini <- names(distanze.ordinate)[1:i]
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   P.AD <- ADbootstrap.test(x.reg,cod.reg,Nsim,index)[2]
   if (P.AD<=limit) t=TRUE
   i=i-1
  }
 }
 if (i <= 1) regione <- NULL
 else regione <- bacini
 return(regione)
}


# -------------------------------------------------------------------- #

roi.st.year <- function (p.ungauged,p.gauged,cod.p,x,cod,test="HW",station.year=500,Nsim=500,index=2) {

 parametri <- rbind(p.ungauged,p.gauged)
 n <- dim(parametri)[1]
 k <- dim(parametri)[2]
 matrice.distanze <- dist(parametri)
 distanze <- as.matrix(matrice.distanze)[-1,1]
 names(distanze) <- cod.p
 distanze.ordinate <- sort(distanze)

 ni <- tapply(x,cod,length)
 sum.ni <- cumsum(ni[names(distanze.ordinate)])
 bacini <- names(sum.ni)[1:(sum(sum.ni<station.year)+1)]

 if (test=="HW") {
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   H.HW <- HW.tests(x.reg,cod.reg,Nsim)
   homtest <- H.HW
 }
 else if (test=="AD") {
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   P.AD <- ADbootstrap.test(x.reg,cod.reg,Nsim,index)
   homtest <- P.AD
 }

 roiout <- list(region=bacini,test=homtest)

 return(roiout)
}
