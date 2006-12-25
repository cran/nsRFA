exampleRFA01 <- function () {

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
         REGIONAL FREQUENCY ANALYSIS OF ANNUAL FLOWS IN PIEMONTE AND VALLE D'AOSTA
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n
  # ------------------------------------------------------------------------------------- #
                           Regionalization of the index-flow
  # ------------------------------------------------------------------------------------- #
  \n")

  readline("\n
  Press Return to continue
  \n")

  cat("\n
  Data loading: \n
  > data(hydroSIMN)
  \n")

  data(hydroSIMN)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  > str(parameters)
  \n")

  str(parameters)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Catchment parameters plot: \n
  > plot(parameters)
  \n")

  plot(parameters)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Catchment parameters plot: \n
  > par(mfrow=c(2,2))
  > plot(parameters[c(\"Xbar\",\"Ybar\")])
  > plot(parameters[c(\"Hm\",\"Ybar\")])
  > plot(parameters[c(\"Xbar\",\"S\")])
  > plot(parameters[c(\"Hm\",\"S\")])
  > par(mfrow=c(1,1))
  \n")

  par(mfrow=c(2,2))
  plot(parameters[c("Xbar","Ybar")])
  plot(parameters[c("Hm","Ybar")])
  plot(parameters[c("Xbar","S")])
  plot(parameters[c("Hm","S")])
  par(mfrow=c(1,1))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Choice of the best regressions:
   create a function using the 'leaps' function of package 'subselect'
   to perform all-possible-regressions: \n
  > bestregressions <- function(dip,ind) {
  +  Y <- as.numeric(dip)
  +  X <- ind
  +  Sy <- var(Y)
  +  Sx <- var(X)
  +  Sxy <- var(X,Y)
  +  Dm.mat <- Sx
  +  Dm.H <- Sxy %*% t(Sxy)/Sy
  +  require(subselect)
  +  Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=3, H=Dm.H, r=1, nsol=3)
  +  Dm.leaps
  +  for(i in 3:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:3),i]])}}
  + }

  > dataregr <- cbind(parameters,(parameters[\"Dm\"])^(1/3),log(parameters[\"Dm\"]),log(parameters[\"Am\"]))
  > names(dataregr) <- c(names(parameters),\"sq3.Dm\",\"ln.Dm\",\"ln.Am\")

  Regressions between Dm and parameters: \n
  > bestregressions(dataregr[,2],dataregr[,-c(1,2,17,18)])
[1] \"Am\"    \"S2000\" \"IT\"
[1] \"S2000\" \"IT\"    \"ln.Am\"
[1] \"Am\"    \"S2000\" \"IB\"
[1] \"S2000\" \"ln.Am\"
[1] \"Am\"    \"S2000\"
[1] \"Hm\"    \"ln.Am\"
[1] \"IT\"
[1] \"IB\"
[1] \"ln.Am\"
  \n")

  dataregr <- cbind(parameters,(parameters["Dm"])^(1/3),log(parameters["Dm"]),log(parameters["Am"]))
  names(dataregr) <- c(names(parameters),"sq3.Dm","ln.Dm","ln.Am")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Regressions between sq3.Dm and parameters: \n
  > bestregressions(dataregr[,17],dataregr[,-c(1,2,17,18)])
[1] \"Hm\"   \"NORD\" \"IB\"
[1] \"S2000\" \"IT\"    \"ln.Am\"
[1] \"Hm\"    \"S2000\" \"ln.Am\"
[1] \"S2000\" \"ln.Am\"
[1] \"Hm\"    \"ln.Am\"
[1] \"Hm\" \"IB\"
[1] \"IT\"
[1] \"IB\"
[1] \"ln.Am\"
  \n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Regressions between ln.Dm and parameters: \n
  > bestregressions(dataregr[,18],dataregr[,-c(1,2,17,18)])
[1] \"Hm\"   \"NORD\" \"IB\"
[1] \"Hm\"    \"NORD\"  \"ln.Am\"
[1] \"Hm\"    \"EST\"   \"ln.Am\"
[1] \"S2000\" \"ln.Am\"
[1] \"Hm\"    \"ln.Am\"
[1] \"Hm\" \"IB\"
[1] \"IT\"
[1] \"IB\"
[1] \"ln.Am\"
  \n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Best regression: \n
  > regr01 <- lm(ln.Dm ~ 1 + Hm + NORD + IB, dataregr)
  > summary(regr01)
  \n")

  regr01 <- lm(ln.Dm ~ 1 + Hm + NORD + IB, dataregr)
  print(summary(regr01))

  cat("\n
  > R2.lm(regr01)
  \n")

  print(R2.lm(regr01))

  cat("\n
  > prt.lm(regr01)
  \n")

  print(prt.lm(regr01))

  cat("\n
  > vif.lm(regr01)
  \n")

  print(vif.lm(regr01))

  cat("\n
  RMSE: \n
  > resid <- exp(regr01$fitted.values) - parameters[,c(\"Dm\")]
  > sqrt(sum((resid)^2)/length(resid))
  \n")

  resid <- exp(regr01$fitted.values) - parameters[,c("Dm")]
  print(sqrt(sum((resid)^2)/length(resid)))

  cat("\n
  RMSEcv: \n
  > predictions <- jackknife1.lm(regr01)
  > resid <- exp(predictions) - parameters[,c(\"Dm\")]
  > sqrt(sum((resid)^2)/length(resid))
  \n")

  predictions <- jackknife1.lm(regr01)
  resid <- exp(predictions) - parameters[,c("Dm")]
  print(sqrt(sum((resid)^2)/length(resid)))

  cat("\n
  > round(cor(regr01$model[-1]),3)
  \n")

  print(round(cor(regr01$model[-1]),3))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Best regression plots: \n
  > par(mfrow=c(2,2))
  > plot(regr01$fitted.values,regr01$residuals,xlab=\"Fitted\",ylab=\"Residuals\")
  > abline(0,0,lty=3)
  > normplot(regr01$residuals,xlab=\"Residuals\")
  > plot(parameters[,c(\"Dm\")],exp(regr01$fitted.values),xlab=\"Originals\",ylab=\"Fitted\")
  > abline(0,1,lty=3)
  > intervals <- predinterval.lm(regr01)
  > intervals <- intervals[order(intervals[,1]),]
  > plot(parameters[,c(\"Dm\")],exp(predictions),xlab=\"Originals\",ylab=\"Predicted\")
  > abline(0,1,lty=3)
  > lines(exp(intervals[,c(1,2)]),lty=2)
  > lines(exp(intervals[,c(1,3)]),lty=2)
  > par(mfrow=c(1,1))
  \n")

  par(mfrow=c(2,2))
  plot(regr01$fitted.values,regr01$residuals,xlab="Fitted",ylab="Residuals")
  abline(0,0,lty=3)
  normplot(regr01$residuals,xlab="Residuals")
  plot(parameters[,c("Dm")],exp(regr01$fitted.values),xlab="Originals",ylab="Fitted")
  abline(0,1,lty=3)
  intervals <- predinterval.lm(regr01)
  intervals <- intervals[order(intervals[,1]),]
  plot(parameters[,c("Dm")],exp(predictions),xlab="Originals",ylab="Predicted")
  abline(0,1,lty=3)
  lines(exp(intervals[,c(1,2)]),lty=2)
  lines(exp(intervals[,c(1,3)]),lty=2)
  par(mfrow=c(1,1))
  
  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Simplest regression: \n
  > regr02 <- lm(sq3.Dm ~ 1 + ln.Am + Hm, dataregr)
  > summary(regr02)
  \n")

  regr02 <- lm(sq3.Dm ~ 1 + ln.Am + Hm, dataregr)
  print(summary(regr02))

  cat("\n
  > R2.lm(regr02)
  \n")

  print(R2.lm(regr02))

  cat("\n
  > prt.lm(regr02)
  \n")

  print(prt.lm(regr02))

  cat("\n
  > vif.lm(regr02)
  \n")

  print(vif.lm(regr02))

  cat("\n
  RMSE: \n
  > resid <- regr02$fitted.values^3 - parameters[,c(\"Dm\")]
  > sqrt(sum((resid)^2)/length(resid))
  \n")

  resid <- regr02$fitted.values^3 - parameters[,c("Dm")]
  print(sqrt(sum((resid)^2)/length(resid)))

  cat("\n
  RMSEcv: \n
  > predictions <- jackknife1.lm(regr02)
  > resid <- predictions^3 - parameters[,c(\"Dm\")]
  > sqrt(sum((resid)^2)/length(resid))
  \n")

  predictions <- jackknife1.lm(regr02)
  resid <- predictions^3 - parameters[,c("Dm")]
  print(sqrt(sum((resid)^2)/length(resid)))

  cat("\n
  > round(cor(regr02$model[-1]),3)
  \n")

  print(round(cor(regr02$model[-1]),3))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Best regression plots: \n
  > par(mfrow=c(2,2))
  > plot(regr02$fitted.values,regr02$residuals,xlab=\"Fitted\",ylab=\"Residuals\")
  > abline(0,0,lty=3)
  > normplot(regr02$residuals,xlab=\"Residuals\")
  > plot(parameters[,c(\"Dm\")],(regr02$fitted.values)^3,xlab=\"Originals\",ylab=\"Fitted\")
  > abline(0,1,lty=3)
  > intervals <- predinterval.lm(regr02)
  > intervals <- intervals[order(intervals[,1]),]
  > plot(parameters[,c(\"Dm\")],predictions^3,xlab=\"Originals\",ylab=\"Predicted\")
  > abline(0,1,lty=3)
  > lines((intervals[,c(1,2)])^3,lty=2)
  > lines((intervals[,c(1,3)])^3,lty=2)
  > par(mfrow=c(1,1))
  \n")

  par(mfrow=c(2,2))
  plot(regr02$fitted.values,regr02$residuals,xlab="Fitted",ylab="Residuals")
  abline(0,0,lty=3)
  normplot(regr02$residuals,xlab="Residuals")
  plot(parameters[,c("Dm")],(regr02$fitted.values)^3,xlab="Originals",ylab="Fitted")
  abline(0,1,lty=3)
  intervals <- predinterval.lm(regr02)
  intervals <- intervals[order(intervals[,1]),]
  plot(parameters[,c("Dm")],predictions^3,xlab="Originals",ylab="Predicted")
  abline(0,1,lty=3)
  lines((intervals[,c(1,2)])^3,lty=2)
  lines((intervals[,c(1,3)])^3,lty=2)
  par(mfrow=c(1,1))


  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                                          THE END
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")
}



# ------------------------------------------------------------------------------------------- #


exampleRFA02 <- function () {

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
         REGIONAL FREQUENCY ANALYSIS OF ANNUAL FLOWS IN PIEMONTE AND VALLE D'AOSTA
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n
  # ------------------------------------------------------------------------------------- #
                           Regionalization of the growth-curve
  # ------------------------------------------------------------------------------------- #
  \n")

  readline("\n
  Press Return to continue
  \n")

  cat("\n
  Data loading: \n
  > data(hydroSIMN)
  \n")

  data(hydroSIMN)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  > str(annualflows)
  \n")

  str(annualflows)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Plot consistency of data series: \n
  > y <- annualflows[\"anno\"][,]
  > cod <- annualflows[\"cod\"][,]
  > consistencyplot(y,cod) 
  \n")


  y <- annualflows["anno"][,]
  cod <- annualflows["cod"][,]
  consistencyplot(y,cod)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  L-moments of data series: \n
  > D <- annualflows[\"dato\"][,]
  > lmSIMN <- data.frame(t(sapply(split(D,cod),Lmoments)))
  > lmSIMN
  \n")
  
  D <- annualflows["dato"][,]
  lmSIMN <- data.frame(t(sapply(split(D,cod),Lmoments)))
  print(lmSIMN)
  
  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  L-moments ratio plots: \n
  > plot(lmSIMN[3:5])
  \n")

  plot(lmSIMN[3:5])

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Choice of sites with more than 15 records: \n
  > Dlen <- tapply(D,cod,length)
  > annualflows15 <- annualflows[unsplit(Dlen,cod)>=15,]
  > parameters15 <- parameters[Dlen>=15,]
  > D15 <- annualflows15[\"dato\"][,]
  > cod15 <- annualflows15[\"cod\"][,]
  \n")

  Dlen <- tapply(D,cod,length)
  annualflows15 <- annualflows[unsplit(Dlen,cod)>=15,]
  parameters15 <- parameters[Dlen>=15,]
  D15 <- annualflows15["dato"][,]
  cod15 <- annualflows15["cod"][,]

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  L-moments ratio plot: \n
  > D15lcv <- tapply(D15,cod15,LCV)
  > D15lca <- tapply(D15,cod15,LCA)
  > plot(D15lca,D15lcv,xlab=\"L-CA\",ylab=\"L-CV\",pch=1,col=1,main=\"L-moments space\",cex.main=1,font.main=1)
  > grid()
  \n")

  D15lcv <- tapply(D15,cod15,LCV)
  D15lca <- tapply(D15,cod15,LCA)
  plot(D15lca,D15lcv,xlab="L-CA",ylab="L-CV",pch=1,col=1,main="L-moments space",cex.main=1,font.main=1)
  grid()

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Choice of the classification variables:
   create a function using the 'leaps' function of package 'subselect' 
   to perform all-possible-regressions: \n
  > bestregressions <- function(dip,ind) {
  +  Y <- as.numeric(dip)
  +  X <- ind
  +  Sy <- var(Y)
  +  Sx <- var(X)
  +  Sxy <- var(X,Y)
  +  Dm.mat <- Sx
  +  Dm.H <- Sxy %*% t(Sxy)/Sy
  +  require(subselect)
  +  Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=3, H=Dm.H, r=1, nsol=3)
  +  Dm.leaps
  +  for(i in 3:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:3),i]])}}
  + }

  > bestregressions(as.numeric(D15lcv),parameters15[,3:16])
[1] \"S2000\" \"Rc\"    \"IB\"
[1] \"Am\"    \"S2000\" \"Rc\"
[1] \"Am\"    \"LLDP\"  \"S2000\"
[1] \"S2000\" \"Ybar\"
[1] \"Hm\"   \"Ybar\"
[1] \"Am\"    \"S2000\"
[1] \"S2000\"
[1] \"Hm\"
[1] \"Ybar\"
  \n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Choice of the classification variables:
   or reasoning with distance matrices:
  > bestregressions(as.numeric(AD.dist(D15,cod15)),data.frame(apply(parameters15[,3:16],2,dist)))
[1] \"Am\"    \"S2000\" \"Ybar\"
[1] \"S2000\" \"EST\"   \"Ybar\"
[1] \"Pm\"    \"S2000\" \"Ybar\"
[1] \"S2000\" \"Ybar\"
[1] \"Hm\"   \"Ybar\"
[1] \"S2000\" \"EST\"
[1] \"S2000\"
[1] \"Hm\"
[1] \"Ybar\"
  \n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  We choose Hm and Ybar as classification variables.
  Mantel test:
  > Y <- AD.dist(D15,cod15)
  > X <- data.frame(apply(parameters15[,c(\"Hm\",\"Ybar\")],2,dist))
  > datamantel <- cbind(as.numeric(Y),X)
  > regrmantel <- lm(Y ~ Hm + Ybar, datamantel)
  > summary(regrmantel)
  \n")

  Y <- AD.dist(D15,cod15)
  X <- data.frame(apply(parameters15[,c("Hm","Ybar")],2,dist))
  datamantel <- cbind(as.numeric(Y),X)
  regrmantel <- lm(Y ~ Hm + Ybar, datamantel)
  print(summary(regrmantel))

  cat("\n
  > mantel.lm(regrmantel, Nperm=100)
  \n")

  print(mantel.lm(regrmantel, Nperm=100))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  We choose Hm and Ybar as classification variables.
  Cluster formation:
  > param <- parameters15[c(\"Hm\",\"Ybar\")]
  > n <- dim(param)[1]; k <- dim(param)[2]
  > param.norm <- (param - matrix(mean(param),nrow=n,ncol=k,byrow=TRUE))/matrix(sd(param),nrow=n,ncol=k,byrow=TRUE)
  > clusters <- traceWminim(param.norm,4);
  > names(clusters) <- parameters15[\"cod\"][,]
  > clusters
  \n")

  param <- parameters15[c("Hm","Ybar")]
  n <- dim(param)[1]; k <- dim(param)[2]
  param.norm <- (param - matrix(mean(param),nrow=n,ncol=k,byrow=TRUE))/matrix(sd(param),nrow=n,ncol=k,byrow=TRUE)
  clusters <- traceWminim(param.norm,4);
  names(clusters) <- parameters15["cod"][,]
  print(clusters)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Plot of clusters:
  > par(mfrow=c(2,2))
  > plot(parameters15[c(\"Hm\",\"Ybar\")],col=clusters,pch=clusters,cex=0.6,
  +   main=\"Clusters in morphologic parameters space\",cex.main=1,font.main=1)
  > points(tapply(parameters15[\"Hm\"][,],clusters,mean),tapply(parameters15[\"Ybar\"][,],clusters,mean),
  +   col=c(1:4),pch=c(1:4))
  > legend(\"topleft\",paste(\"clust \",c(1:4)),col=c(1:4),pch=c(1:4),bty=\"n\")
  > grid()

  > plot(parameters15[c(\"Xbar\",\"Ybar\")],col=clusters,pch=clusters,
  +   main=\"Clusters in geographical space\",cex.main=1,font.main=1)
  > grid()

  > plot(D15lca,D15lcv,xlab=\"L-CA\",ylab=\"L-CV\",pch=clusters,col=clusters,
  +   main=\"Clusters in L-moments space\",cex.main=1,font.main=1)
  > grid()

  > D15lkur <- tapply(D15,cod15,Lkur)
  > plot(D15lca,D15lkur,xlab=\"L-CA\",ylab=\"L-kur\",pch=clusters,col=clusters,
  +   main=\"Clusters in L-moments space\",cex.main=1,font.main=1)
  > grid()
  > par(mfrow=c(1,1))
  \n")

  par(mfrow=c(2,2))
  plot(parameters15[c("Hm","Ybar")],col=clusters,pch=clusters,cex=0.6,
    main="Clusters in morphologic parameters space",cex.main=1,font.main=1)
  points(tapply(parameters15["Hm"][,],clusters,mean),tapply(parameters15["Ybar"][,],clusters,mean),
    col=c(1:4),pch=c(1:4))
  legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
  grid()

  plot(parameters15[c("Xbar","Ybar")],col=clusters,pch=clusters,
    main="Clusters in geographical space",cex.main=1,font.main=1)
  #legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
  grid()

  plot(D15lca,D15lcv,xlab="L-CA",ylab="L-CV",pch=clusters,col=clusters,
    main="Clusters in L-moments space",cex.main=1,font.main=1)
  #legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
  grid()
  
  D15lkur <- tapply(D15,cod15,Lkur)
  plot(D15lca,D15lkur,xlab="L-CA",ylab="L-kur",pch=clusters,col=clusters,
    main="Clusters in L-moments space",cex.main=1,font.main=1)
  #legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
  grid()
  par(mfrow=c(1,1))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Homogeneity test (Cluster 1):
  > fac <- factor(annualflows15[\"cod\"][,],levels=names(clusters[clusters==1]))
  > D15.1 <- annualflows15[!is.na(fac),\"dato\"]
  > cod15.1 <- factor(annualflows15[!is.na(fac),\"cod\"])
  > D15.1med <- tapply(D15.1,cod15.1,mean)
  > D15.1adim <- D15.1/unsplit(D15.1med,cod15.1)
  > regionalplotpos(D15.1adim,cod15.1,xlab=\"cluster1\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > HWs <- HW.tests(D15.1,cod15.1); HWs
  > mtext(paste(\"HW = \",round(HWs[1],2)),3,-1.5,adj=0.05)
  > nomi <- row.names(parameters15[clusters==1,])
  > legend(\"bottomright\",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty=\"n\")
  \n")

  fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==1]))
  D15.1 <- annualflows15[!is.na(fac),"dato"]
  cod15.1 <- factor(annualflows15[!is.na(fac),"cod"])
  D15.1med <- tapply(D15.1,cod15.1,mean)
  D15.1adim <- D15.1/unsplit(D15.1med,cod15.1)
  regionalplotpos(D15.1adim,cod15.1,xlab="cluster1",main="Empirical distributions",cex.main=1,font.main=1)
  HWs <- HW.tests(D15.1,cod15.1); HWs
  mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
  nomi <- row.names(parameters15[clusters==1,])
  legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Homogeneity test (Cluster 2):
  > fac <- factor(annualflows15[\"cod\"][,],levels=names(clusters[clusters==2]))
  > D15.2 <- annualflows15[!is.na(fac),\"dato\"]
  > cod15.2 <- factor(annualflows15[!is.na(fac),\"cod\"])
  > D15.2med <- tapply(D15.2,cod15.2,mean)
  > D15.2adim <- D15.2/unsplit(D15.2med,cod15.2)
  > regionalplotpos(D15.2adim,cod15.2,xlab=\"cluster2\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > HWs <- HW.tests(D15.2,cod15.2); HWs
  > mtext(paste(\"HW = \",round(HWs[1],2)),3,-1.5,adj=0.05)
  > nomi <- row.names(parameters15[clusters==2,])
  > legend(\"bottomright\",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty=\"n\")
  \n")

  fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==2]))
  D15.2 <- annualflows15[!is.na(fac),"dato"]
  cod15.2 <- factor(annualflows15[!is.na(fac),"cod"])
  D15.2med <- tapply(D15.2,cod15.2,mean)
  D15.2adim <- D15.2/unsplit(D15.2med,cod15.2)
  regionalplotpos(D15.2adim,cod15.2,xlab="cluster2",main="Empirical distributions",cex.main=1,font.main=1)
  HWs <- HW.tests(D15.2,cod15.2); HWs
  mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
  nomi <- row.names(parameters15[clusters==2,])
  legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Homogeneity test (Cluster 3):
  > fac <- factor(annualflows15[\"cod\"][,],levels=names(clusters[clusters==3]))
  > D15.3 <- annualflows15[!is.na(fac),\"dato\"]
  > cod15.3 <- factor(annualflows15[!is.na(fac),\"cod\"])
  > D15.3med <- tapply(D15.3,cod15.3,mean)
  > D15.3adim <- D15.3/unsplit(D15.3med,cod15.3)
  > regionalplotpos(D15.3adim,cod15.3,xlab=\"cluster3\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > HWs <- HW.tests(D15.3,cod15.3); HWs
  > mtext(paste(\"HW = \",round(HWs[1],2)),3,-1.5,adj=0.05)
  > nomi <- row.names(parameters15[clusters==3,])
  > legend(\"bottomright\",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty=\"n\")
  \n")

  fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==3]))
  D15.3 <- annualflows15[!is.na(fac),"dato"]
  cod15.3 <- factor(annualflows15[!is.na(fac),"cod"])
  D15.3med <- tapply(D15.3,cod15.3,mean)
  D15.3adim <- D15.3/unsplit(D15.3med,cod15.3)
  regionalplotpos(D15.3adim,cod15.3,xlab="cluster3",main="Empirical distributions",cex.main=1,font.main=1)
  HWs <- HW.tests(D15.3,cod15.3); HWs
  mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
  nomi <- row.names(parameters15[clusters==3,])
  legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Homogeneity test (Cluster 4):
  > fac <- factor(annualflows15[\"cod\"][,],levels=names(clusters[clusters==4]))
  > D15.4 <- annualflows15[!is.na(fac),\"dato\"]
  > cod15.4 <- factor(annualflows15[!is.na(fac),\"cod\"])
  > D15.4med <- tapply(D15.4,cod15.4,mean)
  > D15.4adim <- D15.4/unsplit(D15.4med,cod15.4)
  > regionalplotpos(D15.4adim,cod15.4,xlab=\"cluster4\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > HWs <- HW.tests(D15.4,cod15.4); HWs
  > mtext(paste(\"HW = \",round(HWs[1],2)),3,-1.5,adj=0.05)
  > nomi <- row.names(parameters15[clusters==4,])
  > legend(\"bottomright\",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty=\"n\")
  \n")

  fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==4]))
  D15.4 <- annualflows15[!is.na(fac),"dato"]
  cod15.4 <- factor(annualflows15[!is.na(fac),"cod"])
  D15.4med <- tapply(D15.4,cod15.4,mean)
  D15.4adim <- D15.4/unsplit(D15.4med,cod15.4)
  regionalplotpos(D15.4adim,cod15.4,xlab="cluster4",main="Empirical distributions",cex.main=1,font.main=1)
  HWs <- HW.tests(D15.4,cod15.4); HWs
  mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
  nomi <- row.names(parameters15[clusters==4,])
  legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Homogeneity test (justification of the test choise):
  > lmom4reg <- data.frame(rbind(Lmoments(D15.1adim),Lmoments(D15.2adim),Lmoments(D15.3adim),Lmoments(D15.4adim)))
  > Lspace.HWvsAD()
  > points(lmom4reg[,c(4,3)])
  > text(lmom4reg[,c(4,3)],row.names(lmom4reg),adj=c(-.5,-.5))
  \n")

  lmom4reg <- data.frame(rbind(Lmoments(D15.1adim),Lmoments(D15.2adim),Lmoments(D15.3adim),Lmoments(D15.4adim)))
  Lspace.HWvsAD()
  points(lmom4reg[,c(4,3)])
  text(lmom4reg[,c(4,3)],row.names(lmom4reg),adj=c(-.5,-.5))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Model selection (L-moments ratio diagram):
  > Lmoment.ratio.diagram()
  > points(Lmoments(D15.1adim)[4],Lmoments(D15.1adim)[5],pch=1,col=1)
  > points(Lmoments(D15.2adim)[4],Lmoments(D15.2adim)[5],pch=2,col=2)
  > points(Lmoments(D15.3adim)[4],Lmoments(D15.3adim)[5],pch=3,col=3)
  > points(Lmoments(D15.4adim)[4],Lmoments(D15.4adim)[5],pch=4,col=4)
  > legend(\"bottomleft\",paste(\"clust \",c(1:4)),col=c(1:4),pch=c(1:4),bty=\"n\")
  \n")

  Lmoment.ratio.diagram()
  points(Lmoments(D15.1adim)[4],Lmoments(D15.1adim)[5],pch=1,col=1)
  points(Lmoments(D15.2adim)[4],Lmoments(D15.2adim)[5],pch=2,col=2)
  points(Lmoments(D15.3adim)[4],Lmoments(D15.3adim)[5],pch=3,col=3)
  points(Lmoments(D15.4adim)[4],Lmoments(D15.4adim)[5],pch=4,col=4)
  legend("bottomleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Model selection (Cluster 1):
  > Lmom.1 <- Lmoments(D15.1adim); Lmom.1
  \n")

  Lmom.1 <- Lmoments(D15.1adim); print(Lmom.1)

  cat("\n
  > par1 <- par.gamma(Lmom.1[\"l1\"],Lmom.1[\"l2\"],Lmom.1[\"lca\"])
  > mom2par.gamma(par1$mu,par1$sigma,par1$gamm)
  \n")

  par1 <- par.gamma(Lmom.1["l1"],Lmom.1["l2"],Lmom.1["lca"]);
  print(mom2par.gamma(par1$mu,par1$sigma,par1$gamm))

  cat("\n
  > F1 <- F.gamma(D15.1adim,par1$mu,par1$sigma,par1$gamm)
  > regionalplotpos(D15.1adim,cod15.1,xlab=\"cluster1\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.1adim),sort(F1))
  > AD <- gofP3test(D15.1adim,Nsim=100)
  > mtext(paste(\"P(AD) = \",round(AD[2],3)),3,-1.5,adj=0.05)
  \n")

  F1 <- F.gamma(D15.1adim,par1$mu,par1$sigma,par1$gamm)
  regionalplotpos(D15.1adim,cod15.1,xlab="cluster1",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.1adim),sort(F1))
  AD <- gofP3test(D15.1adim,Nsim=100)
  mtext(paste("P(AD) = ",round(AD[2],3)),3,-1.5,adj=0.05)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Model selection (Cluster 2):
  > Lmom.2 <- Lmoments(D15.2adim); Lmom.2
  \n")

  Lmom.2 <- Lmoments(D15.2adim); print(Lmom.2)

  cat("\n
  > par2 <- par.gamma(Lmom.2[\"l1\"],Lmom.2[\"l2\"],Lmom.2[\"lca\"])
  > mom2par.gamma(par2$mu,par2$sigma,par2$gamm)
  \n")

  par2 <- par.gamma(Lmom.2["l1"],Lmom.2["l2"],Lmom.2["lca"])
  print(mom2par.gamma(par2$mu,par2$sigma,par2$gamm))

  cat("\n
  > F2 <- F.gamma(D15.2adim,par2$mu,par2$sigma,par2$gamm)
  > regionalplotpos(D15.2adim,cod15.2,xlab=\"cluster2\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.2adim),sort(F2))
  > AD <- gofP3test(D15.2adim,Nsim=100); AD
  > mtext(paste(\"P(AD) = \",round(AD[2],3)),3,-1.5,adj=0.05)
  \n")

  F2 <- F.gamma(D15.2adim,par2$mu,par2$sigma,par2$gamm)
  regionalplotpos(D15.2adim,cod15.2,xlab="cluster2",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.2adim),sort(F2))
  AD <- gofP3test(D15.2adim,Nsim=100); AD
  mtext(paste("P(AD) = ",round(AD[2],3)),3,-1.5,adj=0.05)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Model selection (Cluster 3):
  > Lmom.3 <- Lmoments(D15.3adim); Lmom.3 
  \n")

  Lmom.3 <- Lmoments(D15.3adim); print(Lmom.3)

  cat("\n
  > par3 <- par.gamma(Lmom.3[\"l1\"],Lmom.3[\"l2\"],Lmom.3[\"lca\"])
  > mom2par.gamma(par3$mu,par3$sigma,par3$gamm)
  \n")

  par3 <- par.gamma(Lmom.3["l1"],Lmom.3["l2"],Lmom.3["lca"])
  print(mom2par.gamma(par3$mu,par3$sigma,par3$gamm))

  cat("\n
  > F3 <- F.gamma(D15.3adim,par3$mu,par3$sigma,par3$gamm)
  > regionalplotpos(D15.3adim,cod15.3,xlab=\"cluster3\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.3adim),sort(F3))
  > AD <- gofP3test(D15.3adim,Nsim=100); AD
  > mtext(paste(\"P(AD) = \",round(AD[2],3)),3,-1.5,adj=0.05)
  \n")

  F3 <- F.gamma(D15.3adim,par3$mu,par3$sigma,par3$gamm)
  regionalplotpos(D15.3adim,cod15.3,xlab="cluster3",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.3adim),sort(F3))
  AD <- gofP3test(D15.3adim,Nsim=100); AD
  mtext(paste("P(AD) = ",round(AD[2],3)),3,-1.5,adj=0.05)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Model selection (Cluster 4):
  > Lmom.4 <- Lmoments(D15.4adim); Lmom.4 
  \n")

  Lmom.4 <- Lmoments(D15.4adim); print(Lmom.4)

  cat("\n
  > par4 <- par.gamma(Lmom.4[\"l1\"],Lmom.4[\"l2\"],Lmom.4[\"lca\"])
  > mom2par.gamma(par4$mu,par4$sigma,par4$gamm)
  \n")

  par4 <- par.gamma(Lmom.4["l1"],Lmom.4["l2"],Lmom.4["lca"])
  print(mom2par.gamma(par4$mu,par4$sigma,par4$gamm))

  cat("\n
  > F4 <- F.gamma(D15.4adim,par4$mu,par4$sigma,par4$gamm)
  > regionalplotpos(D15.4adim,cod15.4,xlab=\"cluster4\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.4adim),sort(F4))
  > AD <- gofP3test(D15.4adim,Nsim=100); AD
  > mtext(paste(\"P(AD) = \",round(AD[2],3)),3,-1.5,adj=0.05)
  \n")

  F4 <- F.gamma(D15.4adim,par4$mu,par4$sigma,par4$gamm)
  regionalplotpos(D15.4adim,cod15.4,xlab="cluster4",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.4adim),sort(F4))
  AD <- gofP3test(D15.4adim,Nsim=100); AD
  mtext(paste("P(AD) = ",round(AD[2],3)),3,-1.5,adj=0.05)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Comparison between regional growth-curves:
  > Fs <- seq(0.001,0.999,by=.001)
  > q1 <- invF.gamma(Fs,par1$mu,par1$sigma,par1$gamm)
  > q2 <- invF.gamma(Fs,par2$mu,par2$sigma,par2$gamm)
  > q3 <- invF.gamma(Fs,par3$mu,par3$sigma,par3$gamm)
  > q4 <- invF.gamma(Fs,par4$mu,par4$sigma,par4$gamm)

  > lognormplot(c(q1,q2,q3,q4),line=FALSE,type=\"n\")
  > normpoints(q1,type=\"l\",lty=1,col=1)
  > normpoints(q2,type=\"l\",lty=2,col=2)
  > normpoints(q3,type=\"l\",lty=3,col=3)
  > normpoints(q4,type=\"l\",lty=4,col=4)
  > legend(\"bottomright\",paste(\"cluster \",c(1:4)),col=c(1:4),lty=c(1:4),bty=\"n\")
  \n")

  Fs <- seq(0.001,0.999,by=.001)
  q1 <- invF.gamma(Fs,par1$mu,par1$sigma,par1$gamm)
  q2 <- invF.gamma(Fs,par2$mu,par2$sigma,par2$gamm)
  q3 <- invF.gamma(Fs,par3$mu,par3$sigma,par3$gamm)
  q4 <- invF.gamma(Fs,par4$mu,par4$sigma,par4$gamm)

  lognormplot(c(q1,q2,q3,q4),line=FALSE,type="n")
  normpoints(q1,type="l",lty=1,col=1)
  normpoints(q2,type="l",lty=2,col=2)
  normpoints(q3,type="l",lty=3,col=3)
  normpoints(q4,type="l",lty=4,col=4)
  legend("bottomright",paste("cluster ",c(1:4)),col=c(1:4),lty=c(1:4),bty="n")

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                                          THE END
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")
}
