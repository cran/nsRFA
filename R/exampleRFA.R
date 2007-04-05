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
  > par1
  \n")

  par1 <- par.gamma(Lmom.1["l1"],Lmom.1["l2"],Lmom.1["lca"]);
  print(par1)

  cat("\n
  > F1 <- F.gamma(D15.1adim,par1$xi,par1$beta,par1$alfa)
  > regionalplotpos(D15.1adim,cod15.1,xlab=\"cluster1\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.1adim),sort(F1))
  > AD <- gofP3test(D15.1adim,Nsim=100)
  > ADb <- A2_GOFlaio(D15.1adim, dist=\"GAM\")
  > AD2 <- gofGEVtest(D15.1adim,Nsim=100)
  > AD2b <- A2_GOFlaio(D15.1adim, dist=\"GEV\")
  > AD3 <- gofGENLOGIStest(D15.1adim,Nsim=100)
  > AD4 <- gofGENPARtest(D15.1adim,Nsim=100)
  > AD5 <- gofLOGNORMtest(D15.1adim,Nsim=100)
  > mtext(paste(\"P(AD) P3 = \",round(AD[2],3),\", \",round(ADb[2],3)),3,-1.5,adj=0.05)
  > mtext(paste(\"P(AD) GEV = \",round(AD2[2],3),\", \",round(AD2b[2],3)),3,-2.5,adj=0.05)
  > mtext(paste(\"P(AD) GL = \",round(AD3[2],3)),3,-3.5,adj=0.05)
  > mtext(paste(\"P(AD) GP = \",round(AD4[2],3)),3,-4.5,adj=0.05)
  > mtext(paste(\"P(AD) LN = \",round(AD5[2],3)),3,-5.5,adj=0.05)
  \n")

  F1 <- F.gamma(D15.1adim,par1$xi,par1$beta,par1$alfa)
  regionalplotpos(D15.1adim,cod15.1,xlab="cluster1",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.1adim),sort(F1))
  AD <- gofP3test(D15.1adim,Nsim=100)
  ADb <- A2_GOFlaio(D15.1adim, dist="GAM")
  AD2 <- gofGEVtest(D15.1adim,Nsim=100)
  AD2b <- A2_GOFlaio(D15.1adim, dist="GEV")
  AD3 <- gofGENLOGIStest(D15.1adim,Nsim=100)
  AD4 <- gofGENPARtest(D15.1adim,Nsim=100)
  AD5 <- gofLOGNORMtest(D15.1adim,Nsim=100)
  mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
  mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
  mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
  mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
  mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)

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
  > par2
  \n")

  par2 <- par.gamma(Lmom.2["l1"],Lmom.2["l2"],Lmom.2["lca"])
  print(par2)

  cat("\n
  > F2 <- F.gamma(D15.2adim,par2$xi,par2$beta,par2$alfa)
  > regionalplotpos(D15.2adim,cod15.2,xlab=\"cluster2\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.2adim),sort(F2))
  > AD <- gofP3test(D15.2adim,Nsim=100)
  > ADb <- A2_GOFlaio(D15.2adim, dist=\"GAM\")
  > AD2 <- gofGEVtest(D15.2adim,Nsim=100)
  > AD2b <- A2_GOFlaio(D15.2adim, dist=\"GEV\")
  > AD3 <- gofGENLOGIStest(D15.2adim,Nsim=100)
  > AD4 <- gofGENPARtest(D15.2adim,Nsim=100)
  > AD5 <- gofLOGNORMtest(D15.2adim,Nsim=100)
  > mtext(paste(\"P(AD) P3 = \",round(AD[2],3),\", \",round(ADb[2],3)),3,-1.5,adj=0.05)
  > mtext(paste(\"P(AD) GEV = \",round(AD2[2],3),\", \",round(AD2b[2],3)),3,-2.5,adj=0.05)
  > mtext(paste(\"P(AD) GL = \",round(AD3[2],3)),3,-3.5,adj=0.05)
  > mtext(paste(\"P(AD) GP = \",round(AD4[2],3)),3,-4.5,adj=0.05)
  > mtext(paste(\"P(AD) LN = \",round(AD5[2],3)),3,-5.5,adj=0.05)
  \n")

  F2 <- F.gamma(D15.2adim,par2$xi,par2$beta,par2$alfa)
  regionalplotpos(D15.2adim,cod15.2,xlab="cluster2",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.2adim),sort(F2))
  AD <- gofP3test(D15.2adim,Nsim=100)
  ADb <- A2_GOFlaio(D15.2adim, dist="GAM")
  AD2 <- gofGEVtest(D15.2adim,Nsim=100)
  AD2b <- A2_GOFlaio(D15.2adim, dist="GEV")
  AD3 <- gofGENLOGIStest(D15.2adim,Nsim=100)
  AD4 <- gofGENPARtest(D15.2adim,Nsim=100)
  AD5 <- gofLOGNORMtest(D15.2adim,Nsim=100)
  mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
  mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
  mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
  mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
  mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)

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
  > par3
  \n")

  par3 <- par.gamma(Lmom.3["l1"],Lmom.3["l2"],Lmom.3["lca"])
  print(par3)

  cat("\n
  > F3 <- F.gamma(D15.3adim,par3$xi,par3$beta,par3$alfa)
  > regionalplotpos(D15.3adim,cod15.3,xlab=\"cluster3\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.3adim),sort(F3))
  > AD <- gofP3test(D15.3adim,Nsim=100)
  > ADb <- A2_GOFlaio(D15.3adim, dist=\"GAM\")
  > AD2 <- gofGEVtest(D15.3adim,Nsim=100)
  > AD2b <- A2_GOFlaio(D15.3adim, dist=\"GEV\")
  > AD3 <- gofGENLOGIStest(D15.3adim,Nsim=100)
  > AD4 <- gofGENPARtest(D15.3adim,Nsim=100)
  > AD5 <- gofLOGNORMtest(D15.3adim,Nsim=100)
  > mtext(paste(\"P(AD) P3 = \",round(AD[2],3),\", \",round(ADb[2],3)),3,-1.5,adj=0.05)
  > mtext(paste(\"P(AD) GEV = \",round(AD2[2],3),\", \",round(AD2b[2],3)),3,-2.5,adj=0.05)
  > mtext(paste(\"P(AD) GL = \",round(AD3[2],3)),3,-3.5,adj=0.05)
  > mtext(paste(\"P(AD) GP = \",round(AD4[2],3)),3,-4.5,adj=0.05)
  > mtext(paste(\"P(AD) LN = \",round(AD5[2],3)),3,-5.5,adj=0.05)
  \n")

  F3 <- F.gamma(D15.3adim,par3$xi,par3$beta,par3$alfa)
  regionalplotpos(D15.3adim,cod15.3,xlab="cluster3",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.3adim),sort(F3))
  AD <- gofP3test(D15.3adim,Nsim=100)
  ADb <- A2_GOFlaio(D15.3adim, dist="GAM")
  AD2 <- gofGEVtest(D15.3adim,Nsim=100)
  AD2b <- A2_GOFlaio(D15.3adim, dist="GEV")
  AD3 <- gofGENLOGIStest(D15.3adim,Nsim=100)
  AD4 <- gofGENPARtest(D15.3adim,Nsim=100)
  AD5 <- gofLOGNORMtest(D15.3adim,Nsim=100)
  mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
  mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
  mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
  mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
  mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)

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
  > par4
  \n")

  par4 <- par.gamma(Lmom.4["l1"],Lmom.4["l2"],Lmom.4["lca"])
  print(par4)

  cat("\n
  > F4 <- F.gamma(D15.4adim,par4$xi,par4$beta,par4$alfa)
  > regionalplotpos(D15.4adim,cod15.4,xlab=\"cluster4\",main=\"Empirical distributions\",cex.main=1,font.main=1)
  > lines(sort(D15.4adim),sort(F4))
  > AD <- gofP3test(D15.4adim,Nsim=100)
  > ADb <- A2_GOFlaio(D15.4adim, dist=\"GAM\")
  > AD2 <- gofGEVtest(D15.4adim,Nsim=100)
  > AD2b <- A2_GOFlaio(D15.4adim, dist=\"GEV\")
  > AD3 <- gofGENLOGIStest(D15.4adim,Nsim=100)
  > AD4 <- gofGENPARtest(D15.4adim,Nsim=100)
  > AD5 <- gofLOGNORMtest(D15.4adim,Nsim=100)
  > mtext(paste(\"P(AD) P3 = \",round(AD[2],3),\", \",round(ADb[2],3)),3,-1.5,adj=0.05)
  > mtext(paste(\"P(AD) GEV = \",round(AD2[2],3),\", \",round(AD2b[2],3)),3,-2.5,adj=0.05)
  > mtext(paste(\"P(AD) GL = \",round(AD3[2],3)),3,-3.5,adj=0.05)
  > mtext(paste(\"P(AD) GP = \",round(AD4[2],3)),3,-4.5,adj=0.05)
  > mtext(paste(\"P(AD) LN = \",round(AD5[2],3)),3,-5.5,adj=0.05)
  \n")

  F4 <- F.gamma(D15.4adim,par4$xi,par4$beta,par4$alfa)
  regionalplotpos(D15.4adim,cod15.4,xlab="cluster4",main="Empirical distributions",cex.main=1,font.main=1)
  lines(sort(D15.4adim),sort(F4))
  AD <- gofP3test(D15.4adim,Nsim=100)
  ADb <- A2_GOFlaio(D15.4adim, dist="GAM")
  AD2 <- gofGEVtest(D15.4adim,Nsim=100)
  AD2b <- A2_GOFlaio(D15.4adim, dist="GEV")
  AD3 <- gofGENLOGIStest(D15.4adim,Nsim=100)
  AD4 <- gofGENPARtest(D15.4adim,Nsim=100)
  AD5 <- gofLOGNORMtest(D15.4adim,Nsim=100)
  mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
  mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
  mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
  mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
  mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Comparison between regional growth-curves:
  > Fs <- seq(0.001,0.999,by=.001)
  > q1 <- invF.gamma(Fs,par1$xi,par1$beta,par1$alfa)
  > q2 <- invF.gamma(Fs,par2$xi,par2$beta,par2$alfa)
  > q3 <- invF.gamma(Fs,par3$xi,par3$beta,par3$alfa)
  > q4 <- invF.gamma(Fs,par4$xi,par4$beta,par4$alfa)

  > lognormplot(c(q1,q2,q3,q4),line=FALSE,type=\"n\")
  > normpoints(q1,type=\"l\",lty=1,col=1)
  > normpoints(q2,type=\"l\",lty=2,col=2)
  > normpoints(q3,type=\"l\",lty=3,col=3)
  > normpoints(q4,type=\"l\",lty=4,col=4)
  > legend(\"bottomright\",paste(\"cluster \",c(1:4)),col=c(1:4),lty=c(1:4),bty=\"n\")
  \n")

  Fs <- seq(0.001,0.999,by=.001)
  q1 <- invF.gamma(Fs,par1$xi,par1$beta,par1$alfa)
  q2 <- invF.gamma(Fs,par2$xi,par2$beta,par2$alfa)
  q3 <- invF.gamma(Fs,par3$xi,par3$beta,par3$alfa)
  q4 <- invF.gamma(Fs,par4$xi,par4$beta,par4$alfa)

  lognormplot(c(q1,q2,q3,q4),line=FALSE,type="n")
  normpoints(q1,type="l",lty=1,col=1)
  normpoints(q2,type="l",lty=2,col=2)
  normpoints(q3,type="l",lty=3,col=3)
  normpoints(q4,type="l",lty=4,col=4)
  legend("bottomright",paste("cluster ",c(1:4)),col=c(1:4),lty=c(1:4),bty="n")


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Final tests of goodness-of-fit. 
  
  The regionalized distribution adaptation, with mean 
    D.est = exp(7.86 + 0.000291*Hm + 0.0722*NORD - 1.70*IB)
  and parameters of Pearson type III distribution
    xi1 = 0.5458;  beta1 = 0.09534; alfa1 = 4.764;
    xi2 = 0.09843; beta2 = 0.08508; alfa2 = 10.60;
    xi3 = 0.2801;  beta3 = 0.1237;  alfa3 = 5.817;
    xi4 = 0.3496;  beta4 = 0.2093;  alfa4 = 3.107
  is tested throught the Anderson-Darling and the Cramer-von Mises
  goodness-of-fit tests.
  In this case (fully specified distribution) the percentage points
  are:
         .25    .10    .05    .01   .001
    A2 1.248  1.933  2.492  3.880  6.000
    W2 0.209  0.347  0.461  0.743  1.167  \n

  > A2fullspecif <- function (x,xi,beta,alfa) {
  >  x <- sort(x)
  >  n <- length(x)
  >  F <- pgamma((x - xi)/beta, alfa)
  >  F[F<=0] <- 0.000001
  >  F[F>=1] <- 0.999999
  >  -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
  > }
  > W2fullspecif <- function (x,xi,beta,alfa) {
  >  x <- sort(x)
  >  n <- length(x)
  >  F <- pgamma((x - xi)/beta, alfa)
  >  sum((F - seq(1,2*n-1,by=2)/(2*n))^2) + 1/(12*n)
  > }
  > par(mfrow=c(3,2))
  > W2s <- rep(NA,47)
  > A2s <- rep(NA,47)
  > em <- rep(NA,47)
  > gruppi <- rep(NA,47)
  > ferma=6
  > for(i in 1:47) {
  >  codice <- unique(cod)[i]
  >  Hm <- parameters[parameters[\"cod\"]==codice,\"Hm\"]
  >  NORD <- parameters[parameters[\"cod\"]==codice,\"NORD\"]
  >  IB <- parameters[parameters[\"cod\"]==codice,\"IB\"]
  >  Ybar <- parameters[parameters[\"cod\"]==codice,\"Ybar\"]
  >  baricentri.clusterD <- cbind(c(2327.455,1404.200,1892.214,812.125),c(45.77200,46.01480,44.58071,44.47513))
  >  baricentri.clusterD[,1] <- (baricentri.clusterD[,1]-1694)/652.1
  >  baricentri.clusterD[,2] <- (baricentri.clusterD[,2]-45.12)/0.7229
  >  HmYbar <- c((Hm-1694)/652.1,(Ybar-45.12)/0.7229)
  >  distanze <- as.matrix(dist(rbind(HmYbar,baricentri.clusterD)))[-1,1]
  >  gruppo <- which.min(distanze)
  >  em[i] <- exp(7.86 + 0.000291*Hm + 0.0722*NORD - 1.70*IB)
  >  if(gruppo==1){
  >   xi=0.5458; beta=0.09534; alfa=4.764;
  >  }
  >  else if(gruppo==2){
  >   xi=0.09843; beta=0.08508; alfa=10.60;
  >  }
  >  else if(gruppo==3){
  >   xi=0.2801; beta=0.1237; alfa=5.817;
  >  }
  >  else if(gruppo==4){
  >   xi=0.3496; beta=0.2093; alfa=3.107;
  >  }
  >  gruppi[i] <- gruppo
  >  campione <- annualflows[annualflows[\"cod\"]==codice,3]
  >  plotpos(campione,main=row.names(parameters[parameters[\"cod\"]==codice,]),cex.main=1,font.main=1)
  >  lines(sort(campione),pgamma((sort(campione)-em[i]*xi)/(em[i]*beta),alfa))
  >  A2s[i] <- A2fullspecif(campione,em[i]*xi,em[i]*beta,alfa)
  >  W2s[i] <- W2fullspecif(campione,em[i]*xi,em[i]*beta,alfa)
  >  mtext(paste(\"A2 = \",round(A2s[i],3)),3,-1.5,adj=0.05)
  >  mtext(paste(\"W2 = \",round(W2s[i],3)),1,-1.5,adj=0.95)
  >  if(i==ferma) {
  >   ferma <- ferma+6
  >   readline(\"Press Return to continue to the next graphs\")
  >  }
  > }
  > tabella <- cbind(parameters[1:2],round(em),gruppi,A2s,W2s)
  > tabella
  > par(mfrow=c(1,1))
  \n")

  A2fullspecif <- function (x,xi,beta,alfa) {
   x <- sort(x)
   n <- length(x)
   F <- pgamma((x - xi)/beta, alfa)
   F[F<=0] <- 0.000001
   F[F>=1] <- 0.999999
   -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
  } 
  W2fullspecif <- function (x,xi,beta,alfa) {
   x <- sort(x)
   n <- length(x)
   F <- pgamma((x - xi)/beta, alfa)
   sum((F - seq(1,2*n-1,by=2)/(2*n))^2) + 1/(12*n)
  }
  par(mfrow=c(3,2))
  W2s <- rep(NA,47)
  A2s <- rep(NA,47)
  em <- rep(NA,47)
  gruppi <- rep(NA,47)
  ferma=6
  for(i in 1:47) {
   codice <- unique(cod)[i]
   Hm <- parameters[parameters["cod"]==codice,"Hm"]
   NORD <- parameters[parameters["cod"]==codice,"NORD"]
   IB <- parameters[parameters["cod"]==codice,"IB"]
   Ybar <- parameters[parameters["cod"]==codice,"Ybar"]
   baricentri.clusterD <- cbind(c(2327.455,1404.200,1892.214,812.125),c(45.77200,46.01480,44.58071,44.47513))
   baricentri.clusterD[,1] <- (baricentri.clusterD[,1]-1694)/652.1
   baricentri.clusterD[,2] <- (baricentri.clusterD[,2]-45.12)/0.7229
   HmYbar <- c((Hm-1694)/652.1,(Ybar-45.12)/0.7229)
   distanze <- as.matrix(dist(rbind(HmYbar,baricentri.clusterD)))[-1,1]
   gruppo <- which.min(distanze)
   em[i] <- exp(7.86 + 0.000291*Hm + 0.0722*NORD - 1.70*IB)
   if(gruppo==1){
    xi=0.5458; beta=0.09534; alfa=4.764;
   }
   else if(gruppo==2){
    xi=0.09843; beta=0.08508; alfa=10.60;
   }
   else if(gruppo==3){
    xi=0.2801; beta=0.1237; alfa=5.817;
   }
   else if(gruppo==4){
    xi=0.3496; beta=0.2093; alfa=3.107;
   }
   gruppi[i] <- gruppo
   campione <- annualflows[annualflows["cod"]==codice,3]
   plotpos(campione,main=row.names(parameters[parameters["cod"]==codice,]),cex.main=1,font.main=1)
   lines(sort(campione),pgamma((sort(campione)-em[i]*xi)/(em[i]*beta),alfa))
   A2s[i] <- A2fullspecif(campione,em[i]*xi,em[i]*beta,alfa)
   W2s[i] <- W2fullspecif(campione,em[i]*xi,em[i]*beta,alfa)
   mtext(paste("A2 = ",round(A2s[i],3)),3,-1.5,adj=0.05) 
   mtext(paste("W2 = ",round(W2s[i],3)),1,-1.5,adj=0.95)
   if(i==ferma) {
    ferma <- ferma+6
    readline("Press Return to continue to the next graphs")
   }
  }
  #tabella <- cbind(parameters[1],round(tapply(D,cod,mean)),round(em),gruppi,A2s,W2s)
  tabella <- cbind(parameters[1:2],round(em),gruppi,A2s,W2s)
  names(tabella) <- c("cod","tildeDm","hatDm","group","A2","W2")
  print(tabella) 
  par(mfrow=c(1,1))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Dependency on the estimation of the mean.
  > errors <- abs((tabella[,2] - tabella[,3])/tabella[,2])
  > codici <- tabella[,1]
  > plot(A2s,errors,main=\"\",xlab=expression(A^2),ylab=expression(paste(\"|(\",hat(D)[m] - tilde(D)[m],\")/\",tilde(D)[m],\"|\")))
  > grid()
  > abline(v=2.492,lty=2,col=2)
  > abline(v=3.880,lty=2,col=2)
  > text(A2s[A2s>3.880],errors[A2s>3.880],labels=codici[A2s>3.880]-2,pos=1)
  > axis(1,at=c(2.492,3.880),labels=c(\"5\%\",\"1\%\"),col=2,col.axis=2,cex.axis=.8)
  \n")

  errors <- abs((tabella[,2] - tabella[,3])/tabella[,2])
  codici <- tabella[,1]
  plot(A2s,errors,main="",xlab=expression(A^2),ylab=expression(paste("|(",hat(D)[m] - tilde(D)[m],")/",tilde(D)[m],"|")))
  grid()
  abline(v=2.492,lty=2,col=2)
  abline(v=3.880,lty=2,col=2)
  text(A2s[A2s>3.880],errors[A2s>3.880],labels=codici[A2s>3.880]-2,pos=1)
  axis(1,at=c(2.492,3.880),labels=c("5\%","1\%"),col=2,col.axis=2,cex.axis=.8)

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                                          THE END
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")
}






































exampleRFA03 <- function () {

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                           ANALYSIS OF FEH DATA: POOLING GROUPS
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")

  readline("\n
  Press Return to continue
  \n")

  cat("\n
  Data loading: \n
  > data(FEH1000)
  \n")

  data(FEH1000)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Sample moments and L-moments:
  > momenti <- t(sapply(split(am[,4],am[,1]),moments))
  > Lmomenti <- t(sapply(split(am[,4],am[,1]),Lmoments))
  > plot(as.data.frame(momenti))
  > plot(as.data.frame(Lmomenti))
  \n")

  momenti <- t(sapply(split(am[,4],am[,1]),moments))
  Lmomenti <- t(sapply(split(am[,4],am[,1]),Lmoments))
  plot(as.data.frame(momenti),pch=".",cex=2)
  readline("Press Return to continue to the next graphs")
  plot(as.data.frame(Lmomenti),pch=".",cex=2)



  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Criteria used in the FEH to choose stations for pooling groups: 
   n>7
   area, saar and bfihost are known
   urbext<0.025
   area>0.5 \n
  > n <- tapply(am[,3],am[,1],length)
  > urbext <- cd[,\"urbext1990\"]
  > area <- cd[,\"dtm_area\"]
  > cd696 <- cd[(!is.nan(cd[,\"dtm_area\"]))&(!is.nan(cd[,\"saar\"]))&(!is.nan(cd[,\"bfihost\"]))&(n>7)&(urbext<0.025)&(area>0.5),]

  > fac <- factor(am[,\"number\"],levels=cd696[,\"number\"])
  > am696 <- am[!is.na(fac),]
  > nlevels(as.factor(am696[,\"number\"]))
  \n")

  n <- tapply(am[,4],am[,1],length)
  urbext <- cd[,"urbext1990"]
  area <- cd[,"dtm_area"]
  cd696 <- cd[(!is.nan(cd[,"dtm_area"]))&(!is.nan(cd[,"saar"]))&(!is.nan(cd[,"bfihost"]))&(n>7)&(urbext<0.025)&(area>0.5),]
  
  fac <- factor(am[,"number"],levels=cd696[,"number"])
  am696 <- am[!is.na(fac),]
  nlevels(as.factor(am696[,"number"]))


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Figure 16.2 pag.157, FEH Vol.3 \n
  > plot(cd696[c(\"dtm_area\",\"saar\")], pch=\".\", cex=2, log=\"x\")
  > plot(cd696[c(\"dtm_area\",\"bfihost\")], pch=\".\", cex=2, log=\"x\")
  > plot(cd696[c(\"saar\",\"bfihost\")], pch=\".\", cex=2)
  \n")
  
  plot(cd696[c("dtm_area","saar")], pch=".", cex=2, log="x")
  readline("Press Return to continue to the next graphs")
  plot(cd696[c("dtm_area","bfihost")], pch=".", cex=2, log="x")
  readline("Press Return to continue to the next graphs")
  plot(cd696[c("saar","bfihost")], pch=".", cex=2)
  

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Discordancy measure: \n
  > Lmomenti696 <- t(sapply(split(am696[,4],am696[,1]),Lmoments))
  > Di <- discordancy(am696[,\"am\"], am696[,\"number\"])
  > par(mfrow=c(1,2))
  >  plot(Lmomenti696[,c(\"lca\",\"lcv\")],xlab=\"L-CA\",ylab=\"L-CV\",cex=.7); grid()
  >  points(Lmomenti696[(Di>3),c(\"lca\",\"lcv\")],pch=19,cex=.7)
  >  plot(Lmomenti696[,c(\"lca\",\"lkur\")],xlab=\"L-CA\",ylab=\"L-kur\",cex=.7); grid()
  >  points(Lmomenti696[(Di>3),c(\"lca\",\"lkur\")],pch=19,cex=.7)
  > par(mfrow=c(1,1))
  \n")

  Lmomenti696 <- t(sapply(split(am696[,4],am696[,1]),Lmoments))
  Di <- discordancy(am696[,"am"], am696[,"number"])
  par(mfrow=c(1,2))
   plot(Lmomenti696[,c("lca","lcv")],xlab="L-CA",ylab="L-CV",pch=".",cex=2); grid()
   points(Lmomenti696[(Di>3),c("lca","lcv")],pch=19,cex=.7)
   plot(Lmomenti696[,c("lca","lkur")],xlab="L-CA",ylab="L-kur",pch=".",cex=2); grid()
   points(Lmomenti696[(Di>3),c("lca","lkur")],pch=19,cex=.7)
  par(mfrow=c(1,1))


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Region of influence (Table 16.2, pag.164, FEH Vol.3):
  using lnAREA, lnSAAR and BFIHOST to measure distances among sites:\n
  > sd(log(cd696[,\"dtm_area\"])) 
  > sd(log(cd696[,\"saar\"]))     
  > sd(cd696[,\"bfihost\"])       

  > AREAterm <- log(cd696[,\"dtm_area\"])/(sd(log(cd696[,\"dtm_area\"]))*sqrt(2))
  > SAARterm <- log(cd696[,\"saar\"])/sd(log(cd696[,\"saar\"]))
  > BFIHOSTterm <- cd696[,\"bfihost\"]/sd(cd696[,\"bfihost\"])

  > distFEH <- dist(cbind(AREAterm,SAARterm,BFIHOSTterm))

  > roi.cd <- data.frame(cbind(AREAterm,SAARterm,BFIHOSTterm))
  > row.names(roi.cd) <- cd696[,\"number\"]

  > roi01.50year <- new.env()
  > for(i in 1:696) {
  >  print(paste(i,\"/ 696\"))
  >  assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  +         row.names(roi.cd),am696[,\"am\"],am696[,\"number\"],test=\"HW\",station.year=250,Nsim=100), env=roi01.50year)
  > }
  > roi01.50year <- as.list(roi01.50year)

  > estrai.region <- function (x) {x$region}
  > estrai.test <- function (x) {x$test}
  > regioni.50year <- sapply(roi01.50year, estrai.region)
  > test.50year <- sapply(roi01.50year, estrai.test)
  > mL.50year <- mean(sapply(regioni.50year,length)) 
  > mH2.50year <- mean(test.50year[\"H2\",]) 
  > gH2gr2.50year <- sum(test.50year[\"H2\",]>2)/696 
  > gH2gr4.50year <- sum(test.50year[\"H2\",]>4)/696 

  > roi01.100year <- new.env()
  > for(i in 1:696) {
  >  print(paste(i,\"/ 696\"))
  >  assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  +         row.names(roi.cd),am696[,\"am\"],am696[,\"number\"],test=\"HW\",station.year=500,Nsim=100), env=roi01.100year)
  > }
  > roi01.100year <- as.list(roi01.100year)

  > regioni.100year <- sapply(roi01.100year, estrai.region)
  > test.100year <- sapply(roi01.100year, estrai.test)
  > mL.100year <- mean(sapply(regioni.100year,length)) 
  > mH2.100year <- mean(test.100year[\"H2\",]) 
  > gH2gr2.100year <- sum(test.100year[\"H2\",]>2)/696 
  > gH2gr4.100year <- sum(test.100year[\"H2\",]>4)/696 

  > table16.2 <- data.frame(signif(rbind(c(mL.50year,mH2.50year,gH2gr2.50year*100,gH2gr4.50year*100),
  +              c(mL.100year,mH2.100year,gH2gr2.100year*100,gH2gr4.100year*100)),3), row.names=c(\"50-year\",\"100-year\"))
  > names(table16.2) <- c(\"Avg. n sites\",\"m(H2)\",\"% H2>2\",\"% H2>4\")
  > print(table16.2)
  \n")

  sd(log(cd696[,"dtm_area"])) # 1.345515 (vs 1.34)
  sd(log(cd696[,"saar"]))     # 0.38534 (vs 0.38)
  sd(cd696[,"bfihost"])       # 0.1485239 (vs 0.15)

  AREAterm <- log(cd696[,"dtm_area"])/(sd(log(cd696[,"dtm_area"]))*sqrt(2))
  SAARterm <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
  BFIHOSTterm <- cd696[,"bfihost"]/sd(cd696[,"bfihost"])

  distFEH <- dist(cbind(AREAterm,SAARterm,BFIHOSTterm))

  roi.cd <- data.frame(cbind(AREAterm,SAARterm,BFIHOSTterm))
  row.names(roi.cd) <- cd696[,"number"]

  # roi01.50year <- new.env()
  # for(i in 1:696) {
  #  print(paste(i,"/ 696"))
  #  assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  #         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW",station.year=250,Nsim=100), env=roi01.50year)
  # }
  # roi01.50year <- as.list(roi01.50year)

  # estrai.region <- function (x) {x$region}
  # estrai.test <- function (x) {x$test}
  # regioni.50year <- sapply(roi01.50year, estrai.region)
  # test.50year <- sapply(roi01.50year, estrai.test)
  # mL.50year <- mean(sapply(regioni.50year,length)) #  11.2
  # mH2.50year <- mean(test.50year["H2",]) #   1.53
  # gH2gr2.50year <- sum(test.50year["H2",]>2)/696 #   0.34
  # gH2gr4.50year <- sum(test.50year["H2",]>4)/696 #   0.07

  # roi01.100year <- new.env()
  # for(i in 1:696) {
  #  print(paste(i,"/ 696"))
  #  assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  #         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW",station.year=500,Nsim=100), env=roi01.100year)
  # }
  # roi01.100year <- as.list(roi01.100year)

  # regioni.100year <- sapply(roi01.100year, estrai.region)
  # test.100year <- sapply(roi01.100year, estrai.test)
  # mL.100year <- mean(sapply(regioni.100year,length)) #  21.8
  # mH2.100year <- mean(test.100year["H2",]) #   2.19
  # gH2gr2.100year <- sum(test.100year["H2",]>2)/696 #   0.52
  # gH2gr4.100year <- sum(test.100year["H2",]>4)/696 #   0.15

  # table16.2 <- data.frame(signif(rbind(c(mL.50year,mH2.50year,gH2gr2.50year*100,gH2gr4.50year*100),
  #              c(mL.100year,mH2.100year,gH2gr2.100year*100,gH2gr4.100year*100)),3), row.names=c("50-year","100-year"))
  # names(table16.2) <- c("Avg. n sites","m(H2)","% H2>2","% H2>4")
  # print(table16.2)
  table16.2 <- data.frame(signif(rbind(c(11.2,1.53,0.34*100,0.07*100),
               c(21.8,2.19,0.52*100,0.15*100)),3), row.names=c("50-year","100-year"))
  names(table16.2) <- c("Avg. n sites","m(H2)","% H2>2","% H2>4")
  print(table16.2)



  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Example 16.3 pag.164, FEH Vol.3:\n
  > roi01.50year$\"54088\"
  > print(length(roi01.50year$\"54088\"$region))
  [1] 11
  > selez.reg <- factor(am[,\"number\"], levels=roi01.50year$\"54088\"$region)
  > print(length(am[!is.na(selez.reg),\"am\"]))
  [1] 274
  > print(HW.tests(am[!is.na(selez.reg),\"am\"], am[!is.na(selez.reg),\"number\"], 1000))
        H1       H2
  6.564278 4.187761

  > prova28018 <- roi.st.year(roi.cd[\"28018\",],roi.cd,row.names(roi.cd),am696[,\"am\"],
  >                           am696[,\"number\"],test=\"HW\",station.year=250,Nsim=500)
  > selez.reg <- factor(am[,\"number\"], levels=prova28018$region)
  > print(length(am[!is.na(selez.reg),\"am\"]))
  [1] 258
  > print(HW.tests(am[!is.na(selez.reg),\"am\"], am[!is.na(selez.reg),\"number\"], 1000))
          H1         H2
  -0.6000429 -1.0508830

  > Lmomenti696 <- as.data.frame(Lmomenti696)
  > par(mfrow=c(1,2))
  >  plot(Lmomenti696[c(\"lca\",\"lcv\")],xlab=\"L-CA\",ylab=\"L-CV\",cex=.7); grid()
  >  points(Lmomenti696[c(\"54088\"),c(\"lca\",\"lcv\")],pch=19,col=\"red\",cex=1)
  >  points(Lmomenti696[roi01.50year$\"54088\"$region[-1],c(\"lca\",\"lcv\")],pch=19,cex=1)
  >  plot(Lmomenti696[,c(\"lca\",\"lkur\")],xlab=\"L-CA\",ylab=\"L-kur\",cex=.7); grid()
  >  points(Lmomenti696[c(\"28018\"),c(\"lca\",\"lcv\")],pch=19,col=\"red\",cex=1)
  >  points(Lmomenti696[roi01.50year$\"28018\"$region[-1],c(\"lca\",\"lcv\")],pch=19,cex=1)
  > par(mfrow=c(1,1))
  \n")

  # roi01.50year$"54088"
  # print(length(roi01.50year$"54088"$region))   # 11 OK
  # selez.reg <- factor(am[,"number"], levels=roi01.50year$"54088"$region)
  # print(length(am[!is.na(selez.reg),"am"]))   # 274 OK
  # print(HW.tests(am[!is.na(selez.reg),"am"], am[!is.na(selez.reg),"number"], 1000))   # 4.03 vs 4.08  
  prova54088 <- roi.st.year(roi.cd["54088",],roi.cd,row.names(roi.cd),am696[,"am"],
                           am696[,"number"],test="HW",station.year=250,Nsim=500)
  # selez.reg <- factor(am[,"number"], levels=prova54088$region)
  # print(length(am[!is.na(selez.reg),"am"]))   # 258 OK
  # print(HW.tests(am[!is.na(selez.reg),"am"], am[!is.na(selez.reg),"number"], 1000))

  prova28018 <- roi.st.year(roi.cd["28018",],roi.cd,row.names(roi.cd),am696[,"am"],
                            am696[,"number"],test="HW",station.year=250,Nsim=500)
  # selez.reg <- factor(am[,"number"], levels=prova28018$region)
  # print(length(am[!is.na(selez.reg),"am"]))   # 258 OK
  # print(HW.tests(am[!is.na(selez.reg),"am"], am[!is.na(selez.reg),"number"], 1000))   # -1.06 vs -0.98

  Lmomenti696 <- as.data.frame(Lmomenti696)
  par(mfrow=c(1,2))
   plot(Lmomenti696[c("lca","lcv")],xlab="L-CA",ylab="L-CV",pch=".",cex=2,main="54088"); grid()
   points(Lmomenti696[c("54088"),c("lca","lcv")],pch=19,col="red",cex=1)
   points(Lmomenti696[prova54088$region[-1],c("lca","lcv")],pch=19,cex=1)
   plot(Lmomenti696[,c("lca","lkur")],xlab="L-CA",ylab="L-kur",pch=".",cex=2,main="28018"); grid()
   points(Lmomenti696[c("28018"),c("lca","lcv")],pch=19,col="red",cex=1)
   points(Lmomenti696[prova28018$region[-1],c("lca","lcv")],pch=19,cex=1)
  par(mfrow=c(1,1))



  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Figure 16.9 pag.174 (1st part), FEH Vol.3:\n
  > figure16.9a <- function (x,r,cd) {
  +  # x = station of interest (e.g. \"28018\")
  +  # r = output of roi.st.year()
  +
  +  if(!r$region[1]==x) r$region <- c(x,r$region)
  +  row.names(cd) <- cd[,\"number\"]
  +  n <- length(cd[,\"number\"])
  +  cd.r <- cd[r$region,]
  +  par(mfrow=c(2,3))
  +   hist(log(cd[,\"dtm_area\"]),col=\"lightgray\",border=\"lightgray\",
  +        main=\"\",xlab=\"AREA\",axes=FALSE)
  +   axis(1,at=c(log(1),log(10),log(100),log(1000),log(10000)),label=c(\"1\",\"10\",\"100\",\"1000\",\"10000\"))
  +   axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  +   box()
  +   points(cbind(log(cd.r[-1,\"dtm_area\"]),0),pch=19,cex=.7)
  +   points(cbind(log(cd.r[1,\"dtm_area\"]),0),pch=4,cex=2,lwd=2)
  +
  +   hist(cd[,\"saar\"],col=\"lightgray\",border=\"lightgray\",
  +        main=\"\",xlab=\"SAAR\",axes=FALSE)
  +   axis(1)
  +   axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  +   box()
  +   points(cbind(cd.r[-1,\"saar\"],0),pch=19,cex=.7)
  +   points(cbind(cd.r[1,\"saar\"],0),pch=4,cex=2,lwd=2)
  +
  +   hist(cd[,\"bfihost\"],col=\"lightgray\",border=\"lightgray\",
  +        main=\"\",xlab=\"BFIHOST\",axes=FALSE)
  +   axis(1)
  +   axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  +   box()
  +   points(cbind(cd.r[-1,\"bfihost\"],0),pch=19,cex=.7)
  +   points(cbind(cd.r[1,\"bfihost\"],0),pch=4,cex=2,lwd=2)
  +
  +   hist(cd[,\"farl\"],col=\"lightgray\",border=\"lightgray\",
  +        main=\"\",xlab=\"FARL\",axes=FALSE)
  +   axis(1)
  +   axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  +   box()
  +   points(cbind(cd.r[-1,\"farl\"],0),pch=19,cex=.7)
  +   points(cbind(cd.r[1,\"farl\"],0),pch=4,cex=2,lwd=2)
  +
  +   hist(cd[,\"propwet\"],col=\"lightgray\",border=\"lightgray\",
  +        main=\"\",xlab=\"PROPWET\",axes=FALSE)
  +   axis(1)
  +   axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  +   box()
  +   points(cbind(cd.r[-1,\"propwet\"],0),pch=19,cex=.7)
  +   points(cbind(cd.r[1,\"propwet\"],0),pch=4,cex=2,lwd=2)
  +
  +   hist(cd[,\"urbext1990\"],col=\"lightgray\",border=\"lightgray\",
  +        main=\"\",xlab=\"URBEXT\",axes=FALSE)
  +   axis(1)
  +   axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  +   box()
  +   points(cbind(cd.r[-1,\"urbext1990\"],0),pch=19,cex=.7)
  +   points(cbind(cd.r[1,\"urbext1990\"],0),pch=4,cex=2,lwd=2)
  +  par(mfrow=c(1,1))
  +  title(main=x,cex.main=1,font.main=1)
  + }
  > prova40009 <- roi.st.year(roi.cd[\"40009\",],roi.cd,row.names(roi.cd),am696[,\"am\"],
  +                           am696[,\"number\"],test=\"HW\",station.year=500,Nsim=500)
  > figure16.9a(\"40009\",prova40009,cd696)
  \n")

  figure16.9a <- function (x,r,cd) {
   # x = station of interest (e.g. "28018")
   # r = output of roi.st.year()

   if(!r$region[1]==x) r$region <- c(x,r$region)
   row.names(cd) <- cd[,"number"]
   n <- length(cd[,"number"])
   cd.r <- cd[r$region,]
   par(mfrow=c(2,3))
    hist(log(cd[,"dtm_area"]),col="lightgray",border="lightgray",
         main="",xlab="AREA",axes=FALSE)
    axis(1,at=c(log(1),log(10),log(100),log(1000),log(10000)),label=c("1","10","100","1000","10000"))
    axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
    box()
    points(cbind(log(cd.r[-1,"dtm_area"]),0),pch=19,cex=.7)
    points(cbind(log(cd.r[1,"dtm_area"]),0),pch=4,cex=2,lwd=2)
 
    hist(cd[,"saar"],col="lightgray",border="lightgray",
         main="",xlab="SAAR",axes=FALSE)
    axis(1)
    axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
    box()
    points(cbind(cd.r[-1,"saar"],0),pch=19,cex=.7)
    points(cbind(cd.r[1,"saar"],0),pch=4,cex=2,lwd=2)
 
    hist(cd[,"bfihost"],col="lightgray",border="lightgray",
         main="",xlab="BFIHOST",axes=FALSE)
    axis(1)
    axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
    box()
    points(cbind(cd.r[-1,"bfihost"],0),pch=19,cex=.7)
    points(cbind(cd.r[1,"bfihost"],0),pch=4,cex=2,lwd=2)
 
    hist(cd[,"farl"],col="lightgray",border="lightgray",
         main="",xlab="FARL",axes=FALSE)
    axis(1)
    axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
    box()
    points(cbind(cd.r[-1,"farl"],0),pch=19,cex=.7)
    points(cbind(cd.r[1,"farl"],0),pch=4,cex=2,lwd=2)
 
    hist(cd[,"propwet"],col="lightgray",border="lightgray",
         main="",xlab="PROPWET",axes=FALSE)
    axis(1)
    axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
    box()
    points(cbind(cd.r[-1,"propwet"],0),pch=19,cex=.7)
    points(cbind(cd.r[1,"propwet"],0),pch=4,cex=2,lwd=2)
 
    hist(cd[,"urbext1990"],col="lightgray",border="lightgray",
         main="",xlab="URBEXT",axes=FALSE)
    axis(1)
    axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
    box()
    points(cbind(cd.r[-1,"urbext1990"],0),pch=19,cex=.7)
    points(cbind(cd.r[1,"urbext1990"],0),pch=4,cex=2,lwd=2)
   par(mfrow=c(1,1)) 
   title(main=x,cex.main=1,font.main=1)
  }
  prova40009 <- roi.st.year(roi.cd["40009",],roi.cd,row.names(roi.cd),am696[,"am"],
                            am696[,"number"],test="HW",station.year=500,Nsim=500)
  figure16.9a("40009",prova40009,cd696)






  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Figure 16.9 pag.174 (2nd part), FEH Vol.3:
   there are differences because:
   - I plot the empirical growth curves;
   - site 40009 in FEH book has 14 data, while I have 25;
   - book uses POT for the polar plot, I only use annual maximum.\n
  > figure16.9b <- function (x,r,am,cd) {
  +  # x = station of interest (e.g. \"28018\")
  +  # r = output of roi.st.year()
  +
  +  row.names(cd) <- cd[,\"number\"]
  +  n <- length(cd[,\"number\"])
  +  cd.r <- cd[r$region,]
  +  cd.x <- cd[x,]
  +  fac <- factor(am[,\"number\"],levels=cd.r[,\"number\"])
  +  am.r <- am[!is.na(fac),]
  +  fac <- factor(am[,\"number\"],levels=x)
  +  am.x <- am[!is.na(fac),]
  +  am.xr <- rbind(am.x,am.r)
  +  QMED.r <- tapply(am.r[,4],am.r[,1],median)
  +  QMED.x <- median(am.x[,4])
  +  am.r.adim <- am.r; am.r.adim[,4] <- am.r[,4]/unsplit(QMED.r,am.r[,1])
  +  am.x.adim <- am.x; am.x.adim[,4] <- am.x[,4]/QMED.x
  +  lcv <- tapply(am[,4],am[,1],LCV)
  +  lca <- tapply(am[,4],am[,1],LCA)
  +  lkur <- tapply(am[,4],am[,1],Lkur)
  +  lcv.r <- tapply(am.r[,4],am.r[,1],LCV)
  +  lca.r <- tapply(am.r[,4],am.r[,1],LCA)
  +  lkur.r <- tapply(am.r[,4],am.r[,1],Lkur)
  +  lcv.x <- LCV(am.x[,4])
  +  lca.x <- LCA(am.x[,4])
  +  lkur.x <- Lkur(am.x[,4])
  +  days <- as.numeric(format(as.Date(am[,2]),\"%j\"))
  +  days.r <- as.numeric(format(as.Date(am.r[,2]),\"%j\"))
  +  days.x <- as.numeric(format(as.Date(am.x[,2]),\"%j\"))
  +
  +  par(mfrow=c(2,3))
  +   lognormplot(am.r.adim[,4],line=FALSE,xlab=\"Q/QMED\",type=\"n\")
  +   for(i in r$region) {
  +    x <- am.r.adim[am.r.adim[,1]==i,4]
  +    normpoints(x,type=\"l\",col=\"gray\")
  +   }
  +   normpoints(am.r.adim[,4],type=\"l\",lwd=2)
  +   normpoints(am.x.adim[,4],type=\"l\",col=2,lwd=2)
  +
  +   plot(lca,lcv,pch=\".\",cex=2)
  +   points(lca.r,lcv.r,pch=19)
  +   points(lca.x,lcv.x,pch=4,cex=2,lwd=2)
  +
  +   plot(lca,lkur,pch=\".\",cex=2)
  +   points(lca.r,lkur.r,pch=19)
  +   points(lca.x,lkur.x,pch=4,cex=2,lwd=2)
  +
  +   plot(cd[c(\"ihdtm_ngr_x\",\"ihdtm_ngr_y\")],pch=\".\",cex=2,xlab=\"\",ylab=\"\",axes=FALSE)
  +   points(cd.r[c(\"ihdtm_ngr_x\",\"ihdtm_ngr_y\")],pch=19)
  +   points(cd.x[c(\"ihdtm_ngr_x\",\"ihdtm_ngr_y\")],pch=4,cex=2,lwd=2)
  +
  +   consistencyplot (am.r[,3],am.r[,1])
  +
  +   dummy <- seq(0,2*pi,length=100)
  +   plot(cos(dummy),sin(dummy),type=\"l\",xlab=\"\",ylab=\"\",axes=FALSE)
  +   abline(h=0,lty=3); abline(v=0,lty=3)
  +   radd <- days*pi/180
  +   XFLOOD <- tapply(cos(radd),am[,1],mean)
  +   YFLOOD <- tapply(sin(radd),am[,1],mean)
  +   points(XFLOOD,YFLOOD,pch=\".\",cex=2)
  +   radd <- days.r*pi/180
  +   XFLOOD <- tapply(cos(radd),am.r[,1],mean)
  +   YFLOOD <- tapply(sin(radd),am.r[,1],mean)
  +   points(XFLOOD,YFLOOD,pch=19,cex=1)
  +   radd <- days.x*pi/180
  +   XFLOOD <- tapply(cos(radd),am.x[,1],mean)
  +   YFLOOD <- tapply(sin(radd),am.x[,1],mean)
  +   points(XFLOOD,YFLOOD,pch=4,cex=2,lwd=2)
  +   axis(1,at=0,label=\"Oct 1\")
  +   axis(2,at=0,label=\"Jul 1\")
  +   axis(3,at=0,label=\"Apr 1\")
  +   axis(4,at=0,label=\"Jan 1\")
  +  par(mfrow=c(1,1))
  +  title(main=x,cex.main=1,font.main=1)
  + }
  > figure16.9b(\"40009\",prova40009,am696,cd696)
  \n")

  figure16.9b <- function (x,r,am,cd) {
   # x = station of interest (e.g. "28018")
   # r = output of roi.st.year()

   row.names(cd) <- cd[,"number"]
   n <- length(cd[,"number"])
   cd.r <- cd[r$region,]
   cd.x <- cd[x,]
   fac <- factor(am[,"number"],levels=cd.r[,"number"])
   am.r <- am[!is.na(fac),]
   fac <- factor(am[,"number"],levels=x)
   am.x <- am[!is.na(fac),]
   am.xr <- rbind(am.x,am.r)
   QMED.r <- tapply(am.r[,4],am.r[,1],median)
   QMED.x <- median(am.x[,4])
   am.r.adim <- am.r; am.r.adim[,4] <- am.r[,4]/unsplit(QMED.r,am.r[,1])
   am.x.adim <- am.x; am.x.adim[,4] <- am.x[,4]/QMED.x
   lcv <- tapply(am[,4],am[,1],LCV)
   lca <- tapply(am[,4],am[,1],LCA)
   lkur <- tapply(am[,4],am[,1],Lkur)
   lcv.r <- tapply(am.r[,4],am.r[,1],LCV)
   lca.r <- tapply(am.r[,4],am.r[,1],LCA)
   lkur.r <- tapply(am.r[,4],am.r[,1],Lkur)
   lcv.x <- LCV(am.x[,4])
   lca.x <- LCA(am.x[,4])
   lkur.x <- Lkur(am.x[,4])
   days <- as.numeric(format(as.Date(am[,2]),"%j"))
   days.r <- as.numeric(format(as.Date(am.r[,2]),"%j"))
   days.x <- as.numeric(format(as.Date(am.x[,2]),"%j"))

   par(mfrow=c(2,3))
    lognormplot(am.r.adim[,4],line=FALSE,xlab="Q/QMED",type="n")
    for(i in r$region) {
     xxx <- am.r.adim[am.r.adim[,1]==i,4]
     normpoints(xxx,type="l",col="gray")
    }
    normpoints(am.r.adim[,4],type="l",lwd=2)
    normpoints(am.x.adim[,4],type="l",col=2,lwd=2)
 
    plot(lca,lcv,pch=".",cex=2)
    points(lca.r,lcv.r,pch=19)
    points(lca.x,lcv.x,pch=4,cex=2,lwd=2)
 
    plot(lca,lkur,pch=".",cex=2)
    points(lca.r,lkur.r,pch=19)
    points(lca.x,lkur.x,pch=4,cex=2,lwd=2)

    plot(cd[c("ihdtm_ngr_x","ihdtm_ngr_y")],pch=".",cex=2,xlab="",ylab="",axes=FALSE)
    points(cd.r[c("ihdtm_ngr_x","ihdtm_ngr_y")],pch=19)
    points(cd.x[c("ihdtm_ngr_x","ihdtm_ngr_y")],pch=4,cex=2,lwd=2)

    consistencyplot (am.r[,3],am.r[,1])

    dummy <- seq(0,2*pi,length=100)
    plot(cos(dummy),sin(dummy),type="l",xlab="",ylab="",axes=FALSE)
    abline(h=0,lty=3); abline(v=0,lty=3)
    radd <- days*pi/180
    XFLOOD <- tapply(cos(radd),am[,1],mean)
    YFLOOD <- tapply(sin(radd),am[,1],mean)
    points(XFLOOD,YFLOOD,pch=".",cex=2)
    radd <- days.r*pi/180
    XFLOOD <- tapply(cos(radd),am.r[,1],mean)
    YFLOOD <- tapply(sin(radd),am.r[,1],mean)
    points(XFLOOD,YFLOOD,pch=19,cex=1)
    radd <- days.x*pi/180
    XFLOOD <- tapply(cos(radd),am.x[,1],mean)
    YFLOOD <- tapply(sin(radd),am.x[,1],mean)
    points(XFLOOD,YFLOOD,pch=4,cex=2,lwd=2)
    axis(1,at=0,label="Oct 1")
    axis(2,at=0,label="Jul 1")
    axis(3,at=0,label="Apr 1")
    axis(4,at=0,label="Jan 1") 
   par(mfrow=c(1,1))
   title(main=x,cex.main=1,font.main=1)
  }
  figure16.9b("40009",prova40009,am696,cd696)


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Figure 16.9 pag.174 (2nd part), FEH Vol.3:
   there are differences because:
   - I plot the empirical growth curves;
   - book uses POT for the polar plot, I only use annual maximum.\n
  > prova45001 <- roi.st.year(roi.cd[\"45001\",],roi.cd,row.names(roi.cd),am696[,\"am\"],
  +                           am696[,\"number\"],test=\"HW\",station.year=250,Nsim=500)
  > figure16.9a(\"45001\",prova45001,cd696)
  > readline(\"Press Return to continue to the next graphs\")
  > figure16.9b(\"45001\",prova45001,am696,cd696)
  \n")
  
  prova45001 <- roi.st.year(roi.cd["45001",],roi.cd,row.names(roi.cd),am696[,"am"],
                            am696[,"number"],test="HW",station.year=250,Nsim=500)
  figure16.9a("45001",prova45001,cd696)
  readline("Press Return to continue to the next graphs")
  figure16.9b("45001",prova45001,am696,cd696)

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                                          THE END
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")
}






































exampleRFA04 <- function () {

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                           ANALYSIS OF FEH DATA: INDEX-VALUE
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")

  readline("\n
  Press Return to continue
  \n")

  cat("\n
  Data loading: \n
  > data(FEH1000)
  \n")

  data(FEH1000)

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Sites used in model development (pag.102 FEH Vol.3):
   area>0.5 km2
   digital catchment data available
   urbext<0.025\n
  > urbext <- cd[,\"urbext1990\"]
  > area <- cd[,\"dtm_area\"]
  > cd732 <- cd[(!is.nan(cd[,\"dtm_area\"]))&(urbext<0.025)&(area>0.5),]

  > fac <- factor(am[,\"number\"],levels=cd732[,\"number\"])
  > am732 <- am[!is.na(fac),]
  > nlevels(as.factor(am732[,\"number\"]))
  \n")

  urbext <- cd[,"urbext1990"]
  area <- cd[,"dtm_area"]
  cd732 <- cd[(!is.nan(cd[,"dtm_area"]))&(urbext<0.025)&(area>0.5),] # vs 687 - 728 of FEH

  fac <- factor(am[,"number"],levels=cd732[,"number"])
  am732 <- am[!is.na(fac),]
  nlevels(as.factor(am732[,"number"]))


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Index-flood = median:\n
  > QMED <- tapply(am732[,4],am732[,1],median)
  > lnQMED <- log(QMED)
  \n")

  QMED <- tapply(am732[,4],am732[,1],median)
  lnQMED <- log(QMED)

  cat("\n
  Catchment descriptors (fig. 13.1 pag. 104):\n
  > lnAREA <- log(cd732[,\"dtm_area\"])
  > lnDPLBAR <- log(cd732[,\"dplbar\"])
  > lnSPRHOST <- log(cd732[,\"sprhost\"])
  > lnBFIHOST <- log(cd732[,\"bfihost\"])
  > lnSAAR <- log(cd732[,\"saar\"])
  > lnRMED1 <- log(cd732[,\"rmed_1d\"])
  > lnALTBAR <- log(cd732[,\"altbar\"])
  > lnFARL <- log(cd732[,\"farl\"])

  > M <- data.frame(cbind(lnQMED,lnAREA,lnDPLBAR,lnSPRHOST,lnBFIHOST,lnSAAR,lnRMED1,lnALTBAR,lnFARL))
  > print(cor(M))
  > plot(M,pch=\".\",cex=2)
  \n")

  lnAREA <- log(cd732[,"dtm_area"])
  lnDPLBAR <- log(cd732[,"dplbar"])
  lnSPRHOST <- log(cd732[,"sprhost"])
  lnBFIHOST <- log(cd732[,"bfihost"])
  lnSAAR <- log(cd732[,"saar"])
  lnRMED1 <- log(cd732[,"rmed_1d"])
  #lnNWET <- log(cd732[,""])
  lnALTBAR <- log(cd732[,"altbar"])
  lnFARL <- log(cd732[,"farl"])

  M <- data.frame(cbind(lnQMED,lnAREA,lnDPLBAR,lnSPRHOST,lnBFIHOST,lnSAAR,lnRMED1,lnALTBAR,lnFARL))
  print(cor(M))
  plot(M,pch=".",cex=2)


  cat("\n
  Additional variables (pag. 105):\n
  > RESHOST <- cd732[,\"bfihost\"] + 1.30*(cd732[,\"sprhost\"]/100)-0.987
  > lnAREAsq <- lnAREA^2
  > lnSAARsq <- lnSAAR^2

  > M <- data.frame(cbind(M,RESHOST,lnAREAsq,lnSAARsq))
  \n")

  RESHOST <- cd732[,"bfihost"] + 1.30*(cd732[,"sprhost"]/100)-0.987
  lnAREAsq <- lnAREA^2
  lnSAARsq <- lnSAAR^2

  M <- data.frame(cbind(M,RESHOST,lnAREAsq,lnSAARsq))


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Ordinary Least Square models:\n
   create a function using the 'leaps' function of package 'subselect'
   to perform all-possible-regressions: \n
  >   bestregressions <- function(dip,ind) {
  +    Y <- as.numeric(dip)
  +    X <- ind
  +    Sy <- var(Y)
  +    Sx <- var(X)
  +    Sxy <- var(X,Y)
  +    Dm.mat <- Sx
  +    Dm.H <- Sxy %*% t(Sxy)/Sy
  +    require(subselect)
  +    Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=6, H=Dm.H, r=1, nsol=3)
  +    Dm.leaps
  +    for(i in 6:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:6),i]])}}
  +   }
  >
  >   bestregressions(M[,1],M[,-1])
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"lnFARL\"    \"RESHOST\"   \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnBFIHOST\" \"lnSAAR\"    \"lnFARL\"    \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"lnRMED1\"   \"lnFARL\"    \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"lnFARL\"    \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"lnFARL\"    \"RESHOST\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"RESHOST\"   \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"lnFARL\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"    \"RESHOST\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAAR\"
  [1] \"lnAREA\"    \"lnSPRHOST\" \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnBFIHOST\" \"lnSAAR\"
  [1] \"lnAREA\" \"lnSAAR\"
  [1] \"lnAREA\"   \"lnSAARsq\"
  [1] \"lnAREA\"    \"lnSPRHOST\"
  [1] \"lnAREAsq\"
  [1] \"lnAREA\"
  [1] \"lnDPLBAR\"
  \n")


  #bestregressions <- function(dip,ind) {
  # Y <- as.numeric(dip)
  # X <- ind
  # Sy <- var(Y)
  # Sx <- var(X)
  # Sxy <- var(X,Y)
  # Dm.mat <- Sx
  # Dm.H <- Sxy %*% t(Sxy)/Sy
  # require(subselect)
  # Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=6, H=Dm.H, r=1, nsol=3)
  # Dm.leaps
  # for(i in 6:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:6),i]])}}
  #}

  #bestregressions(M[,1],M[,-1])



  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Ordinary Least Square models (graphics and statistics):\n
  > graphics.lm <- function (regr) {
  +  par(mfrow=c(2,2))
  +  plot(regr$fitted.values,regr$residuals,pch=\".\",cex=3,xlab=\"lnQMED Fitted\",ylab=\"lnQMED Residuals\")
  +  abline(0,0,lty=3)
  +  normplot(regr$residuals,pch=\".\",cex=3,xlab=\"lnQMED Residuals\")
  +
  +  plot(regr$fitted.values,lnQMED,pch=\".\",cex=3,xlab=\"lnQMED Originals\",ylab=\"lnQMED Fitted\")
  +  abline(0,1,lty=3)
  +  intervals <- predinterval.lm(regr)
  +  intervals <- intervals[order(intervals[,1]),]
  +  lines(intervals[,c(1,2)],lty=2)
  +  lines(intervals[,c(1,3)],lty=2)
  +  Rsq <- round(R2.lm(regr),3)
  +  rmse <- signif(RMSE.lm(regr),3)
  +  rmsep <- signif(RMSEP(lnQMED,regr$fitted.values)*100,3)
  +  mtext(paste(\"R2 = \",Rsq),3,-1.5,adj=0.05)
  +  mtext(paste(\"RMSE = \",rmse),1,-2.5,adj=0.95)
  +  mtext(paste(\"RMSEP = \",rmsep,\"%\"),1,-1.5,adj=0.95)
  +
  +  plot(exp(regr$fitted.values),exp(lnQMED),pch=\".\",cex=3,xlab=\"QMED Originals\",ylab=\"QMED Fitted\")
  +  abline(0,1,lty=3)
  +  lines(exp(intervals[,c(1,2)]),lty=2)
  +  lines(exp(intervals[,c(1,3)]),lty=2)
  +  Rsq <- round(R2(exp(lnQMED),exp(regr$fitted.values)),3)
  +  rmse <- signif(RMSE(exp(lnQMED),exp(regr$fitted.values)),3)
  +  rmsep <- signif(RMSEP(exp(lnQMED),exp(regr$fitted.values))*100,3)
  +  mtext(paste(\"R2 = \",Rsq),3,-1.5,adj=0.05)
  +  mtext(paste(\"RMSE = \",rmse),1,-2.5,adj=0.95)
  +  mtext(paste(\"RMSEP = \",rmsep,\"%\"),1,-1.5,adj=0.95)
  +
  +  par(mfrow=c(1,1))
  +  title(main=paste(names(regr$coefficients)[-1], collapse=\", \"),cex.main=.7,font.main=1)
  + }

  > graphics.lm(lm(lnQMED ~ lnDPLBAR))
  > graphics.lm(lm(lnQMED ~ lnAREA))
  > graphics.lm(lm(lnQMED ~ lnAREAsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSAARsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSAAR))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnBFIHOST + lnSAAR))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAARsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + RESHOST))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnSAARsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + RESHOST + lnSAARsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + RESHOST))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + lnSAARsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnRMED1 + lnFARL + lnSAARsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnBFIHOST + lnSAAR + lnFARL + lnSAARsq))
  > graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + RESHOST + lnSAARsq))
  \n")


  graphics.lm <- function (regr) {
   par(mfrow=c(2,2), cex=.7)
   plot(regr$fitted.values,regr$residuals,pch=".",cex=3,xlab="lnQMED Fitted",ylab="lnQMED Residuals")
   abline(0,0,lty=3)
   normplot(regr$residuals,pch=".",cex=3,xlab="lnQMED Residuals")
 
   plot(regr$fitted.values,lnQMED,pch=".",cex=3,xlab="lnQMED Originals",ylab="lnQMED Fitted")
   abline(0,1,lty=3)
   intervals <- predinterval.lm(regr)
   intervals <- intervals[order(intervals[,1]),]
   lines(intervals[,c(1,2)],lty=2)
   lines(intervals[,c(1,3)],lty=2)
   Rsq <- signif(R2.lm(regr),3)
   rmse <- signif(RMSE.lm(regr),3)
   rmsep <- signif(RMSEP(lnQMED,regr$fitted.values)*100,3)
   mtext(paste("R2 = ",Rsq),3,-1.5,adj=0.05,cex=.7)
   mtext(paste("RMSE = ",rmse),1,-2.5,adj=0.95,cex=.7)
   mtext(paste("RMSEP = ",rmsep,"%"),1,-1.5,adj=0.95,cex=.7)

   plot(exp(regr$fitted.values),exp(lnQMED),pch=".",cex=3,xlab="QMED Originals",ylab="QMED Fitted")
   abline(0,1,lty=3)
   lines(exp(intervals[,c(1,2)]),lty=2)
   lines(exp(intervals[,c(1,3)]),lty=2)
   Rsq <- signif(R2(exp(lnQMED),exp(regr$fitted.values)),3)
   rmse <- signif(RMSE(exp(lnQMED),exp(regr$fitted.values)),3)
   rmsep <- signif(RMSEP(exp(lnQMED),exp(regr$fitted.values))*100,3)
   mtext(paste("R2 = ",Rsq),3,-1.5,adj=0.05,cex=.7)
   mtext(paste("RMSE = ",rmse),1,-2.5,adj=0.95,cex=.7)
   mtext(paste("RMSEP = ",rmsep,"%"),1,-1.5,adj=0.95,cex=.7)

   par(mfrow=c(1,1),cex=1)
   title(main=paste(names(regr$coefficients)[-1], collapse=", "),cex.main=.7,font.main=1)
  }

  graphics.lm(lm(lnQMED ~ lnDPLBAR))
  readline("Press Return to continue to the next graphs") 
  graphics.lm(lm(lnQMED ~ lnAREA))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREAsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSAARsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSAAR))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnBFIHOST + lnSAAR))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAARsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + RESHOST))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnSAARsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + RESHOST + lnSAARsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + RESHOST))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + lnSAARsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnRMED1 + lnFARL + lnSAARsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnBFIHOST + lnSAAR + lnFARL + lnSAARsq))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + RESHOST + lnSAARsq))



  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                                          THE END
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")
}






































exampleRFA05 <- function () {

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                        ANALYSIS OF FEH DATA: CLASSIFICATION-VARIABLES
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")

  readline("\n
  Press Return to continue
  \n")

  cat("\n
  Data loading: \n
  > data(FEH1000)
  \n")

  data(FEH1000)


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Criteria used in the FEH to choose stations for pooling groups:
   n>7
   area, saar and bfihost are known
   urbext<0.025
   area>0.5 \n
  > n <- tapply(am[,3],am[,1],length)
  > urbext <- cd[,\"urbext1990\"]
  > area <- cd[,\"dtm_area\"]
  > cd696 <- cd[(!is.nan(cd[,\"dtm_area\"]))&(!is.nan(cd[,\"saar\"]))&(!is.nan(cd[,\"bfihost\"]))&(n>7)&(urbext<0.025)&(area>0.5),]

  > fac <- factor(am[,\"number\"],levels=cd696[,\"number\"])
  > am696 <- am[!is.na(fac),]
  > nlevels(as.factor(am696[,\"number\"]))
  \n")

  n <- tapply(am[,4],am[,1],length)
  urbext <- cd[,"urbext1990"]
  area <- cd[,"dtm_area"]
  cd696 <- cd[(!is.nan(cd[,"dtm_area"]))&(!is.nan(cd[,"saar"]))&(!is.nan(cd[,"bfihost"]))&(n>7)&(urbext<0.025)&(area>0.5),]

  fac <- factor(am[,"number"],levels=cd696[,"number"])
  am696 <- am[!is.na(fac),]
  nlevels(as.factor(am696[,"number"]))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  FEH classification variables:
   lnAREA, lnSAAR, BFIHOST. \n
  > lcv <- tapply(am696[,4],am696[,1],LCV)
  > lca <- tapply(am696[,4],am696[,1],LCA)

  > lnAREA <- log(cd696[,\"dtm_area\"])
  > lnSAAR <- log(cd696[,\"saar\"])
  > BFIHOST <- cd696[,\"bfihost\"]

  > lcvVSclassvarFEH <- lm(lcv ~ lnAREA + lnSAAR + BFIHOST)
  > lcaVSclassvarFEH <- lm(lca ~ lnAREA + lnSAAR + BFIHOST)
  
  > graphics.lm <- function (regr) {
  +  par(mfrow=c(1,2), cex=1)
  +   plot(regr$fitted.values,regr$model[,1],pch=\".\",cex=3,xlab=\"Originals\",ylab=\"Fitted\")
  +   abline(0,1,lty=3)
  +   intervals <- predinterval.lm(regr)
  +   intervals <- intervals[order(intervals[,1]),]
  +   lines(intervals[,c(1,2)],lty=2)
  +   lines(intervals[,c(1,3)],lty=2)
  +   Rsq <- signif(R2.lm(regr),3)
  +   rmse <- signif(RMSE.lm(regr),3)
  +   rmsep <- signif(RMSEP(regr$model[,1],regr$fitted.values)*100,3)
  +   mtext(paste(\"R2 = \",Rsq[1]),3,-1.5,adj=0.05,cex=1)
  +   mtext(paste(\"R2adj = \",Rsq[2]),3,-2.5,adj=0.05,cex=1)
  +   mtext(paste(\"RMSE = \",rmse),1,-2.5,adj=0.95,cex=1)
  +   mtext(paste(\"RMSEP = \",rmsep,\"%\"),1,-1.5,adj=0.95,cex=1)
  +
  +   normplot(regr$residuals,pch=\".\",cex=3,xlab=\"Residuals\")
  +  par(mfrow=c(1,1),cex=1)
  +  title(main=paste(names(regr$model[1]),\"~\",paste(names(regr$model[-1]), collapse=\", \")),cex.main=1,font.main=1)
  + }

  > graphics.lm(lcvVSclassvarFEH)
  > graphics.lm(lcaVSclassvarFEH)
  \n")

  lcv <- tapply(am696[,4],am696[,1],LCV)
  lca <- tapply(am696[,4],am696[,1],LCA)

  lnAREA <- log(cd696[,"dtm_area"])
  lnSAAR <- log(cd696[,"saar"])
  BFIHOST <- cd696[,"bfihost"]

  lcvVSclassvarFEH <- lm(lcv ~ lnAREA + lnSAAR + BFIHOST)
  lcaVSclassvarFEH <- lm(lca ~ lnAREA + lnSAAR + BFIHOST)

  graphics.lm <- function (regr) {
   par(mfrow=c(1,2), cex=1)
    plot(regr$fitted.values,regr$model[,1],pch=".",cex=3,xlab="Originals",ylab="Fitted")
    abline(0,1,lty=3)
    intervals <- predinterval.lm(regr)
    intervals <- intervals[order(intervals[,1]),]
    lines(intervals[,c(1,2)],lty=2)
    lines(intervals[,c(1,3)],lty=2)
    Rsq <- signif(R2.lm(regr),3)
    rmse <- signif(RMSE.lm(regr),3)
    rmsep <- signif(RMSEP(regr$model[,1],regr$fitted.values)*100,3)
    mtext(paste("R2 = ",Rsq[1]),3,-1.5,adj=0.05,cex=1)
    mtext(paste("R2adj = ",Rsq[2]),3,-2.5,adj=0.05,cex=1)
    mtext(paste("RMSE = ",rmse),1,-2.5,adj=0.95,cex=1)
    mtext(paste("RMSEP = ",rmsep,"%"),1,-1.5,adj=0.95,cex=1)

    normplot(regr$residuals,pch=".",cex=3,xlab="Residuals")
   par(mfrow=c(1,1),cex=1)
   title(main=paste(names(regr$model[1]),"~",paste(names(regr$model[-1]), collapse=", ")),cex.main=1,font.main=1)
  }

  graphics.lm(lcvVSclassvarFEH)
  readline("Press Return to continue to the next graphs")
  graphics.lm(lcaVSclassvarFEH)



  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Alternative classification variables:
  > NGRX <- cd696[,\"ihdtm_ngr_x\"]
  > NGRY <- cd696[,\"ihdtm_ngr_y\"]
  > AREA <- cd696[,\"dtm_area\"]
  > SAAR <- cd696[,\"saar\"]
  > SPRHOST <- cd696[,\"sprhost\"]
  > FARL <- cd696[,\"farl\"]
  > RMED1D <- cd696[,\"rmed_1d\"]
  > RMED2D <- cd696[,\"rmed_2d\"]
  > RMED1H <- cd696[,\"rmed_1h\"]
  > SMDBAR <- cd696[,\"smdbar\"]
  > PROPWET <- cd696[,\"propwet\"]
  > LDP <- cd696[,\"ldp\"]
  > DPLBAR <- cd696[,\"dplbar\"]
  > ALTBAR <- cd696[,\"altbar\"]
  > DPSBAR <- cd696[,\"dpsbar\"]
  > ASPBAR <- cd696[,\"aspbar\"]
  > ASPVAR <- cd696[,\"aspvar\"]
  > URBEXT <- cd696[,\"urbext1990\"]

  > potCVnat <- data.frame(cbind(NGRX,NGRY,AREA,SAAR,BFIHOST,SPRHOST,FARL,RMED1D,RMED2D,RMED1H,
  >                              SMDBAR,PROPWET,LDP,DPLBAR,ALTBAR,DPSBAR,ASPBAR,ASPVAR,URBEXT))
  > potCVln <- log(potCVnat[-c(17,19)]); names(potCVln) <- paste(\"ln\",names(potCVln),sep=\"\")
  > potCV <- cbind(potCVnat,potCVln)
  > bestregressions <- function(dip,ind) {
  >  Y <- as.numeric(dip)
  >  X <- ind
  >  Sy <- var(Y)
  >  Sx <- var(X)
  >  Sxy <- var(X,Y)
  >  Dm.mat <- Sx
  >  Dm.H <- Sxy %*% t(Sxy)/Sy
  >  require(subselect)
  >  Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=6, H=Dm.H, r=1, nsol=3)
  >  Dm.leaps
  >  for(i in 6:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:6),i]])}}
  > }

  > bestregressions(lcv,potCVnat[-19])
  > bestregressions(lcv,potCVln)
  > bestregressions(lcv,cbind(potCVnat[-19],potCVln[-c(1,4:5,7:13)]))
  > bestregressions(lcv,cbind(potCVnat[-c(1,4:5,7:13,19)],potCVln))
  > bestregressions(lcv,cbind(potCVnat[-c(1:2,19)],potCVln))
  [1] \"RMED1D\"   \"ALTBAR\"   \"lnAREA\"   \"lnSAAR\"   \"lnRMED2D\" \"lnALTBAR\"
  [1] \"RMED1D\"   \"DPLBAR\"   \"ALTBAR\"   \"lnSAAR\"   \"lnRMED2D\" \"lnALTBAR\"
  [1] \"RMED1D\"   \"LDP\"      \"ALTBAR\"   \"lnSAAR\"   \"lnRMED2D\" \"lnALTBAR\"
  [1] \"DPLBAR\"   \"ALTBAR\"   \"lnSAAR\"   \"lnRMED2D\" \"lnALTBAR\"
  [1] \"SAAR\"     \"RMED1D\"   \"lnAREA\"   \"lnSAAR\"   \"lnRMED2D\"
  [1] \"RMED1D\"   \"DPSBAR\"   \"lnAREA\"   \"lnSAAR\"   \"lnRMED2D\"
  [1] \"lnAREA\"   \"lnSAAR\"   \"lnRMED1D\" \"lnRMED2D\"
  [1] \"RMED1D\"   \"lnAREA\"   \"lnSAAR\"   \"lnRMED2D\"
  [1] \"DPLBAR\"   \"ALTBAR\"   \"lnSAAR\"   \"lnRMED2D\"
  [1] \"DPLBAR\"   \"lnSAAR\"   \"lnRMED2D\"
  [1] \"lnAREA\"   \"lnSAAR\"   \"lnRMED2D\"
  [1] \"LDP\"      \"lnSAAR\"   \"lnRMED2D\"
  [1] \"lnAREA\" \"lnSAAR\"
  [1] \"lnSAAR\"   \"lnDPLBAR\"
  [1] \"lnSAAR\" \"lnLDP\"
  [1] \"lnSAAR\"
  [1] \"SMDBAR\"
  [1] \"lnPROPWET\"

  > names(cbind(potCVnat[-c(1:2,19)],potCVln))

  > graphics.lm(lm(lcv ~ lnAREA + lnSAAR + lnRMED2D,data=potCV))
  > readline(\"Press Return to continue to the next graphs\")
  > graphics.lm(lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV))
  > readline(\"Press Return to continue to the next graphs\")
  > graphics.lm(lm(lcv ~ DPLBAR + ALTBAR + lnSAAR + lnRMED2D,data=potCV))
  > readline(\"Press Return to continue to the next graphs\")
  > graphics.lm(lm(lcv ~ RMED1D + ALTBAR + lnAREA + lnSAAR + lnRMED2D + lnALTBAR,data=potCV))

  > prt.lm(lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV))
  > prt.lm(lm(lcv ~ lnAREA + lnSAAR + lnRMED1D + lnRMED2D,data=potCV))
  > prt.lm(lm(lcv ~ DPLBAR + ALTBAR + lnSAAR + lnRMED2D + lnALTBAR,data=potCV))
  > prt.lm(lm(lcv ~ RMED1D + ALTBAR + lnAREA + lnSAAR + lnRMED2D + lnALTBAR,data=potCV))
  \n")

  NGRX <- cd696[,"ihdtm_ngr_x"]
  NGRY <- cd696[,"ihdtm_ngr_y"]
  AREA <- cd696[,"dtm_area"]
  SAAR <- cd696[,"saar"]
  SPRHOST <- cd696[,"sprhost"]
  FARL <- cd696[,"farl"]
  RMED1D <- cd696[,"rmed_1d"]
  RMED2D <- cd696[,"rmed_2d"]
  RMED1H <- cd696[,"rmed_1h"]
  SMDBAR <- cd696[,"smdbar"]
  PROPWET <- cd696[,"propwet"]
  LDP <- cd696[,"ldp"]
  DPLBAR <- cd696[,"dplbar"]
  ALTBAR <- cd696[,"altbar"]
  DPSBAR <- cd696[,"dpsbar"]
  ASPBAR <- cd696[,"aspbar"]
  ASPVAR <- cd696[,"aspvar"]
  URBEXT <- cd696[,"urbext1990"]

  potCVnat <- data.frame(cbind(NGRX,NGRY,AREA,SAAR,BFIHOST,SPRHOST,FARL,RMED1D,RMED2D,RMED1H,
                               SMDBAR,PROPWET,LDP,DPLBAR,ALTBAR,DPSBAR,ASPBAR,ASPVAR,URBEXT))
  potCVln <- log(potCVnat[-c(17,19)]); names(potCVln) <- paste("ln",names(potCVln),sep="")
  potCV <- cbind(potCVnat,potCVln)

  # bestregressions <- function(dip,ind) {
  #  Y <- as.numeric(dip)
  #  X <- ind
  #  Sy <- var(Y)
  #  Sx <- var(X)
  #  Sxy <- var(X,Y)
  #  Dm.mat <- Sx
  #  Dm.H <- Sxy %*% t(Sxy)/Sy
  #  require(subselect)
  #  Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=6, H=Dm.H, r=1, nsol=3)
  #  Dm.leaps
  #  for(i in 6:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:6),i]])}}
  # }

  # bestregressions(lcv,potCVnat[-19])
  # bestregressions(lcv,potCVln)
  # bestregressions(lcv,cbind(potCVnat[-19],potCVln[-c(1,4:5,7:13)]))
  # bestregressions(lcv,cbind(potCVnat[-c(1,4:5,7:13,19)],potCVln))
  # bestregressions(lcv,cbind(potCVnat[-c(1:2,19)],potCVln))
  names(cbind(potCVnat[-c(1:2,19)],potCVln))

  graphics.lm(lm(lcv ~ lnAREA + lnSAAR + lnRMED2D,data=potCV))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lcv ~ DPLBAR + ALTBAR + lnSAAR + lnRMED2D,data=potCV))
  readline("Press Return to continue to the next graphs")
  graphics.lm(lm(lcv ~ RMED1D + ALTBAR + lnAREA + lnSAAR + lnRMED2D + lnALTBAR,data=potCV))

  print(prt.lm(lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV)))
  print(prt.lm(lm(lcv ~ lnAREA + lnSAAR + lnRMED1D + lnRMED2D,data=potCV)))
  print(prt.lm(lm(lcv ~ DPLBAR + ALTBAR + lnSAAR + lnRMED2D + lnALTBAR,data=potCV)))
  print(prt.lm(lm(lcv ~ RMED1D + ALTBAR + lnAREA + lnSAAR + lnRMED2D + lnALTBAR,data=potCV)))


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Distance matrices and Mantel test approach: best models\n
  > distancesAD <- AD.dist(am696[,4],am696[,1])     # if run, it takes 1 hour
  > distancesCVnat <- data.frame(apply(potCVnat,2,dist))
  > distancesCVln <- data.frame(apply(potCVln,2,dist))
  > distancesCV <- data.frame(apply(cbind(potCVnat[-c(1:2,19)],potCVln),2,dist))

  > bestregressions(as.numeric(distancesAD),distancesCVnat)
  > bestregressions(as.numeric(distancesAD),distancesCVln)
  > bestregressions(as.numeric(distancesAD),distancesCV)
  [1] \"FARL\"     \"SMDBAR\"   \"DPLBAR\"   \"lnNGRY\"   \"lnSAAR\"   \"lnRMED2D\"
  [1] \"FARL\"   \"RMED2D\" \"DPLBAR\" \"lnNGRX\" \"lnNGRY\" \"lnSAAR\"
  [1] \"FARL\"     \"RMED2D\"   \"lnNGRX\"   \"lnNGRY\"   \"lnSAAR\"   \"lnDPLBAR\"
  [1] \"FARL\"     \"RMED2D\"   \"lnNGRY\"   \"lnSAAR\"   \"lnDPLBAR\"
  [1] \"FARL\"   \"RMED2D\" \"DPLBAR\" \"lnNGRY\" \"lnSAAR\"
  [1] \"RMED2D\"   \"lnNGRY\"   \"lnSAAR\"   \"lnFARL\"   \"lnDPLBAR\"
  [1] \"FARL\"   \"RMED2D\" \"DPLBAR\" \"lnSAAR\"
  [1] \"FARL\"     \"RMED2D\"   \"lnSAAR\"   \"lnDPLBAR\"
  [1] \"RMED2D\" \"DPLBAR\" \"lnSAAR\" \"lnFARL\"
  [1] \"FARL\"   \"RMED2D\" \"lnSAAR\"
  [1] \"RMED2D\" \"lnSAAR\" \"lnFARL\"
  [1] \"RMED2D\" \"DPLBAR\" \"lnSAAR\"
  [1] \"RMED2D\" \"lnSAAR\"
  [1] \"lnSAAR\"   \"lnRMED2D\"
  [1] \"RMED1D\" \"lnSAAR\"
  [1] \"SMDBAR\"
  [1] \"lnSAAR\"
  [1] \"lnPROPWET\"
  \n")

  # distancesAD <- AD.dist(am696[,4],am696[,1])                                   # 1 hour
  # distancesCVnat <- data.frame(apply(potCVnat,2,dist))                          
  # distancesCVln <- data.frame(apply(potCVln,2,dist))                            
  # distancesCV <- data.frame(apply(cbind(potCVnat[-c(1:2,19)],potCVln),2,dist))  

  # bestregressions(as.numeric(distancesAD),distancesCVnat)  
  # bestregressions(as.numeric(distancesAD),distancesCVln)   
  # bestregressions(as.numeric(distancesAD),distancesCV)     
  

  # graphics.lm(lm(distancesAD ~ FARL + RMED2D + lnSAAR,data=distancesCV))
  # graphics.lm(lm(distancesAD ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=distancesCV))
  # graphics.lm(lm(distancesAD ~ DPLBAR + lnSAAR + lnRMED2D,data=distancesCV))

  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Distance matrices and Mantel test approach: mantel test\n
  > regr1CV <- lm(distancesAD ~ SMDBAR,data=distancesCV)
  > regr2CV <- lm(distancesAD ~ RMED2D + lnSAAR,data=distancesCV)
  > regr3CV <- lm(distancesAD ~ FARL + RMED2D + lnSAAR,data=distancesCV)
  > regr4CV <- lm(distancesAD ~ FARL + RMED2D + DPLBAR + lnSAAR,data=distancesCV)
  > regr5CV <- lm(distancesAD ~ FARL + RMED2D + lnNGRY + lnSAAR + lnDPLBAR,data=distancesCV)
  > regr6CV <- lm(distancesAD ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=distancesCV)

  > mantel.lm(regr1CV, Nperm=100) # if run, it takes some time
  P.SMDBAR
         1
  > mantel.lm(regr2CV, Nperm=100)
  P.RMED2D P.lnSAAR
         1        1
  > mantel.lm(regr3CV, Nperm=100)
  P.FARL P.RMED2D P.lnSAAR
       1        1        1
  > mantel.lm(regr4CV, Nperm=100)
  P.FARL P.RMED2D P.DPLBAR P.lnSAAR
       1        1        1        1
  > mantel.lm(regr5CV, Nperm=100)
  P.FARL   P.RMED2D   P.lnNGRY   P.lnSAAR P.lnDPLBAR
       1          1          1          1          1
  > mantel.lm(regr6CV, Nperm=100)
  P.FARL   P.SMDBAR   P.DPLBAR   P.lnNGRY   P.lnSAAR P.lnRMED2D
       1          1          1          1          1          1
  > mantel.lm(lm(distancesAD ~ RMED2D + RMED1D,data=distancesCV), Nperm=100)
  P.RMED2D P.RMED1D
      1.00     0.91
  > prt.lm(lm(lcv ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=potCV))
          FARL       SMDBAR       DPLBAR       lnNGRY       lnSAAR     lnRMED2D
  2.076164e-01 6.356903e-01 8.685326e-10 4.455110e-01 3.362767e-13 2.852258e-08
  \n")

  # regr1CV <- lm(distancesAD ~ SMDBAR,data=distancesCV)
  # regr2CV <- lm(distancesAD ~ RMED2D + lnSAAR,data=distancesCV)
  # regr3CV <- lm(distancesAD ~ FARL + RMED2D + lnSAAR,data=distancesCV)
  # regr4CV <- lm(distancesAD ~ FARL + RMED2D + DPLBAR + lnSAAR,data=distancesCV)
  # regr5CV <- lm(distancesAD ~ FARL + RMED2D + lnNGRY + lnSAAR + lnDPLBAR,data=distancesCV)
  # regr6CV <- lm(distancesAD ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=distancesCV)

  # mantel.lm(regr1CV, Nperm=100)   # 2 min 
  # mantel.lm(regr2CV, Nperm=100)   # 4 min
  # mantel.lm(regr3CV, Nperm=100)   # 6 min
  # mantel.lm(regr4CV, Nperm=100)   # 6 min
  # mantel.lm(regr5CV, Nperm=100)   # 6 min
  # mantel.lm(regr6CV, Nperm=100)   # 6 min
  # mantel.lm(lm(distancesAD ~ RMED2D + RMED1D,data=distancesCV), Nperm=100)
  # prt.lm(lm(lcv ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=potCV))



  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Distance matrices and Mantel test approach: regression statistics\n
  > R2.lm(regr1CV)
          R2      adjR2
  0.01794498 0.01794091
  > R2.lm(regr2CV)
          R2      adjR2
  0.02638443 0.02637638
  > R2.lm(regr3CV)
          R2      adjR2
  0.03162555 0.03161354
  > R2.lm(regr4CV)
          R2      adjR2
  0.03612941 0.03611347
  > R2.lm(regr5CV)
          R2      adjR2
  0.03838762 0.03836774
  > R2.lm(regr6CV)
          R2      adjR2
  0.03957432 0.03955050

  > RMSEP(regr1CV$model[,1],regr1CV$fitted.values)*100
  [1] 104.1658
  > RMSEP(regr2CV$model[,1],regr2CV$fitted.values)*100
  [1] 103.6993
  > RMSEP(regr3CV$model[,1],regr3CV$fitted.values)*100
  [1] 103.4915
  > RMSEP(regr4CV$model[,1],regr4CV$fitted.values)*100
  [1] 103.4556
  > RMSEP(regr5CV$model[,1],regr5CV$fitted.values)*100
  [1] 103.0755
  > RMSEP(regr6CV$model[,1],regr6CV$fitted.values)*100
  [1] 102.9498
  \n")

  # R2.lm(regr1CV)
  # R2.lm(regr2CV)
  # R2.lm(regr3CV)
  # R2.lm(regr4CV)
  # R2.lm(regr5CV)
  # R2.lm(regr6CV)

  # RMSEP(regr1CV$model[,1],regr1CV$fitted.values)*100
  # RMSEP(regr2CV$model[,1],regr2CV$fitted.values)*100
  # RMSEP(regr3CV$model[,1],regr3CV$fitted.values)*100
  # RMSEP(regr4CV$model[,1],regr4CV$fitted.values)*100
  # RMSEP(regr5CV$model[,1],regr5CV$fitted.values)*100
  # RMSEP(regr6CV$model[,1],regr6CV$fitted.values)*100


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Formation of groups and homogeneity.
  Classification variables:
   0) lnAREA, lnSAAR, BFIHOST (same weights + lnAREA/sqrt(2))\n
  > term1 <- log(cd696[,\"dtm_area\"])/(sd(log(cd696[,\"dtm_area\"]))*sqrt(2))
  > term2 <- log(cd696[,\"saar\"])/sd(log(cd696[,\"saar\"]))
  > term3 <- cd696[,\"bfihost\"]/sd(cd696[,\"bfihost\"])

  > quantile(term1)
          0%        25%        50%        75%       100%
  0.03555657 2.18981948 2.67119010 3.09257122 4.64171915
  > quantile(term2)
        0%      25%      50%      75%     100%
  16.36074 17.35620 18.02938 18.77770 21.15735
  > quantile(term3)
        0%      25%      50%      75%     100%
  1.528373 2.693169 3.137541 3.724989 6.517468
  > roi.cd <- data.frame(cbind(term1,term2,term3))
  > row.names(roi.cd) <- cd696[,\"number\"]

  > roi00.50year <- new.env()
  > for(i in 1:696) {
  + print(paste(i,\"/ 696\"))
  + assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  +         row.names(roi.cd),am696[,\"am\"],am696[,\"number\"],test=\"HW and AD\",station.year=250,Nsim=100), env=roi00.50year)
  + }
  > roi00.50year <- as.list(roi00.50year)

  > estrai.region <- function (x) {x$region}
  > estrai.test <- function (x) {x$test}
  > regioni.50year.0 <- sapply(roi00.50year, estrai.region)
  > test.50year.0 <- sapply(roi00.50year, estrai.test)
  > mL.50year.0 <- mean(sapply(regioni.50year.0,length))
  > mH1.50year.0 <- mean(test.50year.0[\"H1\",])
  > mH2.50year.0 <- mean(test.50year.0[\"H2\",])
  > mAD.50year.0 <- median(test.50year.0[\"P\",])
  > gH1gr2.50year.0 <- sum(test.50year.0[\"H1\",]>2)/696
  > gH1gr4.50year.0 <- sum(test.50year.0[\"H1\",]>4)/696
  > gH2gr2.50year.0 <- sum(test.50year.0[\"H2\",]>2)/696
  > gH2gr4.50year.0 <- sum(test.50year.0[\"H2\",]>4)/696
  > gADgr99.50year.0 <- sum(test.50year.0[\"P\",]>.99)/696
  > gADgr95.50year.0 <- sum(test.50year.0[\"P\",]>.95)/696

  > table.0 <- signif(c(mL.50year.0,mH1.50year.0,mH2.50year.0,mAD.50year.0,
  +            gH1gr2.50year.0*100,gH1gr4.50year.0*100,gH2gr2.50year.0*100,gH2gr4.50year.0*100,
  +            gADgr95.50year.0*100,gADgr99.50year.0*100),3)
  > names(table.0) <- c(\"Avg. n sites\",\"m(H1)\",\"m(H2)\",\"med(p(AD))\",\"% H1>2\",\"% H1>4\",
  +                     \"% H2>2\",\"% H2>4\",\"% p(AD)>0.95\",\"% p(AD)>0.99\")
  > print(table.0)
  \n")

  # term1 <- log(cd696[,"dtm_area"])/(sd(log(cd696[,"dtm_area"]))*sqrt(2))
  # term2 <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
  # term3 <- cd696[,"bfihost"]/sd(cd696[,"bfihost"])

  # roi.cd <- data.frame(cbind(term1,term2,term3))
  # row.names(roi.cd) <- cd696[,"number"]

  # roi00.50year <- new.env()
  # for(i in 1:696) {
  # print(paste(i,"/ 696"))
  # assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  #         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi00.50year)
  # }
  # roi00.50year <- as.list(roi00.50year)

  # estrai.region <- function (x) {x$region}
  # estrai.test <- function (x) {x$test}
  # regioni.50year.0 <- sapply(roi00.50year, estrai.region)
  # test.50year.0 <- sapply(roi00.50year, estrai.test)
  # mL.50year.0 <- mean(sapply(regioni.50year.0,length))
  # mH1.50year.0 <- mean(test.50year.0["H1",])
  # mH2.50year.0 <- mean(test.50year.0["H2",])
  # mAD.50year.0 <- median(test.50year.0["P",])
  # gH1gr2.50year.0 <- sum(test.50year.0["H1",]>2)/696
  # gH1gr4.50year.0 <- sum(test.50year.0["H1",]>4)/696
  # gH2gr2.50year.0 <- sum(test.50year.0["H2",]>2)/696
  # gH2gr4.50year.0 <- sum(test.50year.0["H2",]>4)/696
  # gADgr99.50year.0 <- sum(test.50year.0["P",]>.99)/696
  # gADgr95.50year.0 <- sum(test.50year.0["P",]>.95)/696

  # table.0 <- signif(c(mL.50year.0,mH1.50year.0,mH2.50year.0,mAD.50year.0,
  #            gH1gr2.50year.0*100,gH1gr4.50year.0*100,gH2gr2.50year.0*100,gH2gr4.50year.0*100,
  #            gADgr95.50year.0*100,gADgr99.50year.0*100),3)
  # names(table.0) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
  #                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  # print(table.0)

  table.0 <- c(11.20,2.87,1.53,0.99,62.50,25.10,33.80,5.89,68.20,46.00)
  names(table.0) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  print(table.0)


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Formation of groups and homogeneity. 
  Classification variables:
   1) DPLBAR, lnSAAR, lnRMED2D (with the same weight)\n
  > term1 <- cd696[,\"dplbar\"]/sd(cd696[,\"dplbar\"])
  > term2 <- log(cd696[,\"saar\"])/sd(log(cd696[,\"saar\"]))
  > term3 <- log(cd696[,\"rmed_2d\"])/sd(log(cd696[,\"rmed_2d\"]))

  > quantile(term1)
          0%        25%        50%        75%       100%
  0.06404583 0.56222690 0.90338333 1.45325051 7.57763327
  > quantile(term2)
        0%      25%      50%      75%     100%
  16.36074 17.35620 18.02938 18.77770 21.15735
  > quantile(term3)
        0%      25%      50%      75%     100%
  13.79527 14.99042 15.55170 16.24327 19.08796

  > roi.cd <- data.frame(cbind(term1,term2,term3))
  > row.names(roi.cd) <- cd696[,\"number\"]

  > roi.st.year(roi.cd[\"54088\",],roi.cd,row.names(roi.cd),am696[,\"am\"],
  +             am696[,\"number\"],test=\"HW and AD\",station.year=250,Nsim=500)

  > roi01.50year <- new.env()
  > for(i in 1:696) {
  + print(paste(i,\"/ 696\"))
  + assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  +     row.names(roi.cd),am696[,\"am\"],am696[,\"number\"],test=\"HW and AD\",station.year=250,Nsim=100), env=roi01.50year)
  + }
  > roi01.50year <- as.list(roi01.50year)

  > estrai.region <- function (x) {x$region}
  > estrai.test <- function (x) {x$test}
  > regioni.50year.1 <- sapply(roi01.50year, estrai.region)
  > test.50year.1 <- sapply(roi01.50year, estrai.test)
  > mL.50year.1 <- mean(sapply(regioni.50year.1,length))
  > mH1.50year.1 <- mean(test.50year.1[\"H1\",])
  > mH2.50year.1 <- mean(test.50year.1[\"H2\",])
  > mAD.50year.1 <- median(test.50year.1[\"P\",])
  > gH1gr2.50year.1 <- sum(test.50year.1[\"H1\",]>2)/696
  > gH1gr4.50year.1 <- sum(test.50year.1[\"H1\",]>4)/696
  > gH2gr2.50year.1 <- sum(test.50year.1[\"H2\",]>2)/696
  > gH2gr4.50year.1 <- sum(test.50year.1[\"H2\",]>4)/696
  > gADgr99.50year.1 <- sum(test.50year.1[\"P\",]>.99)/696
  > gADgr95.50year.1 <- sum(test.50year.1[\"P\",]>.95)/696

  > table.1 <- signif(c(mL.50year.1,mH1.50year.1,mH2.50year.1,mAD.50year.1,
  +            gH1gr2.50year.1*100,gH1gr4.50year.1*100,gH2gr2.50year.1*100,gH2gr4.50year.1*100,
  +            gADgr95.50year.1*100,gADgr99.50year.1*100),3)
  > names(table.1) <- c(\"Avg. n sites\",\"m(H1)\",\"m(H2)\",\"med(p(AD))\",\"% H1>2\",\"% H1>4\",
  +                    \"% H2>2\",\"% H2>4\",\"% p(AD)>0.95\",\"% p(AD)>0.99\")
  > print(table.1)
  \n")

  # term1 <- cd696[,"dplbar"]/sd(cd696[,"dplbar"])
  # term2 <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
  # term3 <- log(cd696[,"rmed_2d"])/sd(log(cd696[,"rmed_2d"]))

  # roi.cd <- data.frame(cbind(term1,term2,term3))
  # row.names(roi.cd) <- cd696[,"number"]

  # roi.st.year(roi.cd["54088",],roi.cd,row.names(roi.cd),am696[,"am"],
  #             am696[,"number"],test="HW and AD",station.year=250,Nsim=500)

  # roi01.50year <- new.env()
  # for(i in 1:696) {
  # print(paste(i,"/ 696"))
  # assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  #         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi01.50year)
  # }
  # roi01.50year <- as.list(roi01.50year)

  # estrai.region <- function (x) {x$region}
  # estrai.test <- function (x) {x$test}
  # regioni.50year.1 <- sapply(roi01.50year, estrai.region)
  # test.50year.1 <- sapply(roi01.50year, estrai.test)
  # mL.50year.1 <- mean(sapply(regioni.50year.1,length))
  # mH1.50year.1 <- mean(test.50year.1["H1",])
  # mH2.50year.1 <- mean(test.50year.1["H2",])
  # mAD.50year.1 <- median(test.50year.1["P",])
  # gH1gr2.50year.1 <- sum(test.50year.1["H1",]>2)/696
  # gH1gr4.50year.1 <- sum(test.50year.1["H1",]>4)/696
  # gH2gr2.50year.1 <- sum(test.50year.1["H2",]>2)/696
  # gH2gr4.50year.1 <- sum(test.50year.1["H2",]>4)/696
  # gADgr99.50year.1 <- sum(test.50year.1["P",]>.99)/696
  # gADgr95.50year.1 <- sum(test.50year.1["P",]>.95)/696

  # table.1 <- signif(c(mL.50year.1,mH1.50year.1,mH2.50year.1,mAD.50year.1,
  #            gH1gr2.50year.1*100,gH1gr4.50year.1*100,gH2gr2.50year.1*100,gH2gr4.50year.1*100,
  #            gADgr95.50year.1*100,gADgr99.50year.1*100),3)
  # names(table.1) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
  #                    "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  # print(table.1)
  table.1 <- c(11.30,2.94,1.67,0.99,60.90,27.00,38.20,8.62,65.80,44.70)
  names(table.1) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  print(table.1)





  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Formation of groups and homogeneity.
  Classification variables:
   2) DPLBAR, lnSAAR, lnRMED2D (weighted with regression coefficients)\n
  > coeff <- lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV)$coefficients[-1]
  > print(coeff)

  > term1 <- cd696[,\"dplbar\"]*coeff[1]
  > term2 <- log(cd696[,\"saar\"])*coeff[2]
  > term3 <- log(cd696[,\"rmed_2d\"])*coeff[3]

  > quantile(term1)
            0%          25%          50%          75%         100%
  -0.120140237 -0.023040685 -0.014322768 -0.008913875 -0.001015420
  > quantile(term2)
         0%       25%       50%       75%      100%
  -1.755994 -1.558490 -1.496383 -1.440511 -1.357891
  > quantile(term3)
         0%       25%       50%       75%      100%
  0.6600229 0.7172041 0.7440577 0.7771454 0.9132472

  > roi.cd <- data.frame(cbind(term1,term2,term3))
  > row.names(roi.cd) <- cd696[,\"number\"]

  > roi02.50year <- new.env()
  > for(i in 1:696) {
  + print(paste(i,\"/ 696\"))
  + assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  +     row.names(roi.cd),am696[,\"am\"],am696[,\"number\"],test=\"HW and AD\",station.year=250,Nsim=100), env=roi02.50year)
  + }
  > roi02.50year <- as.list(roi02.50year)

  > estrai.region <- function (x) {x$region}
  > estrai.test <- function (x) {x$test}
  > regioni.50year.2 <- sapply(roi02.50year, estrai.region)
  > test.50year.2 <- sapply(roi02.50year, estrai.test)
  > mL.50year.2 <- mean(sapply(regioni.50year.2,length))
  > mH1.50year.2 <- mean(test.50year.2[\"H1\",])
  > mH2.50year.2 <- mean(test.50year.2[\"H2\",])
  > mAD.50year.2 <- median(test.50year.2[\"P\",])
  > gH1gr2.50year.2 <- sum(test.50year.2[\"H1\",]>2)/696
  > gH1gr4.50year.2 <- sum(test.50year.2[\"H1\",]>4)/696
  > gH2gr2.50year.2 <- sum(test.50year.2[\"H2\",]>2)/696
  > gH2gr4.50year.2 <- sum(test.50year.2[\"H2\",]>4)/696
  > gADgr99.50year.2 <- sum(test.50year.2[\"P\",]>.99)/696
  > gADgr95.50year.2 <- sum(test.50year.2[\"P\",]>.95)/696

  > table.2 <- signif(c(mL.50year.2,mH1.50year.2,mH2.50year.2,mAD.50year.2,
  +            gH1gr2.50year.2*100,gH1gr4.50year.2*100,gH2gr2.50year.2*100,gH2gr4.50year.2*100,
  +            gADgr95.50year.2*100,gADgr99.50year.2*100),3)
  > names(table.2) <- c(\"Avg. n sites\",\"m(H1)\",\"m(H2)\",\"med(p(AD))\",\"% H1>2\",\"% H1>4\",
  +                     \"% H2>2\",\"% H2>4\",\"% p(AD)>0.95\",\"% p(AD)>0.99\")
  > print(table.2)
  \n")

  coeff <- lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV)$coefficients[-1]
  print(coeff)

  # term1 <- cd696[,"dplbar"]*coeff[1]
  # term2 <- log(cd696[,"saar"])*coeff[2]
  # term3 <- log(cd696[,"rmed_2d"])*coeff[3]

  # roi.cd <- data.frame(cbind(term1,term2,term3))
  # row.names(roi.cd) <- cd696[,"number"]

  # roi02.50year <- new.env()
  # for(i in 1:696) {
  # print(paste(i,"/ 696"))
  # assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  #         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi02.50year)
  # }
  # roi02.50year <- as.list(roi02.50year)

  # estrai.region <- function (x) {x$region}
  # estrai.test <- function (x) {x$test}
  # regioni.50year.2 <- sapply(roi02.50year, estrai.region)
  # test.50year.2 <- sapply(roi02.50year, estrai.test)
  # mL.50year.2 <- mean(sapply(regioni.50year.2,length))
  # mH1.50year.2 <- mean(test.50year.2["H1",])
  # mH2.50year.2 <- mean(test.50year.2["H2",])
  # mAD.50year.2 <- median(test.50year.2["P",])
  # gH1gr2.50year.2 <- sum(test.50year.2["H1",]>2)/696
  # gH1gr4.50year.2 <- sum(test.50year.2["H1",]>4)/696
  # gH2gr2.50year.2 <- sum(test.50year.2["H2",]>2)/696
  # gH2gr4.50year.2 <- sum(test.50year.2["H2",]>4)/696
  # gADgr99.50year.2 <- sum(test.50year.2["P",]>.99)/696
  # gADgr95.50year.2 <- sum(test.50year.2["P",]>.95)/696

  # table.2 <- signif(c(mL.50year.2,mH1.50year.2,mH2.50year.2,mAD.50year.2,
  #            gH1gr2.50year.2*100,gH1gr4.50year.2*100,gH2gr2.50year.2*100,gH2gr4.50year.2*100,
  #            gADgr95.50year.2*100,gADgr99.50year.2*100),3)
  # names(table.2) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
  #                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  # print(table.2)

  table.2 <- c(11.20,2.95,1.64,0.99,57.20,28.60,37.90,9.34,66.40,45.10) # weigths doesnt change the result
  names(table.2) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  print(table.2)


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Formation of groups and homogeneity.
  Classification variables:
   3) FARL, RMED2D, lnSAAR (same weight)\n
  > term1 <- cd696[,\"farl\"]/sd(cd696[,\"farl\"])
  > term2 <- cd696[,\"rmed_2d\"]/sd(cd696[,\"rmed_2d\"])
  > term3 <- log(cd696[,\"saar\"])/sd(log(cd696[,\"saar\"]))

  > quantile(term1)
        0%      25%      50%      75%     100%
  12.79836 20.77582 21.24136 21.45746 21.50261
  > quantile(term2)
        0%      25%      50%      75%     100%
  2.165229 2.925077 3.368881 4.009372 8.203663
  > quantile(term3)
        0%      25%      50%      75%     100%
  16.36074 17.35620 18.02938 18.77770 21.15735

  > roi.cd <- data.frame(cbind(term1,term2,term3))
  > row.names(roi.cd) <- cd696[,\"number\"]

  > roi03.50year <- new.env()
  > for(i in 1:696) {
  + print(paste(i,\"/ 696\"))
  + assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  +         row.names(roi.cd),am696[,\"am\"],am696[,\"number\"],test=\"HW and AD\",station.year=250,Nsim=100), env=roi03.50year)
  + }
  > roi03.50year <- as.list(roi03.50year)

  > estrai.region <- function (x) {x$region}
  > estrai.test <- function (x) {x$test}
  > regioni.50year.3 <- sapply(roi03.50year, estrai.region)
  > test.50year.3 <- sapply(roi03.50year, estrai.test)
  > mL.50year.3 <- mean(sapply(regioni.50year.3,length))
  > mH1.50year.3 <- mean(test.50year.3[\"H1\",])
  > mH2.50year.3 <- mean(test.50year.3[\"H2\",])
  > mAD.50year.3 <- median(test.50year.3[\"P\",])
  > gH1gr2.50year.3 <- sum(test.50year.3[\"H1\",]>2)/696
  > gH1gr4.50year.3 <- sum(test.50year.3[\"H1\",]>4)/696
  > gH2gr2.50year.3 <- sum(test.50year.3[\"H2\",]>2)/696
  > gH2gr4.50year.3 <- sum(test.50year.3[\"H2\",]>4)/696
  > gADgr99.50year.3 <- sum(test.50year.3[\"P\",]>.99)/696
  > gADgr95.50year.3 <- sum(test.50year.3[\"P\",]>.95)/696

  > table.3 <- signif(c(mL.50year.3,mH1.50year.3,mH2.50year.3,mAD.50year.3,
  +            gH1gr2.50year.3*100,gH1gr4.50year.3*100,gH2gr2.50year.3*100,gH2gr4.50year.3*100,
  +            gADgr95.50year.3*100,gADgr99.50year.3*100),3)
  > names(table.3) <- c(\"Avg. n sites\",\"m(H1)\",\"m(H2)\",\"med(p(AD))\",\"% H1>2\",\"% H1>4\",
  +                     \"% H2>2\",\"% H2>4\",\"% p(AD)>0.95\",\"% p(AD)>0.99\")
  > print(table.3)
  \n")

  # term1 <- cd696[,"farl"]/sd(cd696[,"farl"])
  # term2 <- cd696[,"rmed_2d"]/sd(cd696[,"rmed_2d"])
  # term3 <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))

  # roi.cd <- data.frame(cbind(term1,term2,term3))
  # row.names(roi.cd) <- cd696[,"number"]

  # roi03.50year <- new.env()
  # for(i in 1:696) {
  # print(paste(i,"/ 696"))
  # assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  #         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi03.50year)
  # }
  # roi03.50year <- as.list(roi03.50year)

  # estrai.region <- function (x) {x$region}
  # estrai.test <- function (x) {x$test}
  # regioni.50year.3 <- sapply(roi03.50year, estrai.region)
  # test.50year.3 <- sapply(roi03.50year, estrai.test)
  # mL.50year.3 <- mean(sapply(regioni.50year.3,length))
  # mH1.50year.3 <- mean(test.50year.3["H1",])
  # mH2.50year.3 <- mean(test.50year.3["H2",])
  # mAD.50year.3 <- median(test.50year.3["P",])
  # gH1gr2.50year.3 <- sum(test.50year.3["H1",]>2)/696
  # gH1gr4.50year.3 <- sum(test.50year.3["H1",]>4)/696
  # gH2gr2.50year.3 <- sum(test.50year.3["H2",]>2)/696
  # gH2gr4.50year.3 <- sum(test.50year.3["H2",]>4)/696
  # gADgr99.50year.3 <- sum(test.50year.3["P",]>.99)/696
  # gADgr95.50year.3 <- sum(test.50year.3["P",]>.95)/696

  # table.3 <- signif(c(mL.50year.3,mH1.50year.3,mH2.50year.3,mAD.50year.3,
  #            gH1gr2.50year.3*100,gH1gr4.50year.3*100,gH2gr2.50year.3*100,gH2gr4.50year.3*100,
  #            gADgr95.50year.3*100,gADgr99.50year.3*100),3)
  # names(table.3) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
  #                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  # print(table.3)

  table.3 <- c(11.00,3.37,1.60,0.99,67.20,35.10,34.60,8.05,67.40,46.30)
  names(table.3) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  print(table.3)



  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Formation of groups and homogeneity.
  Classification variables:
   0) lnAREA, lnSAAR, BFIHOST (same weight + lnAREA/sqrt(2))
   1) DPLBAR, lnSAAR, lnRMED2D (with the same weight)
   2) DPLBAR, lnSAAR, lnRMED2D (weighted with regression coefficients)
   3) FARL, RMED2D, lnSAAR (with the same weight)\n
  > table.0123 <- rbind(table.0,table.1,table.2,table.3)
  > print(table.0123)
  \n")

  table.0123 <- rbind(table.0,table.1,table.2,table.3)
  print(table.0123)


  readline("\n
  # ------------------------------------------------------------------------------------- #\n
  Press Return to continue
  \n")

  cat("\n
  Formation of groups and homogeneity (other cases).
  Classification variables:
  > poolingroupsum <- function(roi.cd) {
  +  roi.50year <- new.env()
  +  for(i in 1:696) {
  +  print(paste(i,\"/ 696\"))
  +  assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
  +          row.names(roi.cd),am696[,\"am\"],am696[,\"number\"],test=\"HW and AD\",station.year=250,Nsim=100), env=roi.50year)
  +  }
  +  roi.50year <- as.list(roi.50year)
  +
  +  estrai.region <- function (x) {x$region}
  +  estrai.test <- function (x) {x$test}
  +  regioni.50year <- sapply(roi.50year, estrai.region)
  +  test.50year <- sapply(roi.50year, estrai.test)
  +  mL.50year <- mean(sapply(regioni.50year,length))
  +  mH1.50year <- mean(test.50year[\"H1\",])
  +  mH2.50year <- mean(test.50year[\"H2\",])
  +  mAD.50year <- median(test.50year[\"P\",])
  +  gH1gr2.50year <- sum(test.50year[\"H1\",]>2)/696
  +  gH1gr4.50year <- sum(test.50year[\"H1\",]>4)/696
  +  gH2gr2.50year <- sum(test.50year[\"H2\",]>2)/696
  +  gH2gr4.50year <- sum(test.50year[\"H2\",]>4)/696
  +  gADgr99.50year <- sum(test.50year[\"P\",]>.99)/696
  +  gADgr95.50year <- sum(test.50year[\"P\",]>.95)/696
  +
  +  table <- signif(c(mL.50year,mH1.50year,mH2.50year,mAD.50year,
  +             gH1gr2.50year*100,gH1gr4.50year*100,gH2gr2.50year*100,gH2gr4.50year*100,
  +             gADgr95.50year*100,gADgr99.50year*100),3)
  +  names(table) <- c(\"Avg. n sites\",\"m(H1)\",\"m(H2)\",\"med(p(AD))\",\"% H1>2\",\"% H1>4\",
  +                      \"% H2>2\",\"% H2>4\",\"% p(AD)>0.95\",\"% p(AD)>0.99\")
  +  table
  + }

  > termFARL <- cd696[,\"farl\"]/sd(cd696[,\"farl\"])
  > termSMDBAR <- cd696[,\"smdbar\"]/sd(cd696[,\"smdbar\"])
  > termDPLBAR <- cd696[,\"dplbar\"]/sd(cd696[,\"dplbar\"])
  > termlnNGRY <- log(cd696[,\"ihdtm_ngr_y\"])/sd(log(cd696[,\"ihdtm_ngr_y\"]))
  > termlnSAAR <- log(cd696[,\"saar\"])/sd(log(cd696[,\"saar\"]))
  > termlnRMED2D <- log(cd696[,\"rmed_2d\"])/sd(log(cd696[,\"rmed_2d\"]))
  > termRMED1D <- cd696[,\"rmed_1d\"]/sd(cd696[,\"rmed_1d\"])
  > termALTBAR <- cd696[,\"altbar\"]/sd(cd696[,\"altbar\"])
  > termlnAREA <- log(cd696[,\"dtm_area\"])/sd(log(cd696[,\"dtm_area\"]))
  > termlnALTBAR <- log(cd696[,\"altbar\"])/sd(log(cd696[,\"altbar\"]))

  > table.4 <- poolingroupsum(data.frame(cbind(termFARL,termSMDBAR,termDPLBAR,termlnNGRY,termlnSAAR,termlnRMED2D),
  >            row.names=cd696[,\"number\"]))
  > table.5 <- poolingroupsum(data.frame(cbind(termRMED1D,termALTBAR,termlnAREA,termlnSAAR,termlnRMED2D,termlnALTBAR),
  >            row.names=cd696[,\"number\"]))

  \nClassification variables:
   0) lnAREA, lnSAAR, BFIHOST (same weight + lnAREA/sqrt(2))
   1) DPLBAR, lnSAAR, lnRMED2D (with the same weight)
   2) DPLBAR, lnSAAR, lnRMED2D (weighted with regression coefficients)
   3) FARL, RMED2D, lnSAAR (with the same weight)
   4) FARL, SMDBAR, DPLBAR, lnNGRY, lnSAAR, lnRMED2D (with the same weight)
   5) RMED1D, ALTBAR, lnAREA, lnSAAR, lnRMED2D, lnALTBAR (with the same weight)\n
  > print(rbind(table.4,table.5))
  > print(rbind(table.0123,table.4,table.5))
  \n")

  poolingroupsum <- function(roi.cd) {
   roi.50year <- new.env()
   for(i in 1:696) {
   print(paste(i,"/ 696"))
   assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
           row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi.50year)
   }
   roi.50year <- as.list(roi.50year)

   estrai.region <- function (x) {x$region}
   estrai.test <- function (x) {x$test}
   regioni.50year <- sapply(roi.50year, estrai.region)
   test.50year <- sapply(roi.50year, estrai.test)
   mL.50year <- mean(sapply(regioni.50year,length))
   mH1.50year <- mean(test.50year["H1",])
   mH2.50year <- mean(test.50year["H2",])
   mAD.50year <- median(test.50year["P",])
   gH1gr2.50year <- sum(test.50year["H1",]>2)/696
   gH1gr4.50year <- sum(test.50year["H1",]>4)/696
   gH2gr2.50year <- sum(test.50year["H2",]>2)/696
   gH2gr4.50year <- sum(test.50year["H2",]>4)/696
   gADgr99.50year <- sum(test.50year["P",]>.99)/696
   gADgr95.50year <- sum(test.50year["P",]>.95)/696
 
   table <- signif(c(mL.50year,mH1.50year,mH2.50year,mAD.50year,
              gH1gr2.50year*100,gH1gr4.50year*100,gH2gr2.50year*100,gH2gr4.50year*100,
              gADgr95.50year*100,gADgr99.50year*100),3)
   names(table) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                       "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
   table
  }

  termFARL <- cd696[,"farl"]/sd(cd696[,"farl"])
  termSMDBAR <- cd696[,"smdbar"]/sd(cd696[,"smdbar"])
  termDPLBAR <- cd696[,"dplbar"]/sd(cd696[,"dplbar"])
  termlnNGRY <- log(cd696[,"ihdtm_ngr_y"])/sd(log(cd696[,"ihdtm_ngr_y"]))
  termlnSAAR <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
  termlnRMED2D <- log(cd696[,"rmed_2d"])/sd(log(cd696[,"rmed_2d"]))
  termRMED1D <- cd696[,"rmed_1d"]/sd(cd696[,"rmed_1d"])
  termALTBAR <- cd696[,"altbar"]/sd(cd696[,"altbar"])
  termlnAREA <- log(cd696[,"dtm_area"])/sd(log(cd696[,"dtm_area"]))
  termlnALTBAR <- log(cd696[,"altbar"])/sd(log(cd696[,"altbar"]))

  # table.4 <- poolingroupsum(data.frame(cbind(termFARL,termSMDBAR,termDPLBAR,termlnNGRY,termlnSAAR,termlnRMED2D),
  #            row.names=cd696[,"number"]))
  # table.5 <- poolingroupsum(data.frame(cbind(termRMED1D,termALTBAR,termlnAREA,termlnSAAR,termlnRMED2D,termlnALTBAR),
  #            row.names=cd696[,"number"]))

  table.4 <- c(11.10,2.68,1.56,0.98,61.20,24.10,35.20,6.03,61.90,40.40)
  names(table.4) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
  table.5 <- c(11.1,2.91,1.77,0.99,62.5,24.9,39.9,9.05,65.9,44.8)
  names(table.4) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")

  print(rbind(table.4,table.5))
  print(rbind(table.0123,table.4,table.5))

  cat("\n\n
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
                                          THE END
  # ------------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------------- #
  \n")
}
