consistencyplot <- function (x,cod,cex.axis=.8,mark.code=TRUE,...) {

 # INPUT:
 # x = years
 # cod = code

 cod <- as.factor(cod)
 n <- nlevels(cod)
 codici <- levels(cod)
 ni <- tapply(x,cod,length)

 dimensioni=cex.axis
  
 time <- c(min(x),max(x))
 code <- c(1,n)
 plot(time,code,type="n",axes=FALSE,...)
 axis(1)
 if (mark.code==TRUE) axis(2,at=c(n:1),labels=codici,las=1,cex.axis=dimensioni)
 box(); grid()
 for (i in n:1) {
  xi <- x[cod==codici[i]]
  points(xi,rep(n-i+1,ni[i]),pch=45)
 }
}
