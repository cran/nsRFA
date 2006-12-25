gofP3test <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.gamma(Lmom.x["l1"],Lmom.x["l2"],Lmom.x["lca"])
 F <- F.gamma(x,par$mu,par$sigma,par$gamm)
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  x.sim <- rand.gamma(n,par$mu,par$sigma,par$gamm)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.gamma(Lmom.xsim["l1"],Lmom.xsim["l2"],Lmom.xsim["lca"])
  F <- F.gamma(x.sim,par.sim$mu,par.sim$sigma,par.sim$gamm)
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }
 
 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# ----------------------------------------------------------------------------- #

gofNORMtest <- function(x) {

 # Francesco Laio

 fw2 <- function (x) {
  if (x < 1.2) {
   fw2 <- ((exp(-(1/16)/x)*besselK((1/16)/x,1/4)+1.11803*exp(-(25/16)/x)*besselK((25/16)/x,1/4))/(x^0.5))/pi
  }
  else fw2 <- 1
  return(fw2)
 }

 b <- sort(x)
 n <-length(x)
 c0=0.851; b0=0.116; csi0=0.0403; eps1=1.2; eps2=0.2

 Mgaus <- c(mean(b),sum((b-mean(b))^2 / n)^0.5)
 Fgaus <- pmax(pmin(pnorm(b,Mgaus[1],Mgaus[2]),0.999999999),0.00000001)

 # Anderson-Darling

 c1=1.147; b1=0.229; csi1=0.167
 c1corr <- c1*(1+0.5/n); b1corr <- b1*(1-0.2/n); csi1corr <- csi1*(1+0.3/n)

 Agaus <- -n-sum((2*(1:n)-1)*log(Fgaus)+(2*n+1-2*(1:n))*log(1-Fgaus))/n

 if (Agaus <= eps1*csi1corr) {
   z3 <- max((csi0+b0*((eps1-1)*csi1corr/b1corr)^(c1corr/c0))/((eps1-eps2)*csi1corr)*(Agaus-eps2*csi1corr),0.00001)
 }
 else {
   z3 <- csi0+b0*((Agaus-csi1corr)/b1corr)^(c1corr/c0)
 }

 pA <- fw2(z3)
 
 output <- c(Agaus,pA); names(output) <- c("A2","P")
 return(output)
}

