# 2005-09-23, Alberto Viglione
#
# A.9) Pearson type III distribution

# gamm=0 => normal, gamm=2 => exponential, gamm=-2 => reverse exponential distribution

f.gamma <- function (x,mu,sigma,gamm) {

  if (gamm > 0) {
    alfa <- 4/(gamm^2)
    beta <- 0.5 * sigma * abs(gamm)
    xi <- mu - 2*sigma/gamm
    #f <- dgamma((x - xi)/beta, alfa)
    f <- ((x - xi)^(alfa - 1) * exp(-(x - xi)/beta))/(beta^alfa * gamma(alfa))
  }
  else if (gamm < 0) {
    alfa <- 4/(gamm^2)
    beta <- 0.5 * sigma * abs(gamm)
    xi <- mu - 2*sigma/gamm
    #f <- dgamma((xi - x)/beta, alfa) 
    f <- ((xi - x)^(alfa - 1) * exp(-(xi - x)/beta))/(beta^alfa * gamma(alfa))
  }
  else if (gamm==0) {
    f <- sigma^(-1) * dnorm((x - mu)/sigma)	# Normal
  }

  return(f)
}

F.gamma <- function (x,mu,sigma,gamm) {

  if (gamm > 0) {
    alfa <- 4/(gamm^2)
    beta <- 0.5 * sigma * abs(gamm)
    xi <- mu - 2*sigma/gamm
    F <- pgamma((x - xi)/beta, alfa)
  }
  else if (gamm < 0) {
    alfa <- 4/(gamm^2)
    beta <- 0.5 * sigma * abs(gamm)
    xi <- mu - 2*sigma/gamm
    F <- 1 - pgamma((xi - x)/beta, alfa) 
  }
  else if (gamm==0) {
    F <- pnorm((x - mu)/sigma)	# Normal
  }

  return(F)
}

invF.gamma <- function (F,mu,sigma,gamm) {

  if ((F < 0) || (F > 1)) {
    stop("F must be between 0 and 1")
  } 

  if (gamm > 0) {
    alfa <- 4/(gamm^2)
    beta <- 0.5 * sigma * abs(gamm)
    xi <- mu - 2*sigma/gamm
    x.st <- qgamma(F, alfa)
    x <- x.st*beta + xi
  }
  else if (gamm < 0) {
    alfa <- 4/(gamm^2)
    beta <- 0.5 * sigma * abs(gamm)
    xi <- mu - 2*sigma/gamm
    x.st <- qgamma(F, alfa)
    x <- xi - x.st*beta
  }
  else if (gamm==0) {
    x.st <- qnorm(F)	# Normal
    x <- mu + sigma*x.st
  }

  return(x)
}

Lmom.gamma <- function(mu,sigma,gamm) {

  A0 = 0.32573501
  A1 = 0.16869150
  A2 = 0.078327243
  A3 = -0.0029120539
  B1 = 0.46697102
  B2 = 0.24255406
  C0 = 0.12260172
  C1 = 0.053730130
  C2 = 0.043384378
  C3 = 0.011101277
  D1 = 0.18324466
  D2 = 0.20166036
  E1 = 2.3807576
  E2 = 1.5931792
  E3 = 0.11618371
  F1 = 5.1533299
  F2 = 7.1425260
  F3 = 1.9745056
  G1 = 2.1235833
  G2 = 4.1670213
  G3 = 3.1925299
  H1 = 9.0551443
  H2 = 0.26649995
  H3 = 0.26193668
  
  quanti <- length(mu)
  lambda1 <- rep(NA,quanti)
  lambda2 <- rep(NA,quanti)
  tau3 <- rep(NA,quanti)
  tau4 <- rep(NA,quanti)
  for (i in 1:quanti) {
    if (gamm[i] > 0) {
      alfa <- 4/(gamm[i]^2)
      beta <- 0.5 * sigma[i] * abs(gamm[i])
      xi <- mu[i] - 2*sigma[i]/gamm[i]
      lambda1[i] <- xi + alfa*beta
      lambda2[i] <- pi^(-0.5) *beta * gamma(alfa + 0.5)/gamma(alfa)
      # tau3 <- 6 * pbeta(1/3,alfa,2*alfa) - 3
      if (alfa >= 1) {
        tau3[i] <- alfa^(-0.5) * (A0 + A1*alfa^(-1) + A2*alfa^(-2) + A3*alfa^(-3))/(1 + B1*alfa^(-1) + B2*alfa^(-2))
        tau4[i] <- (C0 + C1*alfa^(-1) + C2*alfa^(-2) + C3*alfa^(-3))/(1 + D1*alfa^(-1) + D2*alfa^(-2))
      }
      else if (alfa < 1) {
        tau3[i] <- (1 + E1*alfa + E2*alfa^2 + E3*alfa^3)/(1 + F1*alfa + F2*alfa^2 + F3*alfa^3)
        tau4[i] <- (1 + G1*alfa + G2*alfa^2 + G3*alfa^3)/(1 + H1*alfa + H2*alfa^2 + H3*alfa^3)
      }
    }
    else if (gamm[i] < 0) {
      alfa <- 4/(gamm[i]^2)
      beta <- 0.5 * sigma[i] * abs(gamm[i])
      xi <- mu[i] - 2*sigma[i]/gamm[i]
      lambda1[i] <- -(-xi + alfa*beta)
      lambda2[i] <- pi^(-0.5) *beta * gamma(alfa + 0.5)/gamma(alfa)
      if (alfa >= 1) {
        tau3[i] <- -(alfa^(-0.5) * (A0 + A1*alfa^(-1) + A2*alfa^(-2) + A3*alfa^(-3))/(1 + B1*alfa^(-1) + B2*alfa^(-2)))
        tau4[i] <- (C0 + C1*alfa^(-1) + C2*alfa^(-2) + C3*alfa^(-3))/(1 + D1*alfa^(-1) + D2*alfa^(-2))
      }
      else if (alfa < 1) {
        tau3[i] <- -(1 + E1*alfa + E2*alfa^2 + E3*alfa^3)/(1 + F1*alfa + F2*alfa^2 + F3*alfa^3)
        tau4[i] <- (1 + G1*alfa + G2*alfa^2 + G3*alfa^3)/(1 + H1*alfa + H2*alfa^2 + H3*alfa^3)
      }
    }
    else if (gamm[i]==0) {
      lambda1[i] <- mu[i]
      lambda2[i] <- sigma[i]*pi^(-0.5)
      tau3[i] <- 0
      tau4[i] <- 30*pi^(-1)*atan(sqrt(2)) - 9	# 0.1226
    }
  }
  
  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)
  
  return(output)
}

par.gamma <- function(lambda1,lambda2,tau3) {

  lambda1 <- as.numeric(lambda1)
  lambda2 <- as.numeric(lambda2)
  tau3 <- as.numeric(tau3)
  
  quanti <- length(tau3)
  mu <- rep(NA,quanti)
  sigma <- rep(NA,quanti)
  gamm <- rep(NA,quanti)

  for (i in 1:quanti) {
    if (tau3[i]==0) {
      mu[i] <- lambda1[i]
      sigma[i] <- pi^(0.5) * lambda2[i]
      gamm[i] <- 0
    }
    else {
      if ((abs(tau3[i]) > 0)&&(abs(tau3[i]) < 1/3)) {
        z <- 3*pi*tau3[i]^2
        alfa <- (1 + 0.2906*z)/(z + 0.1882*z^2 + 0.0442*z^3)
      }
      else if ((abs(tau3[i]) >= 1/3)&&(abs(tau3[i]) < 1)) {
        z <- 1 - abs(tau3[i])
        alfa <- (0.36067*z - 0.59567*z^2 + 0.25361*z^3)/(1 - 2.78861*z + 2.56096*z^2 - 0.77045*z^3)
      }

      gamm[i] <- 2*alfa^(-0.5) * sign(tau3[i])
      sigma[i] <- suppressWarnings(lambda2[i]*pi^(0.5) * alfa^(0.5) * gamma(alfa)/gamma(alfa + 0.5))
      mu[i] <- lambda1[i]
    }
  }
  output <- list(mu=mu, sigma=sigma, gamm=gamm)
 
  return(output)
}

rand.gamma <- function(numerosita,mu,sigma,gamm) {

  F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  x <- invF.gamma(F,mu,sigma,gamm)

  return(x)
}

mom2par.gamma <- function(mu,sigma,gamm) {
  
  if(gamm==0) {stop("The distribution is Normal")}
  else {
   alpha <- 4/(gamm^2)
   beta <- 0.5*sigma*abs(gamm) 
   xi <- mu - 2*sigma/gamm
  }
  
  output <- list(alpha=alpha, beta=beta, xi=xi)

  return(output)
}
