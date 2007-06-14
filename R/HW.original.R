HW.original <- function (data, cod, Nsim=500) {
 # x = data
 # cod = codes
 #
 # PARAMETERS OF ROUTINE:
 # NSITES * INPUT* NUMBER OF SITES IN REGION
 # NAMES * INPUT* CHARACTER*12 ARRAY OF LENGTH NSITES. SITE NAMES.
 # LEN * INPUT* ARRAY OF LENGTH NSITES. RECORD LENGTHS AT EACH SITE.
 # XMOM * INPUT* ARRAY OF DIMENSION (5,NSITES). ARRAY CONTAINING
 #   THE FIRST 5 SAMPLE L-MOMENTS FOR EACH SITE, IN THE
 #   ORDER MEAN, L-CV, L-SKEWNESS, L-KURTOSIS, T-5, I.E
 #   XMOM(I,J) CONTAINS THE I'TH L-MOMENT FOR SITE J.
 #   N.B. XMOM(2,.) CONTAINS L-CV, NOT THE USUAL L-2!
 # A * INPUT* ) PARAMETERS OF
 # B * INPUT* ) PLOTTING POSITION.
 #   NOTE: A AND B SHOULD BE THE SAME AS THE VALUES USED
 #   TO CALCULATE THE MOMENTS IN THE XMOM ARRAY.
 # SEED * INPUT* SEED FOR RANDOM NUMBER GENERATOR. SHOULD BE A WHOLE
 #   NUMBER IN THE RANGE 2D0 TO 2147483647D0.
 # NSIM * INPUT* NUMBER OF SIMULATED WORLDS FOR HETEROGENEITY AND
 #   GOODNESS-OF-FIT TESTS.
 #   NOTE: NSIM=0 WILL FORCE RETURN AT COMPLETION OF
 #   OUTLIER TEST. NSIM=1 WILL SUPPRESS CALCULATION OF
 #   H AND Z STATISTICS, BUT PARAMETER AND QUANTILE
 #   ESTIMATES WILL BE FOUND.
 # NPROB * INPUT* NUMBER OF QUANTILES TO BE CALCULATED
 # PROB * INPUT* ARRAY OF LENGTH NPROB. PROBABILITIES FOR WHICH
 #   QUANTILES ARE TO BE CALCULATED.
 # KPRINT * INPUT* OUTPUT FLAG. SHOULD BE SET TO
 #   0 TO SUPPRESS OUTPUT
 #   1 TO PRINT OUTPUT
 # KOUT * INPUT* CHANNEL TO WHICH OUTPUT IS DIRECTED
 # RMOM *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE REGIONAL
 #   WEIGHTED AVERAGE L-MOMENT RATIOS.
 # D *OUTPUT* ARRAY OF LENGTH NSITES. ON EXIT, CONTAINS THE
 #   DISCORDANCY MEASURE (D STATISTIC) FOR EACH SITE.
 # VOBS *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE REGIONAL
 #   OBSERVED VALUES OF 3 HETEROGENEITY STATISTICS:
 #   (1) WEIGHTED S.D. OF L-CVS;
 #   (2) AVERAGE OF L-CV/L-SKEW DISTANCES;
 #   (3) AVERAGE OF L-SKEW/L-KURTOSIS DISTANCES.
 # VBAR *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE MEAN OF THE
 #   SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
 # VSD *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE S.D. OF THE
 #   SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
 # H *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS HETEROGENEITY
 #   MEASURES (H STATISTICS), I.E. H=(VOBS-VBAR)/VSD.
 # Z *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS GOODNESS-OF-FIT
 #   MEASURES (Z STATISTICS) FOR 5 DISTRIBUTIONS:
 #   (1) GEN. LOGISTIC, (2) GEN. EXTREME VALUE,
 #   (3) GEN. NORMAL, (4) PEARSON TYPE III,
 #   (5) GEN. PARETO.
 # PARA *OUTPUT* ARRAY OF DIMENSION (5,6). ON EXIT, IF NSIM.GE.1,
 #   CONTAINS PARAMETERS OF GROWTH CURVES FITTED BY THE
 #   ABOVE 5 DISTRIBUTIONS, PLUS WAKEBY.
 # RPARA *OUTPUT* PARAMETERS OF REGIONAL KAPPA DISTRIBUTION

 x <- data

 # Input
 NSITES <- as.integer(length(unique(cod)))
 NAMES <- unique(cod)
 LEN <- as.integer(tapply(x, cod, length))
 XMOM <- array(rbind(sapply(split(x,cod),Lmoments)[c(1,3,4,5),], rep(0,NSITES)), dim=c(5,NSITES))
 A <- 0.0  # -0.35
 B <- 0.0
 # FOR UNBIASED ESTIMATES (OF THE LAMBDA'S) SET A=B=ZERO. OTHERWISE,
 # PLOTTING-POSITION ESTIMATORS ARE USED, BASED ON THE PLOTTING POSITION
 # (J+A)/(N+B) FOR THE J'TH SMALLEST OF N OBSERVATIONS. FOR EXAMPLE,
 # A=-0.35D0 AND B=0.0D0 YIELDS THE ESTIMATORS RECOMMENDED BY
 # HOSKING ET AL. (1985, TECHNOMETRICS) FOR THE GEV DISTRIBUTION.
 SEED <- .Random.seed[1]
 NSIM <- as.integer(Nsim)
 PROB <- as.double(c(.01, .02, .05, .1, .2, .5, .9, .95, .99, .999))
 NPROB <- as.integer(length(PROB))
 KPRINT <- as.integer(0)
 KOUT <- as.integer(6)   # ???

 # Output initialization
 RMOM <- double(length=5)   # as.numeric(rep(1, 5))
 D <- double(length=NSITES)   # rep(1, NSITES)
 VOBS <- double(length=3)   # rep(1, 3)
 VBAR <- double(length=3)   # rep(1, 3)
 VSD <- double(length=3)   # rep(1, 3)
 H <- double(length=3)   # rep(1, 3)
 Z <- double(length=5)   # rep(1, 5)
 PARA <- matrix(1, nrow=5, ncol=6)
 RPARA <- double(length=4)

 out <- .Fortran("REGTST", NSITES=NSITES, NAMES=NAMES, LEN=LEN, XMOM=XMOM, A=A, B=B,
                           SEED=SEED, NSIM=NSIM, NPROB=NPROB, PROB=PROB, KPRINT=KPRINT,
                           KOUT=KOUT, RMOM=RMOM, D=D, VOBS=VOBS, VBAR=VBAR, VSD=VSD, H=H, Z=Z,
                           PARA=PARA, RPARA=RPARA, PACKAGE="nsRFA")
 class(out) <- c("HWorig")
 return(out)
}


.First.lib <- function(libname, pkgname) {
 library.dynam("nsRFA", pkgname, libname)
}

# ------------------------------------------------------------------- #

print.HWorig <- function(x, ...) {
 # x = object of classHWorig

 cat("\n\nRESULTS OBTAINED WITH THE FORTRAN ROUTINE OF HOSKING\n\n")

 names <- x$NAMES
 len <- x$LEN
 xmom <- x$XMOM
 d <- x$D
 tabella01 <- data.frame(cbind(len,t(xmom)[,2:4],d), row.names=names)
 names(tabella01) <- c("n", "L-CV", "L-SKEW", "L-KURT", "D(I)")
 print(tabella01)

 cat("\n\nPARAMETERS OF REGIONAL KAPPA DISTRIBUTION:\n")
 print(x$RPARA)


 cat("\n\n***** HETEROGENEITY MEASURES *****\n")
 cat("\nNUMBER OF SIMULATIONS = ", x$NSIM, "\n", sep="")

 cat("\nOBSERVED S.D. OF GROUP L-CV     = ", x$VOBS[1], "\n", sep="")
 cat("SIM. MEAN OF S.D. OF GROUP L-CV = ", x$VBAR[1], "\n", sep="")
 cat("SIM. S.D. OF S.D. OF GROUP L-CV = ", x$VSD[1], "\n", sep="")
 cat("STANDARDIZED TEST VALUE H(1)    = ", x$H[1], "\n", sep="")

 cat("\nOBSERVED AVE. OF L-CV / L-SKEW DISTANCE  = ", x$VOBS[2], "\n", sep="")
 cat("SIM. MEAN OF AVE. L-CV / L-SKEW DISTANCE = ", x$VBAR[2], "\n", sep="")
 cat("SIM. S.D. OF AVE. L-CV / L-SKEW DISTANCE = ", x$VSD[2], "\n", sep="")
 cat("STANDARDIZED TEST VALUE H(2)             = ", x$H[2], "\n", sep="")

 cat("\nOBSERVED AVE. OF L-SKEW/L-KURT DISTANCE  = ", x$VOBS[3], "\n", sep="")
 cat("SIM. MEAN OF AVE. L-SKEW/L-KURT DISTANCE = ", x$VBAR[3], "\n", sep="")
 cat("SIM. S.D. OF AVE. L-SKEW/L-KURT DISTANCE = ", x$VSD[3], "\n", sep="")
 cat("STANDARDIZED TEST VALUE H(3)             = ", x$H[3], "\n", sep="")


 cat("\n\n***** GOODNESS-OF-FIT MEASURES *****\n")
 cat("\nNUMBER OF SIMULATIONS = ", x$NSIM, "\n", sep="")

 cat("\nGEN. LOGISTIC        Z VALUE = ", x$Z[1], "\n", sep="")
 cat("GEN. EXTREME VALUE   Z VALUE = ", x$Z[2], "\n", sep="")
 cat("GEN. NORMAL          Z VALUE = ", x$Z[3], "\n", sep="")
 cat("PEARSON TYPE III     Z VALUE = ", x$Z[4], "\n", sep="")
 cat("GEN. PARETO          Z VALUE = ", x$Z[5], "\n", sep="")

 cat("\nPARAMETER ESTIMATES FOR DISTRIBUTIONS\n")
 cat("GEN. LOGISTIC        =", x$PARA[1:3,1], "\n", sep=" ")
 cat("GEN. EXTREME VALUE   =", x$PARA[1:3,2], "\n", sep=" ")
 cat("GEN. NORMAL          =", x$PARA[1:3,3], "\n", sep=" ")
 cat("PEARSON TYPE III     =", x$PARA[1:3,4], "\n", sep=" ")
 cat("GEN. PARETO          =", x$PARA[1:3,5], "\n", sep=" ")
 cat("WAKEBY               =", x$PARA[1:5,6], "\n", sep=" ")
}
