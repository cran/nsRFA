R2 <- function(x, y) {

 # INPUT
 # x = observed values
 # y = estimated values

 n <- length(x)
 if (!length(y)==n) stop("R2: x and y must have the same length")
 SST <- sum(x^2) - n*mean(x)^2
 SSRes <- sum((x-y)^2)
 R2 <- 1-SSRes/SST

 return(R2)
}


# ------------------------------------------------------------------ #

RMSE <- function(x, y) {

 # INPUT
 # x = observed values
 # y = estimated values

 n <- length(x)
 if (!length(y)==n) stop("RMSE: x and y must have the same length")
 res <- x-y
 RMSE <- sqrt(sum((res)^2)/n)

 return(RMSE)
}


# ------------------------------------------------------------------ #

MAE <- function(x, y) {

 # INPUT
 # INPUT
 # x = observed values
 # y = estimated values

 n <- length(x)
 if (!length(y)==n) stop("MAE: x and y must have the same length")
 res <- x-y
 MAE <- sum(abs(res))/n

 return(MAE)
}


# ------------------------------------------------------------------ #

RMSEP <- function(x, y) {

 # INPUT
 # x = observed values
 # y = estimated values

 n <- length(x)
 if (!length(y)==n) stop("RMSE: x and y must have the same length")
 res <- (x-y)/x
 RMSEP <- sqrt(sum((res)^2)/n)

 return(RMSEP)
}


# ------------------------------------------------------------------ #

MAEP <- function(x, y) {

 # INPUT
 # INPUT
 # x = observed values
 # y = estimated values

 n <- length(x)
 if (!length(y)==n) stop("MAE: x and y must have the same length")
 res <- (x-y)/x
 MAEP <- sum(abs(res))/n

 return(MAEP)
}

