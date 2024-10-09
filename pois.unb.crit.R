# Defining the control limits and randomization probabilities of ARL unbiased c-chart

# function G
G <- Vectorize( function(x, l0) sum( dpois(0:x, l0) * (0:x) )/l0 , "x" )


# function G^-1
iG1 <- Vectorize( function(g, l0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qpois(g, l0)
      if ( G(x, l0) < g  ) {
        while ( G(x, l0) < g ) x <- x + 1
      } else {
        while ( G(x, l0) >= g & x > 0 ) x <- x - 1
        if ( G(x, l0) < g ) x <- x + 1
      }
    }
  }
  x
}, "g")


# function ~G^-1
iG2 <- Vectorize( function(g, l0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qpois(g, l0)
      if ( G(x, l0) <= g  ) {
        while ( G(x, l0) <= g ) x <- x + 1
      } else {
        while ( G(x, l0) > g & x > 0 ) x <- x - 1
        if ( G(x, l0) <= g ) x <- x + 1
      }
    }
  }
  x
}, "g")


# function F
F <- function(x, l0) ppois(x, l0)


# function F^-1
iF1 <- function(f, l0) qpois(f, l0)


# function ~F^-1
iF2 <- Vectorize( function(f, l0) {
  if ( f < 0 | f > 1 ) {
    stop("error")
  } else {
    if ( f == 1 ) {
      x <- Inf
    } else {
      x <- qpois(f, l0)
      if ( F(x, l0) <= f  ) {
        while ( F(x, l0) <= f ) x <- x + 1
      } else {
        while ( F(x, l0) > f & x > 0 ) x <- x - 1
        if ( F(x, l0) <= f ) x <- x + 1
      }
    }
  }
  x
}, "f")


# bounds
c.bounds <- function(alpha, l0) {
  Lmax1 <- iF2(alpha, l0)
  Lmax2 <- iG2(alpha, l0)
  Umin1 <- iF1(1-alpha, l0)
  Umin2 <- iG1(1-alpha, l0)
  Lmax  <- min(Lmax1, Lmax2)
  Umin  <- max(Umin1, Umin2)

  Lmin1 <- iF1(max(F(Umin-1,l0)-1+alpha,0), l0)
  Lmin2 <- iG1(max(G(Umin-1,l0)-1+alpha,0), l0)
  Umax1 <- iF2(min(F(Lmax,l0)+1-alpha,1), l0)
  Umax2 <- iG2(min(G(Lmax,l0)+1-alpha,1), l0)
  Lmin  <- max(Lmin1, Lmin2)
  Umax  <- min(Umax1, Umax2)

  data.frame(l0, Lmin, Lmax, Umin, Umax)
}

# randomization constants
gammas <- function(cL, cU, l0, alpha) {
  gU <- ( l0*( alpha - G(cL-1,l0) - 1+G(cU,l0) ) - cL*( alpha - F(cL-1,l0) - 1+F(cU,l0) ) ) / (dpois(cU,l0) * (cU-cL))
  gL <- ( cU*( alpha - F(cL-1,l0) - 1+F(cU,l0) ) - l0*( alpha - G(cL-1,l0) - 1+G(cU,l0) ) ) / (dpois(cL,l0) * (cU-cL))
  gg <- data.frame(gL, gU)
  gg
}


# identify admissible constants
search.GandC <- function(l0, alpha, OUTPUT=FALSE) {
  CB <- c.bounds(alpha, l0)
  results <- NULL
  for ( cL in CB$Lmin:CB$Lmax ) {
    for ( cU in CB$Umin:CB$Umax ) {
      GG <- gammas(cL, cU, l0, alpha)
      good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
      if ( OUTPUT ) cat(paste("cL =", cL, ", cU =", cU, ", gL =", GG[1], ", gU =", GG[2], "\n"))
      if ( good ) {
        results <- rbind(results, data.frame(l0, CB$Lmin, CB$Lmax, cL, CB$Umin, CB$Umax, cU, gL=GG[1], gU=GG[2]))
        break
      }
    }
  }
  results
}

#' Control Limits and Randomization Probabilities of the ARL unbiased c-chart function
#'
#' @param mu0 target mean
#' @param arl0 average run length in control or supremum of the ARL function
#' @param delta shift in the process mean
#'
#' @return it returns four values, respectively, the lower and upper control limits (LCL and UCL) the lower and upper randomization probabilities
#' @export
#'
#' @author something
#'
#' @references Paulino, S., Morais, M.C., and Knoth, S. (2016). An ARL-unbiased c-chart. Quality and Reliability Engineering International 32, 2847-2858.
#' @references Morais, M.C. and Knoth, S. (2020b). c-Charts. Wiley StatsRef: Statistics Reference Online.
#' @examples ## In control process, Morais 2016, Table II
#' mu0 <- 19
#' arl0 <- 370.4
#' delta <- 0
#' pois.unb.crit(mu0,arl0,delta)
#'
pois.unb.crit <- function(mu0,arl0,delta){
  results <- search.GandC(mu0+delta, (arl0)^(-1))
  cL4  <- results$cL
  cU4  <- results$cU
  gL4  <- results$gL
  gU4  <- results$gU
  result <- paste("LCL:", cL4, " ", "UCL:", cU4," ", "gl:", round(gL4,4), " ", "gU:", round(gU4,4))
  return(result)
}

mu0 <- 19
arl0 <- 370.4
delta <- 0
pois.unb.crit(mu0,arl0,delta)

####### PARA TABELA TESE #######

# Define your pois.unb.crit function
pois.unb.crit <- function(mu0, arl0, delta) {
  results <- search.GandC(mu0 + delta, (arl0)^(-1))
  cL4 <- results$cL
  cU4 <- results$cU
  gL4 <- results$gL
  gU4 <- results$gU
  result <- paste("mu0:", mu0, "(LCL, UCL):", cL4, cU4, "(gL, gU):", round(gL4, 6), round(gU4, 6))
  return(result)
}

# Vector of mu0 values
mu0_values <- c(0.05, 0.1, 0.5, 1, 5, 10, 15, 20, 25, 30)

# Specify arl0 and delta values
arl0 <- 370.4
delta <- 0 

# Loop through each mu0 value and print the results
for (mu0 in mu0_values) {
  print(pois.unb.crit(mu0, arl0, delta))
}


