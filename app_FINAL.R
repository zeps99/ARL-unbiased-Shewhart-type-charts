library(shiny)
library(shinycssloaders)
library(DT)  # For rendering tables

#############################
#    *****  *****  ***      #
#      *      *    *   *    #
#      *      *    *   *    #
#    *****  *****  ***      #
#############################


#########################
### c-chart (Poisson) ###
#########################

# function G
Gc <- Vectorize( function(x, l0) sum( dpois(0:x, l0) * (0:x) )/l0 , "x" )

# function G^-1
iG1c <- Vectorize( function(g, l0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qpois(g, l0)
      if ( Gc(x, l0) < g  ) {
        while ( Gc(x, l0) < g ) x <- x + 1
      } else {
        while ( Gc(x, l0) >= g & x > 0 ) x <- x - 1
        if ( Gc(x, l0) < g ) x <- x + 1
      }
    }
  }
  x
}, "g")


# function ~G^-1
iG2c <- Vectorize( function(g, l0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qpois(g, l0)
      if ( Gc(x, l0) <= g  ) {
        while ( Gc(x, l0) <= g ) x <- x + 1
      } else {
        while ( Gc(x, l0) > g & x > 0 ) x <- x - 1
        if ( Gc(x, l0) <= g ) x <- x + 1
      }
    }
  }
  x
}, "g")


# function F
Fc <- function(x, l0) ppois(x, l0)


# function F^-1
iF1c <- function(f, l0) qpois(f, l0)


# function ~F^-1
iF2c <- Vectorize( function(f, l0) {
  if ( f < 0 | f > 1 ) {
    stop("error")
  } else {
    if ( f == 1 ) {
      x <- Inf
    } else {
      x <- qpois(f, l0)
      if ( Fc(x, l0) <= f  ) {
        while ( Fc(x, l0) <= f ) x <- x + 1
      } else {
        while ( Fc(x, l0) > f & x > 0 ) x <- x - 1
        if ( Fc(x, l0) <= f ) x <- x + 1
      }
    }
  }
  x
}, "f")


# bounds
c.boundsc <- function(alpha, l0) {
  Lmax1 <- iF2c(alpha, l0)
  Lmax2 <- iG2c(alpha, l0)
  Umin1 <- iF1c(1-alpha, l0)
  Umin2 <- iG1c(1-alpha, l0)
  Lmax  <- min(Lmax1, Lmax2)
  Umin  <- max(Umin1, Umin2)
  
  Lmin1 <- iF1c(max(Fc(Umin-1,l0)-1+alpha,0), l0)
  Lmin2 <- iG1c(max(Gc(Umin-1,l0)-1+alpha,0), l0)
  Umax1 <- iF2c(min(Fc(Lmax,l0)+1-alpha,1), l0)
  Umax2 <- iG2c(min(Gc(Lmax,l0)+1-alpha,1), l0)
  Lmin  <- max(Lmin1, Lmin2)
  Umax  <- min(Umax1, Umax2)
  
  data.frame(l0, Lmin, Lmax, Umin, Umax)
}

# randomization constants
gammasc <- function(cL, cU, l0, alpha) {
  gU <- ( l0*( alpha - Gc(cL-1,l0) - 1+Gc(cU,l0) ) - cL*( alpha - Fc(cL-1,l0) - 1+Fc(cU,l0) ) ) / (dpois(cU,l0) * (cU-cL))
  gL <- ( cU*( alpha - Fc(cL-1,l0) - 1+Fc(cU,l0) ) - l0*( alpha - Gc(cL-1,l0) - 1+Gc(cU,l0) ) ) / (dpois(cL,l0) * (cU-cL))
  gg <- data.frame(gL, gU)
  gg
}


# identify admissible constants
search.GandC_c <- function(l0, alpha, OUTPUT=FALSE) {
  CB <- c.boundsc(alpha, l0)
  results <- NULL
  for ( cL in CB$Lmin:CB$Lmax ) {
    for ( cU in CB$Umin:CB$Umax ) {
      GG <- gammasc(cL, cU, l0, alpha)
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


# ARL function with randomization
arl2c <- Vectorize(function(cL, cU, gL, gU, lambda) {
  aL <- Fc(cL-1, lambda) + gL*dpois(cL, lambda)
  aU <- 1 - Fc(cU, lambda) + gU*dpois(cU, lambda)
  L  <- 1 / ( aL + aU )
  L
}, "lambda")


#############################
#### np-chart (Binomial) ####
#############################

Fnp <- function(x,n,p0) pbinom(x,n,p0)

Pnp <- function(y,n,p0) dbinom(y,n,p0)

initialcL_ <- function(alpha, n, p0) {
  if (alpha < 0 | alpha > 1) {
    stop("error")
  } else {
    if (alpha == 1) {
      x <- Inf
    } else {
      inv <- qbinom(alpha, n, p0)
      if (pbinom(inv, n, p0) > alpha) {
        x <- inv
      } else {
        x <- inv + 1
      }
    }
  }
  x
}

initialcL <- Vectorize("initialcL_", "alpha")

# Randomization probabilities

gammasnp <- function(cL,cU,n,p0,alpha) {
  a<-Pnp(cL,n,p0)
  b<-Pnp(cU,n,p0)
  c<-(cL*Pnp(cL,n,p0))
  d<-(cU*Pnp(cU,n,p0))
  e<-(alpha-1+Fnp(cU,n,p0)-Fnp(cL-1,n,p0))
  f<-(n*p0*(alpha-1+Fnp(cU-1,n-1,p0)-Fnp(cL-2,n-1,p0)))
  gL<-((d*e-b*f)/(a*d-b*c))
  gU<-((a*f-c*e)/(a*d-b*c))
  gg <- data.frame(gL, gU)
  gg
}

search.GandC_np <- function(n,p0,alpha) {
  results <- NULL
  cL <- initialcL(alpha,n,p0)
  q <- (1-alpha+Fnp(cL-1,n,p0))
  cU <- qbinom(q,n,p0)
  GG <- gammasnp(cL,cU,n,p0,alpha)
  good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
  cL <- (cL+1)
  while (cL>0) {
    cL <- (cL-1)
    q <- (1-alpha+Fnp(cL-1,n,p0))
    cU <- qbinom(q,n,p0)
    GG <- gammasnp(cL,cU,n,p0,alpha)
    good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
    if ( good ) {
      sinal <- emitirsinal(n,p0,cU,cL,GG[1],GG[2])
      results <- data.frame(n,p0, cL, cU, gL=GG[1], gU=GG[2],prob=sinal[1],ARL=sinal[2])
      break
    } else {
      cU<-(cU+1)
      GG <- gammasnp(cL,cU,n,p0,alpha)
      good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
    }
    if ( good ) {
      sinal <- emitirsinal(n,p0,cU,cL,GG[1],GG[2])
      results <- data.frame(n,p0, cL, cU, gL=GG[1], gU=GG[2],prob=sinal[1],ARL=sinal[2])
      break
    }
  }
  results
}

# ARL calculation

emitirsinal <- function(n,p0,cU,cL,gL,gU) {
  prob <- (1-Fnp(cU,n,p0)+Fnp(cL-1,n,p0)+gL*Pnp(cL,n,p0)+gU*Pnp(cU,n,p0))
  arl <- (1/prob)
  emissao <- data.frame(prob,arl)
  emissao
}

# arl function

arl2np <- function(n,p0,cU,cL,gL,gU) {
  prob <- (1-Fnp(cU,n,p0)+Fnp(cL-1,n,p0)+gL*Pnp(cL,n,p0)+gU*Pnp(cU,n,p0))
  arl <- (1/prob)
  return(arl)
}


binom.unb.arl <- function(n,p0,arl0, delta) {
  alpha <- 1/arl0
  
  result <- search.GandC_np(n,p0,alpha)
  
  cL1 <- result$cL
  cU1 <- result$cU
  gL1 <- result$gL
  gU1 <- result$gU
  
  # p = p0 + delta
  ARL <- arl2np(n,p0+delta,cU1,cL1,gL1,gU1)
  return(ARL)
}


################################
### g-chart (Geometric) NOVO ###
################################

# Function P, F and etc
Pg <- function(x, p0) {(1-p0)^x * p0} 

Fg <- function(x, p0) {1 - (1-p0)^(x+1)}

iF1g <- function(f, p0) {
  if (f == p0) {
    stop("error")
  } else {
    ceiling(log(1-f)/log(1-p0)) -1
  }
} 

iF2g <- function(f, p0) {
  if (f == p0) {
    stop("error")
  } else {
    ceiling(log(1-f)/log(1-p0)) -1
  }
}


# Function G and etc
G_g <- function(x, p0) {1 - (1-p0)^x * (1+p0*x)}

Gg <- Vectorize("G_g", "x")

iG1_g <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( Gg(x, p0) < g  ) {
        while ( Gg(x, p0) < g ) x <- x + 1
      } else {
        while ( Gg(x, p0) >= g & x > 0 ) x <- x - 1
        if ( Gg(x, p0) < g ) x <- x + 1
      }
    }
  }
  x
}
iG1g <- Vectorize("iG1_g", "g")

iG2_g <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( Gg(x, p0) <= g  ) {
        while ( Gg(x, p0) <= g ) x <- x + 1
      } else {
        while ( Gg(x, p0) > g & x > 0 ) x <- x - 1
        if ( Gg(x, p0) <= g ) x <- x + 1
      }
    }
  }
  x
}
iG2g <- Vectorize("iG2_g", "g")

# Bounds

c.boundsg <- function(alpha, p0) {
  Lmax1 <- iF2g(alpha, p0)
  Lmax2 <- iG2g(alpha, p0)
  Umin1 <- iF1g(1-alpha, p0)
  Umin2 <- iG1g(1-alpha, p0)
  Lmax <- min(Lmax1, Lmax2)
  Umin <- max(Umin1, Umin2)
  
  Lmin1 <- iF1g(max(Fg(Umin-1,p0)-1+alpha,0), p0)
  Lmin2 <- iG1g(max(Gg(Umin-1,p0)-1+alpha,0), p0)
  Umax1 <- iF2g(min(Fg(Lmax,p0)+1-alpha,1), p0)
  Umax2 <- iG2g(min(Gg(Lmax,p0)+1-alpha,1), p0)
  Lmin <- max(Lmin1, Lmin2)
  Umax <- min(Umax1, Umax2)
  
  data.frame(p0, Lmin, Lmax, Umin, Umax)
}


# Randomization Constants
gammasg <- function(cL, cU, p0, alpha) {
  a <- (1-p0)^cL * p0
  b <- (1-p0)^cU * p0
  c <- (cL*(1-p0)^cL * p0)
  d <- (cU*(1-p0)^cU * p0)
  e <- (alpha-1+(1-p0)^(cL)-(1-p0)^(1+cU))
  f <- ((alpha-1)*((1/p0)-1) + (1/p0)*((1-p0)^(cL)*(1+(-1+cL)*p0) + (1-p0)^(cU)*(-1+p0)*(1+p0*cU)))
  gL <- ((d*e-b*f)/(a*d-b*c))
  gU <- ((a*f-c*e)/(a*d-b*c))
  gg <- data.frame(gL, gU)
  gg
}

# Identify Admissible Constants
search.GandC_g <- function(p0, alpha, OUTPUT=FALSE) {
  CB <- c.boundsg(alpha, p0)
  results <- NULL
  
  # Create vectors of possible values for cL and cU
  cL_values <- CB$Lmin:CB$Lmax
  cU_values <- CB$Umin:CB$Umax
  
  newlist <- list()
  
  for (i in cL_values) {
    for (j in cU_values) {
      newlist <- c(newlist, list(c(i, j)))
    }
  }
  
  # Number of possible combinations
  num_combinations <- length(newlist)
  
  while (is.null(results)) {
    
    # Randomly sample indices of the pair (cL,cU) without replacement
    random_index <- sample(num_combinations, 1, replace = FALSE)
    random_pair <- newlist[[random_index]]
    
    cL <- random_pair[1]
    cU <- random_pair[2]
    
    # Remove the selected pair from newlist
    newlist <- newlist[-random_index]
    
    # Recalculate number of combinations after removing one pair
    num_combinations <- length(newlist)
    
    GG <- gammasg(cL, cU, p0, alpha)
    if (0 <= GG[2] & GG[2] <= 1) {
      if (0 <= GG[1] & GG[1] <= 1) {
        if (OUTPUT) cat(paste("cL =", cL, ", cU =", cU, ", gL =", GG[1], ", gU =", GG[2], "\n"))
        results <- rbind(results, data.frame(p0, CB$Lmin, CB$Lmax, cL, CB$Umin, CB$Umax, cU, gL = GG[1], gU = GG[2]))
        break
      }
    }
  }
  results
}

geom.unb.arl <- function(l0,arl0, l){
  results <- search.GandC_g(l0, (arl0)^(-1))
  cL4  <- results$cL
  cU4  <- results$cU
  gL4  <- results$gL
  gU4  <- results$gU
  
  # l <- l0 + delta
  aL <- (1 - (1-l)^(cL4)) + gL4*((1-l)^cL4 * l) 
  aU <- (1-l)^(cU4+1) + gU4*((1-l)^cU4 * l)
  ARL  <- 1 / ( aL + aU )
  ARL
}

################################ 
#     *    ***      *  *  *    #
#   *   *  *   *   * * *   *   #
#   *****  ***     *   *   *   #
#   *   *  *   *    *  *  *    #
################################

#######################
### Poisson INAR(1) ###
#######################

# iid
iid.arl <- function(L, U, mu) {
  pC  <- ppois(U, mu) - ppois(L-1, mu)
  ARL <- 1/(1-pC)
  ARL
}

# INAR(1) w/o randomization, ARL
inar1.arl <- function(L, U, lambda, beta) {
  i    <- L:U
  d    <- U-L+1
  qij_ <- function(i, j) {
    m  <- 0 : min(i,j)
    sum ( dbinom(m, i, beta) * dpois(j-m, lambda) )
  }
  qij  <- Vectorize(qij_)
  Q    <- outer(i, i, qij)
  one  <- rep(1, d)
  I    <- diag(1, d)
  cARL <- solve(I-Q, one)
  ARL  <- 1 + sum( dpois(i, lambda/(1-beta) ) * cARL )
  ARL
}

# INAR(1) w/ randomization, ARL
inar1.arl2 <- function(L, U, gL, gU, lambda, beta) {
  i     <- L:U
  d     <- U-L+1
  qij_  <- function(i, j) {
    m   <- 0 : min(i,j)
    sum ( dbinom(m, i, beta) * dpois(j-m, lambda) )
  }
  qij   <- Vectorize(qij_)
  Q     <- outer(i, i, qij)
  Q[,1] <- (1-gL) * Q[,1]
  Q[,d] <- (1-gU) * Q[,d]
  one   <- rep(1, d)
  I     <- diag(1, d)
  cARL  <- solve(I-Q, one)
  DP    <- dpois(i, lambda/(1-beta))
  DP[1] <- (1-gL) * DP[1]
  DP[d] <- (1-gU) * DP[d]
  ARL   <- 1 + sum( DP * cARL )
  ARL
}

# INAR(1) w/o randomization, L given, U searched for
inar1.get.U <- function(L, ARL0, lambda, beta, OUTPUT=FALSE) {
  U   <- ceiling( lambda / ( 1 - beta ) )
  ARL <- inar1.arl(L, U, lambda, beta)
  if ( OUTPUT ) cat(paste(U, "\t", ARL, "\t", ARL0, "\n"))
  dARL <- 1
  while ( ARL < ARL0 & dARL > 1e-6 ) {
    ARL.old <- ARL
    U   <- U + 1
    ARL <- inar1.arl(L, U, lambda, beta)
    dARL <- abs(ARL - ARL.old)
    if ( OUTPUT ) cat(paste(U, "\t", ARL, "\t", ARL0, "\n"))
  }
  if ( ARL < ARL0 ) U <- NA
  U
}

# INAR(1) w/o randomization, (L,U) with minimal L so that ARL = ARL(lambda) is increasing at given lambda
inar1.get.UL <- function(ARL0, lambda, beta, OUTPUT=FALSE) {
  L <- 0 
  U <- inar1.get.U(L, ARL0, lambda, beta, OUTPUT=OUTPUT) 
  lARL <- inar1.arl(L, U, lambda-1e-3, beta)
  ARL  <- inar1.arl(L, U, lambda, beta)
  while ( ARL < lARL ) {
    L    <- L + 1
    U    <- inar1.get.U(L, ARL0, lambda, beta, OUTPUT=OUTPUT) 
    if ( is.na(U) ) {
      L    <- L - 1
      U    <- inar1.get.U(L, ARL0, lambda, beta, OUTPUT=OUTPUT) 
      break
    }
    lARL <- inar1.arl(L, U, lambda-1e-2, beta)
    ARL  <- inar1.arl(L, U, lambda, beta)
    if ( OUTPUT ) cat(paste(L, "\t", U, "\t", lARL, "\t", ARL, "\t\t", ARL0, "\n"))
  }
  c(L, U)
}

# INAR(1) w/ randomization, L, U, gL given, gU searched for
inar1.get.gU <- function(L, U, gL, ARL0, lambda, beta, OUTPUT=FALSE) {
  minARL <- inar1.arl2(L, U, gL, 1, lambda, beta)
  maxARL <- inar1.arl2(L, U, gL, 0, lambda, beta)
  if ( OUTPUT ) cat(paste(minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gU1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gU1 > 0.1 ) {
      gU2  <- gU1
      ARL2 <- ARL1
      gU1  <- gU1 - .1
      ARL1 <- inar1.arl2(L, U, gL, gU1, lambda, beta)
    }
    if ( gU1 < .05 ) {
      gU1  <- 0.05
      ARL1 <- inar1.arl2(L, U, gL, gU1, lambda, beta)
    }
    # secant rule
    gU.error <- 1
    L.error  <- 1
    while ( gU.error > 1e-10  &  L.error > 1e-9 ) {
      gU3  <- gU1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gU2 - gU1)
      gU3 <- max(0, gU3)
      gU3 <- min(1, gU3)
      ARL3 <- inar1.arl2(L, U, gL, gU3, lambda, beta)
      if ( OUTPUT ) cat(paste(gU3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gU1 <- gU2; gU2 <- gU3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gU.error <- abs(gU2 - gU1)
    }
  } else {
    gU3 <- NA                      # old
    if ( minARL > ARL0 ) gU3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gU3 <- -1 # maximal reachable ARL too small
  }
  gU3
}

# INAR(1) w/ randomization, L, U, gU given, gL searched for
inar1.get.gL <- function(L, U, gU, ARL0, lambda, beta, OUTPUT=FALSE) {
  minARL <- inar1.arl2(L, U, 1, gU, lambda, beta)
  maxARL <- inar1.arl2(L, U, 0, gU, lambda, beta)
  if ( OUTPUT ) cat(paste(minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gL1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gL1 > .1 ) {
      gL2  <- gL1
      ARL2 <- ARL1
      gL1  <- gL1 - .1
      ARL1 <- inar1.arl2(L, U, gL1, gU, lambda, beta)
    }
    if ( gL1 < 0.05 ) {
      gL1  <- 0.05
      ARL1 <- inar1.arl2(L, U, gL1, gU, lambda, beta)
    }
    # secant rule
    gL.error <- 1
    L.error  <- 1
    while ( gL.error > 1e-10  &  L.error > 1e-9 ) {
      gL3  <- gL1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gL2 - gL1)
      gL3 <- max(0, gL3)
      gL3 <- min(1, gL3)
      ARL3 <- inar1.arl2(L, U, gL3, gU, lambda, beta)
      if ( OUTPUT ) cat(paste(gL3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gL1 <- gL2; gL2 <- gL3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gL.error <- abs(gL2 - gL1)
    }
  } else {
    gL3 <- NA                      # old
    if ( minARL > ARL0 ) gL3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gL3 <- -1 # maximal reachable ARL too small
  }
  gL3
}

# some helpers 
inar1.min.gL <- function(L, U, ARL0, lambda, beta, OUTPUT=FALSE) {
  gL  <- 0
  gU  <- inar1.get.gU(L, U, gL, ARL0, lambda, beta, OUTPUT=OUTPUT)
  if ( -1.5 < gU & gU < 0 ) gL <- NA
  if ( gU < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gU < -1.5 & gL < 1-1e-10 ) {
          gL <- gL + 10^(-dig)
          gU <- inar1.get.gU(L, U, gL, ARL0, lambda, beta, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("gL =", gL, ",\tgU =", gU, "\n"))
        }
        if ( gU < -1.5 ) {
          gL <- NA
          break
        }
      } else {
        while ( -1.5 < gU  &  gL > 1e-10 ) {
          gL <- gL - 10^(-dig)
          gU <- inar1.get.gU(L, U, gL, ARL0, lambda, beta, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("gL =", gL, ",\tgU =", gU, "\n"))
        }
      }
    }
  }
  gL
}

inar1.min.gU <- function(L, U, ARL0, lambda, beta, OUTPUT=FALSE) {
  gU  <- 0
  gL  <- inar1.get.gL(L, U, gU, ARL0, lambda, beta, OUTPUT=OUTPUT)
  if ( -1.5 < gL & gL < 0 ) gU <- NA
  if ( gL < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gL < -1.5 & gU < 1 ) {
          gU <- gU + 10^(-dig)
          gL <- inar1.get.gL(L, U, gU, ARL0, lambda, beta, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("gU =", gU, ",\tgL =", gL, "\n"))
        }
        if ( gL < -1.5 ) {
          gU <- NA
          break
        }
      } else {
        while ( -1.5 < gL & gU > 0 ) {
          gU <- gU - 10^(-dig)
          gL <- inar1.get.gL(L, U, gU, ARL0, lambda, beta, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("gU =", gU, ",\tgL =", gL, "\n"))
        }
      }
    }
  }
  gU
}

# INAR(1) w/ randomization, get all 4!
inar1.get.UL2 <- function(ARL0, lambda, beta, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4) {
  if ( target=="lambda" ) {
    l1 <- lambda - eps
    l2 <- lambda + eps
    b1 <- beta
    b2 <- beta
  }
  if ( target=="beta" ) {
    l1 <- lambda
    l2 <- lambda
    b1 <- beta - eps
    b2 <- beta + eps
  }
  LU  <- inar1.get.UL(ARL0, lambda, beta, OUTPUT=OUTPUT)
  L2  <- LU[1]
  U2  <- LU[2]
  if ( L2 > 0 ) {
    L1  <- L2 - 1
    U1  <- inar1.get.U(L1, ARL0, lambda, beta, OUTPUT=OUTPUT)
  } else {
    L1 <- 0
    U1 <- U2
    U2 <- U1 + 1
  }
  #if ( U1==U2 ) U2 <- U2 + 1
  U2 <- U2 + 1
  
  for ( L in L1:L2 ) {
    for ( U in U1:U2 ) {
      gL1  <- inar1.min.gL(L, U, ARL0, lambda, beta, OUTPUT=OUTPUT)
      if ( is.na(gL1) ) next
      gU1  <- inar1.get.gU(L, U, gL1, ARL0, lambda, beta)
      lARL1 <- inar1.arl2(L, U, gL1, gU1, l1, b1)
      rARL1 <- inar1.arl2(L, U, gL1, gU1, l2, b2)
      dratio1 <- ( rARL1 - lARL1 ) / ( 2*eps )
      
      gU2 <- inar1.min.gU(L, U, ARL0, lambda, beta, OUTPUT=OUTPUT)
      if ( is.na(gU2) ) next
      gL2 <- inar1.get.gL(L, U, gU2, ARL0, lambda, beta)
      lARL2 <- inar1.arl2(L, U, gL2, gU2, l1, b1)
      rARL2 <- inar1.arl2(L, U, gL2, gU2, l2, b2)
      dratio2 <- ( rARL2 - lARL2 ) / ( 2*eps )
      
      if ( OUTPUT ) cat( paste("L =", L, ",\t U =", U, ",\tdr1 =", dratio1, ",\tdr2 =", dratio2, "\n") )
      
      if ( dratio1 * dratio2 < 0 ) {
        # secant rule
        gL.error <- 1
        dr.error <- 1
        while ( gL.error > 1e-10  &  dr.error > delta ) {
          gL3 <- gL1 + (0 - dratio1)/(dratio2 - dratio1)*(gL2 - gL1)
          gU3 <- inar1.get.gU(L, U, gL3, ARL0, lambda, beta)
          lARL3 <- inar1.arl2(L, U, gL3, gU3, l1, b1)
          rARL3 <- inar1.arl2(L, U, gL3, gU3, l2, b2)
          dratio3 <- ( rARL3 - lARL3 ) / ( 2*eps )
          if ( OUTPUT ) cat(paste(gL3, "\t", dratio3, "\n"))
          gL1 <- gL2; gL2 <- gL3
          dratio1 <- dratio2; dratio2 <- dratio3
          dr.error  <- abs(dratio2)
          gL.error <- abs(gL2 - gL1)
        }
        L0 <- L
        U0 <- U
        L <- L2 + 1
        U <- U2 + 1
      }
      if ( U > U2 ) break
    }
    if ( L > L2 ) break
  }
  data.frame(L=L0, gL=gL3, U=U0, gU=gU3)
}

########################
#### binomial AR(1) ####
########################

# ARL of an np-chart (without randomization) falsely assuming independent
# output when dealing with AR(1) binomial counts 
ar1.arl <- function(L, U, n, p, rho) {
  if (rho<max(-p/(1-p),-(1-p)/p) || rho>1) break
  i    <- L:U
  d    <- U-L+1
  qij_ <- function(i, j) {
    m  <- max(0,j+i-n) : min(i,j)
    beta <- p*(1-rho)
    alpha <- beta+rho
    sum ( dbinom(m, i, alpha) * dbinom(j-m, n-i, beta) )
  }
  qij  <- Vectorize(qij_)
  Q    <- outer(i, i, qij)
  one  <- rep(1, d)
  I    <- diag(1, d)
  cARL <- solve(I-Q, one)
  ARL  <- 1 + sum( dbinom(i, n, p ) * cARL )
  ARL
}

# ARL of an np-chart with randomization when dealing with AR(1) binomial counts 
ar1.arl2 <- function(L, U, gL, gU, n, p, rho) {
  if (rho<max(-p/(1-p),-(1-p)/p) || rho>1) return(NA)   ## alterei de break para return(NA)
  i    <- L:U
  d    <- U-L+1
  qij_ <- function(i, j) {
    m  <- max(0,j+i-n) : min(i,j)
    beta <- p*(1-rho)
    alpha <- beta+rho
    sum ( dbinom(m, i, alpha) * dbinom(j-m, n-i, beta) )
  }
  qij   <- Vectorize(qij_)
  Q     <- outer(i, i, qij)
  Q[,1] <- (1-gL) * Q[,1]
  Q[,d] <- (1-gU) * Q[,d]
  one   <- rep(1, d)
  I     <- diag(1, d)
  cARL  <- solve(I-Q, one)
  DP    <- dbinom(i, n, p )
  DP[1] <- (1-gL) * DP[1]
  DP[d] <- (1-gU) * DP[d]
  ARL   <- 1 + sum( DP * cARL )
  ARL
}

# np-chart without randomization for binomial AR(1) counts
# L given, searched for U
ar1.get.U <- function(L, ARL0, n, p, rho, OUTPUT=FALSE) {
  U   <-  ceiling(n*p)
  ARL <- ar1.arl(L, U, n, p, rho)
  if ( OUTPUT ) cat(paste(L, "\t", U, "\t", ARL, "\t", ARL0, "\n"))
  dARL <- 1
  while ( ARL < ARL0 & dARL > 1e-6 & U<n) {
    ARL.old <- ARL
    U   <- U + 1
    ARL <- ar1.arl(L, U, n, p, rho)
    dARL <- abs(ARL - ARL.old)
    if ( OUTPUT ) cat(paste(L, "\t", U, "\t", ARL, "\t", ARL0, "\n"))
  }
  if ( ARL < ARL0 ) U <- NA
  U
}

# np-chart without randomization for binomial AR(1) counts
# (L, U) with minimal L so that ARL = ARL(p) is  increasing 
# at given p
ar1.get.UL <- function(ARL0, n, p, rho, OUTPUT=FALSE) {
  L <- 0 
  U <- ar1.get.U(L, ARL0, n, p, rho, OUTPUT=OUTPUT) 
  lARL <- ar1.arl(L, U, n, p-1e-3, rho)
  ARL  <- ar1.arl(L, U, n, p, rho)
  while ( ARL < lARL ) {
    L    <- L + 1
    U    <- ar1.get.U(L, ARL0, n, p, rho, OUTPUT=OUTPUT) 
    if ( is.na(U) ) {
      L    <- L - 1
      U    <- ar1.get.U(L, ARL0, n, p, rho, OUTPUT=OUTPUT) 
      break
    }
    lARL <- ar1.arl(L, U, n, p-1e-2, rho)
    ARL  <- ar1.arl(L, U, n, p, rho)
    if ( OUTPUT ) cat(paste(L, "\t", U, "\t", lARL, "\t", ARL, "\t\t", ARL0, "\n"))
  }
  c(L, U)
}

# np-chart with randomization for binomial AR(1) counts
# L, U, gL given, searched for gU
ar1.get.gU <- function(L, U, gL, ARL0, n, p, rho, OUTPUT=FALSE) {
  minARL <- ar1.arl2(L, U, gL, 1, n, p, rho)
  maxARL <- ar1.arl2(L, U, gL, 0, n, p, rho)
  if ( OUTPUT ) cat(paste(minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gU1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gU1 > 0.1 ) {
      gU2  <- gU1
      ARL2 <- ARL1
      gU1  <- gU1 - .1
      ARL1 <- ar1.arl2(L, U, gL, gU1, n, p, rho)
    }
    if ( gU1 < .05 ) {
      gU1  <- 0.05
      ARL1 <- ar1.arl2(L, U, gL, gU1, n, p, rho)
    }
    # secant rule
    gU.error <- 1
    L.error  <- 1
    while ( gU.error > 1e-10  &  L.error > 1e-9 ) {
      gU3  <- gU1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gU2 - gU1)
      gU3 <- max(0, gU3)
      gU3 <- min(1, gU3)
      ARL3 <- ar1.arl2(L, U, gL, gU3, n, p, rho)
      if ( OUTPUT ) cat(paste(gU3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gU1 <- gU2; gU2 <- gU3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gU.error <- abs(gU2 - gU1)
    }
  } else {
    gU3 <- NA                     
    if ( minARL > ARL0 ) gU3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gU3 <- -1 # maximal reachable ARL too small
  }
  gU3
}

# np-chart with randomization for binomial AR(1) counts
# L, U, gU given, searched for gL
ar1.get.gL <- function(L, U, gU, ARL0, n, p, rho, OUTPUT=FALSE) {
  minARL <- ar1.arl2(L, U, 1, gU, n, p, rho)
  maxARL <- ar1.arl2(L, U, 0, gU, n, p, rho)
  if ( OUTPUT ) cat(paste(minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gL1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gL1 > .1 ) {
      gL2  <- gL1
      ARL2 <- ARL1
      gL1  <- gL1 - .1
      ARL1 <- ar1.arl2(L, U, gL1, gU, n, p, rho)
    }
    if ( gL1 < 0.05 ) {
      gL1  <- 0.05
      ARL1 <- ar1.arl2(L, U, gL1, gU, n, p, rho)
    }
    # secant rule
    gL.error <- 1
    L.error  <- 1
    while ( gL.error > 1e-10  &  L.error > 1e-9 ) {
      gL3  <- gL1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gL2 - gL1)
      gL3 <- max(0, gL3)
      gL3 <- min(1, gL3)
      ARL3 <- ar1.arl2(L, U, gL3, gU, n, p, rho)
      if ( OUTPUT ) cat(paste(gL3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gL1 <- gL2; gL2 <- gL3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gL.error <- abs(gL2 - gL1)
    }
  } else {
    gL3 <- NA                    
    if ( minARL > ARL0 ) gL3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gL3 <- -1 # maximal reachable ARL too small
  }
  gL3
}

# Some helpers 
ar1.min.gL <- function(L, U, ARL0, n, p, rho, OUTPUT=FALSE) {
  gL  <- 0
  gU  <- ar1.get.gU(L, U, gL, ARL0, n, p, rho, OUTPUT=OUTPUT)
  if ( -1.5 < gU & gU < 0 ) gL <- NA
  if ( gU < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gU < -1.5 & gL < 1-1e-10 ) {
          gL <- gL + 10^(-dig)
          gU <- ar1.get.gU(L, U, gL, ARL0, n, p, rho, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("gL =", gL, ",\tgU =", gU, "\n"))
        }
        if ( gU < -1.5 ) {
          gL <- NA
          break
        }
      } else {
        while ( -1.5 < gU  &  gL > 1e-10 ) {
          gL <- gL - 10^(-dig)
          gU <- ar1.get.gU(L, U, gL, ARL0, n, p, rho, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("gL =", gL, ",\tgU =", gU, "\n"))
        }
      }
    }
  }
  gL
}

ar1.min.gU <- function(L, U, ARL0, n, p, rho, OUTPUT=FALSE) {
  gU  <- 0
  gL  <- ar1.get.gL(L, U, gU, ARL0, n, p, rho, OUTPUT=OUTPUT)
  if ( -1.5 < gL & gL < 0 ) gU <- NA
  if ( gL < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gL < -1.5 & gU < 1 ) {
          gU <- gU + 10^(-dig)
          gL <- ar1.get.gL(L, U, gU, ARL0, n, p, rho, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("gU =", gU, ",\tgL =", gL, "\n"))
        }
        if ( gL < -1.5 ) {
          gU <- NA
          break
        }
      } else {
        while ( -1.5 < gL & gU > 0 ) {
          gU <- gU - 10^(-dig)
          gL <- ar1.get.gL(L, U, gU, ARL0, n, p, rho, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("gU =", gU, ",\tgL =", gL, "\n"))
        }
      }
    }
  }
  gU
}

# np-chart with randomization for binomial AR(1) counts get 
# all 4 constants (L, U, gL, gU)
ar1.get.UL2 <- function(ARL0, n, p, rho, target="p", OUTPUT=FALSE, eps=1e-6, delta=1e-4) {
  if ( target=="p" ) {
    l1 <- p - eps
    l2 <- p + eps
    b1 <- rho
    b2 <- rho
  }
  if ( target=="rho" ) {
    l1 <- p 
    l2 <- p
    b1 <- rho - eps
    b2 <- rho + eps
  }
  LU  <- ar1.get.UL(ARL0, n, p, rho, OUTPUT=OUTPUT)
  L2  <- LU[1]
  U2  <- LU[2]
  if ( L2 > 0 ) {
    L1  <- L2 - 1
    U1  <- ar1.get.U(L1, ARL0, n, p, rho, OUTPUT=OUTPUT)
  } else {
    L1 <- 0
    U1 <- U2
    U2 <- U1 + 1
  }
  U2 <- U2 + 1
  
  for ( L in L1:L2 ) {
    for ( U in U1:U2 ) {
      gL1  <- ar1.min.gL(L, U, ARL0, n, p, rho, OUTPUT=OUTPUT)
      if ( is.na(gL1) ) next
      gU1  <- ar1.get.gU(L, U, gL1, ARL0, n, p, rho)
      lARL1 <- ar1.arl2(L, U, gL1, gU1, n, l1, b1)
      rARL1 <- ar1.arl2(L, U, gL1, gU1, n, l2, b2)
      dratio1 <- ( rARL1 - lARL1 ) / ( 2*eps )
      
      gU2 <- ar1.min.gU(L, U, ARL0, n, p, rho, OUTPUT=OUTPUT)
      if ( is.na(gU2) ) next
      gL2 <- ar1.get.gL(L, U, gU2, ARL0, n, p, rho)
      lARL2 <- ar1.arl2(L, U, gL2, gU2, n, l1, b1)
      rARL2 <- ar1.arl2(L, U, gL2, gU2, n, l2, b2)
      dratio2 <- ( rARL2 - lARL2 ) / ( 2*eps )
      
      if ( OUTPUT ) cat( paste("L =", L, ",\t U =", U, ",\tdr1 =", dratio1, ",\tdr2 =", dratio2, "\n") )
      
      if ( dratio1 * dratio2 < 0 ) {
        # secant rule
        gL.error <- 1
        dr.error <- 1
        while ( gL.error > 1e-10  &  dr.error > eps ) {
          gL3 <- gL1 + (0 - dratio1)/(dratio2 - dratio1)*(gL2 - gL1)
          gU3 <- ar1.get.gU(L, U, gL3, ARL0, n, p, rho)
          lARL3 <- ar1.arl2(L, U, gL3, gU3, n, l1, b1)
          rARL3 <- ar1.arl2(L, U, gL3, gU3, n, l2, b2)
          dratio3 <- ( rARL3 - lARL3 ) / ( 2*eps )
          if ( OUTPUT ) cat(paste(gL3, "\t", dratio3, "\n"))
          gL1 <- gL2; gL2 <- gL3
          dratio1 <- dratio2; dratio2 <- dratio3
          L.error  <- abs(dratio2)
          gL.error <- abs(gL2 - gL1)
        }
        L0 <- L
        U0 <- U
        L <- L2 + 1
        U <- U2 + 1
      }
      if ( U > U2 ) break
    }
    if ( L > L2 ) break
  }
  data.frame(L=L0, gL=gL3, U=U0, gU=gU3)
}

# Function to calculate the dynamic range of rho
calculate_rho_range <- function(p) {
  lower_bound <- max(-p / (1 - p), -(1 - p) / p)
  c(lower_bound, 1)
}


#########################
### geometric INAR(1) ###
#########################

ginar1.arl <- function(L, U, lambda, beta) {
  if (beta < 0 || beta > 1) break
  i    <- L:U
  d    <- U-L+1
  # qij_  <- function(i, j) {
  # 	m  <- 0 : min(i,j)
  # 	sum ( dbinom(m, i, beta) * pdf.innov(j-m, lambda, beta) )
  # }
  qij_  <- function(i, j) {
    if(j<=i){
      m <- 0 : j
      sum (dbinom(m, i, beta)*dgeom(i-m, lambda) * (1-beta)*(1-lambda)^(j-i)
           + beta*dbinom(j, i, beta)/(j+1))
    }
    else{
      m <- 0 : i
      sum (dbinom(m, i, beta)*dgeom(i-m, lambda) * (1-beta)*(1-lambda)^(j-i))
    }
  }
  qij  <- Vectorize(qij_)
  Q    <- outer(i, i, qij)
  one  <- rep(1, d)
  I    <- diag(1, d)
  cARL <- solve(I-Q, one)
  ARL  <- 1 + sum( dgeom(i,lambda) * cARL )
  ARL
}


# ARL of a modified chart with randomization for GINAR(1) counts

ginar1.arl2 <- function(L, U, gL, gU, lambda, beta) {
  i     <- L:U
  d     <- U-L+1
  qij_  <- function(i, j) {
    if(j<=i){
      m <- 0 : j
      sum (dbinom(m, i, beta)*dgeom(i-m, lambda) * (1-beta)*(1-lambda)^(j-i)
           + beta*dbinom(j, i, beta)/(j+1))
    }
    else{
      m <- 0 : i
      sum (dbinom(m, i, beta)*dgeom(i-m, lambda) * (1-beta)*(1-lambda)^(j-i))
    }
  }
  qij   <- Vectorize(qij_)
  Q     <- outer(i, i, qij)
  Q[,1] <- (1-gL) * Q[,1]
  Q[,d] <- (1-gU) * Q[,d]
  one   <- rep(1, d)
  I     <- diag(1, d)
  cARL  <- solve(I-Q, one)
  DP    <- dgeom(i, lambda)
  DP[1] <- (1-gL) * DP[1]
  DP[d] <- (1-gU) * DP[d]
  ARL   <- 1 + sum( DP * cARL )
  ARL
}


# Modified chart without randomization for GINAR(1) counts
# L given, U searched for

ginar1.get.U <- function(L, ARL0, lambda, beta, OUTPUT=FALSE) {
  #  U   <- ceiling( lambda / ( 1 - lambda ) )  # Mean of a geometric(lambda)
  U   <- qgeom(1-1/(2*ARL0), lambda)   # (1-1/(2*ARL0))^th quantile of a geom(lambda)
  #  U <- search.GandC(lambda, 1/ARL0)$cU  # A possible solution
  ARL <- ginar1.arl(L, U, lambda, beta)
  if ( OUTPUT ) cat(paste("ginar1.get.U ", U, "\t", ARL, "\t", ARL0, "\n"))
  dARL <- 1
  while ( ARL < ARL0 & dARL > 1e-6 ) {
    ARL.old <- ARL
    U   <- U + 1
    ARL <- ginar1.arl(L, U, lambda, beta)
    dARL <- abs(ARL - ARL.old)
    if ( OUTPUT ) cat(paste(U, "\t", ARL, "\t", ARL0, "\n"))
  }
  if ( ARL < ARL0 ) U <- NA
  U
}


# Modified chart without randomization for GINAR(1) counts
# (L,U) with minimal L so that ARL = ARL(lambda) is
# increasing at given lambda

ginar1.get.UL <- function(ARL0, lambda, beta, OUTPUT=FALSE) {
  L <- 0 
  U <- ginar1.get.U(L, ARL0, lambda, beta, OUTPUT=OUTPUT) 
  lARL <- ginar1.arl(L, U, lambda-1e-3, beta)
  ARL  <- ginar1.arl(L, U, lambda, beta)
  while ( ARL < lARL ) {
    L    <- L + 1
    U    <- ginar1.get.U(L, ARL0, lambda, beta, OUTPUT=OUTPUT) 
    if ( is.na(U) ) {
      L    <- L - 1
      U    <- ginar1.get.U(L, ARL0, lambda, beta, OUTPUT=OUTPUT) 
      break
    }
    lARL <- ginar1.arl(L, U, lambda-1e-2, beta)
    ARL  <- ginar1.arl(L, U, lambda, beta)
    if ( OUTPUT ) cat(paste("ginar1.get.UL ", L, "\t", U, "\t", lARL, "\t", ARL, "\t\t", ARL0, "\n"))
  }
  c(L, U)
}


# Modified chart with randomization for GINAR(1) counts
# L, U, gL given, gU searched for

ginar1.get.gU <- function(L, U, gL, ARL0, lambda, beta, OUTPUT=FALSE) {
  minARL <- ginar1.arl2(L, U, gL, 1, lambda, beta)
  maxARL <- ginar1.arl2(L, U, gL, 0, lambda, beta)
  if ( OUTPUT ) cat(paste("ginar1.get.gU ", minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gU1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gU1 > 0.1 ) {
      gU2  <- gU1
      ARL2 <- ARL1
      gU1  <- gU1 - .1
      ARL1 <- ginar1.arl2(L, U, gL, gU1, lambda, beta)
    }
    if ( gU1 < .05 ) {
      gU1  <- 0.05
      ARL1 <- ginar1.arl2(L, U, gL, gU1, lambda, beta)
    }
    # secant rule
    gU.error <- 1
    L.error  <- 1
    while ( gU.error > 1e-10  &  L.error > 1e-9 ) {
      gU3  <- gU1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gU2 - gU1)
      gU3 <- max(0, gU3)
      gU3 <- min(1, gU3)
      ARL3 <- ginar1.arl2(L, U, gL, gU3, lambda, beta)
      if ( OUTPUT ) cat(paste("ginar1.get.gU ", gU3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gU1 <- gU2; gU2 <- gU3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gU.error <- abs(gU2 - gU1)
    }
  } else {
    gU3 <- NA                      # old
    if ( minARL > ARL0 ) gU3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gU3 <- -1 # maximal reachable ARL too small
  }
  gU3
}


# Modified chart with randomization for GINAR(1) counts
# L, U, gU given, gL searched for

ginar1.get.gL <- function(L, U, gU, ARL0, lambda, beta, OUTPUT=FALSE) {
  minARL <- ginar1.arl2(L, U, 1, gU, lambda, beta)
  maxARL <- ginar1.arl2(L, U, 0, gU, lambda, beta)
  if ( OUTPUT ) cat(paste("ginar1.get.gL ", minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gL1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gL1 > .1 ) {
      gL2  <- gL1
      ARL2 <- ARL1
      gL1  <- gL1 - .1
      ARL1 <- ginar1.arl2(L, U, gL1, gU, lambda, beta)
    }
    if ( gL1 < 0.05 ) {
      gL1  <- 0.05
      ARL1 <- ginar1.arl2(L, U, gL1, gU, lambda, beta)
    }
    # secant rule
    gL.error <- 1
    L.error  <- 1
    while ( gL.error > 1e-10  &  L.error > 1e-9 ) {
      gL3  <- gL1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gL2 - gL1)
      gL3 <- max(0, gL3)
      gL3 <- min(1, gL3)
      ARL3 <- ginar1.arl2(L, U, gL3, gU, lambda, beta)
      if ( OUTPUT ) cat(paste("ginar1.get.gL ", gL3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gL1 <- gL2; gL2 <- gL3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gL.error <- abs(gL2 - gL1)
    }
  } else {
    gL3 <- NA                      # old
    if ( minARL > ARL0 ) gL3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gL3 <- -1 # maximal reachable ARL too small
  }
  gL3
}


# Some helpers

ginar1.min.gL <- function(L, U, ARL0, lambda, beta, OUTPUT=FALSE) {
  gL  <- 0
  gU  <- ginar1.get.gU(L, U, gL, ARL0, lambda, beta, OUTPUT=OUTPUT)
  if ( -1.5 < gU & gU < 0 ) gL <- NA
  if ( gU < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gU < -1.5 & gL < 1-1e-10 ) {
          gL <- gL + 10^(-dig)
          gU <- ginar1.get.gU(L, U, gL, ARL0, lambda, beta, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("ginar1.min.gL ", "gL =", gL, ",\tgU =", gU, "\n"))
        }
        if ( gU < -1.5 ) {
          gL <- NA
          break
        }
      } else {
        while ( -1.5 < gU  &  gL > 1e-10 ) {
          gL <- gL - 10^(-dig)
          gU <- ginar1.get.gU(L, U, gL, ARL0, lambda, beta, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("ginar1.min.gL ", "gL =", gL, ",\tgU =", gU, "\n"))
        }
      }
    }
  }
  gL
}


ginar1.min.gU <- function(L, U, ARL0, lambda, beta, OUTPUT=FALSE) {
  gU  <- 0
  gL  <- ginar1.get.gL(L, U, gU, ARL0, lambda, beta, OUTPUT=OUTPUT)
  if ( -1.5 < gL & gL < 0 ) gU <- NA
  if ( gL < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gL < -1.5 & gU < 1 ) {
          gU <- gU + 10^(-dig)
          gL <- ginar1.get.gL(L, U, gU, ARL0, lambda, beta, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("ginar1.min.gU ", "gU =", gU, ",\tgL =", gL, "\n"))
        }
        if ( gL < -1.5 ) {
          gU <- NA
          break
        }
      } else {
        while ( -1.5 < gL & gU > 0 ) {
          gU <- gU - 10^(-dig)
          gL <- ginar1.get.gL(L, U, gU, ARL0, lambda, beta, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("ginar1.min.gU ", "gU =", gU, ",\tgL =", gL, "\n"))
        }
      }
    }
  }
  gU
}


# Modified chart with randomization for GINAR(1) counts
# get all 4 constants!

ginar1.get.UL2 <- function(ARL0, lambda, beta, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4) {
  if ( target=="lambda" ) {
    l1 <- lambda - eps
    l2 <- lambda + eps
    b1 <- beta
    b2 <- beta
  }
  if ( target=="beta" ) {
    l1 <- lambda
    l2 <- lambda
    b1 <- beta - eps
    b2 <- beta + eps
  }
  LU  <- ginar1.get.UL(ARL0, lambda, beta, OUTPUT=OUTPUT)
  L2  <- LU[1]
  U2  <- LU[2]
  if ( L2 > 0 ) {
    L1  <- L2 - 1
    U1  <- ginar1.get.U(L1, ARL0, lambda, beta, OUTPUT=OUTPUT)
  } else {
    L1 <- 0
    U1 <- U2
    U2 <- U1 + 1
  }
  #if ( U1==U2 ) U2 <- U2 + 1
  #  U2 <- 100
  #  U2 <- 3*qgeom(1-1/(2*ARL0), lambda)  # A possible solution
  #  U2 <- 2*qgeom(1-1/(2*ARL0), lambda)  # A possible solution
  #  U2 <- search.GandC(lambda, 1/ARL0)$cU
  #  U2 <- 739  #p0=0.01
  U2 <- 145  #p0=0.05
  
  for ( L in L1:L2 ) {
    for ( U in U1:U2 ) {
      gL1  <- ginar1.min.gL(L, U, ARL0, lambda, beta, OUTPUT=OUTPUT)
      if ( is.na(gL1) ) next
      gU1  <- ginar1.get.gU(L, U, gL1, ARL0, lambda, beta)
      lARL1 <- ginar1.arl2(L, U, gL1, gU1, l1, b1)
      rARL1 <- ginar1.arl2(L, U, gL1, gU1, l2, b2)
      dratio1 <- ( rARL1 - lARL1 ) / ( 2*eps )
      
      gU2 <- ginar1.min.gU(L, U, ARL0, lambda, beta, OUTPUT=OUTPUT)
      if ( is.na(gU2) ) next
      gL2 <- ginar1.get.gL(L, U, gU2, ARL0, lambda, beta)
      lARL2 <- ginar1.arl2(L, U, gL2, gU2, l1, b1)
      rARL2 <- ginar1.arl2(L, U, gL2, gU2, l2, b2)
      dratio2 <- ( rARL2 - lARL2 ) / ( 2*eps )
      
      if ( OUTPUT ) cat( paste("ginar1.get.UL2 ", "L =", L, ",\t U =", U, ",\tdr1 =", dratio1, ",\tdr2 =", dratio2, "\n") )
      
      if ( dratio1 * dratio2 < 0 ) {
        # secant rule
        gL.error <- 1
        dr.error <- 1
        while ( gL.error > 1e-10  &  dr.error > delta ) {
          gL3 <- gL1 + (0 - dratio1)/(dratio2 - dratio1)*(gL2 - gL1)
          gU3 <- ginar1.get.gU(L, U, gL3, ARL0, lambda, beta)
          lARL3 <- ginar1.arl2(L, U, gL3, gU3, l1, b1)
          rARL3 <- ginar1.arl2(L, U, gL3, gU3, l2, b2)
          dratio3 <- ( rARL3 - lARL3 ) / ( 2*eps )
          if ( OUTPUT ) cat(paste("ginar1.get.UL2 ", gL3, "\t", dratio3, "\n"))
          gL1 <- gL2; gL2 <- gL3
          dratio1 <- dratio2; dratio2 <- dratio3
          dr.error  <- abs(dratio2)
          gL.error <- abs(gL2 - gL1)
        }
        L0 <- L
        U0 <- U
        L <- L2 + 1
        U <- U2 + 1
      }
      if ( U > U2 ) break
    }
    if ( L > L2 ) break
  }
  data.frame(L=L0, gL=gL3, U=U0, gU=gU3)
}


#########################
####### SHINY APP #######
#########################

ui <- fluidPage(
  # Add custom CSS for tab color change
  tags$head(
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"),  # Load Font Awesome 5
    tags$style(HTML("
    /* Change the background color of the navbar */
    .navbar-nav > li > a {
      background-color: #f8f9fa;  /* Default color for tabs */
      color: #333333;  /* Text color */
    }
    
    /* Change the background and text color when the tab is hovered */
    .navbar-nav > li > a:hover {
      background-color: #e7e7e7;
      color: #000000;
    }
    
    /* Change the background and text color of the active tab */
    .navbar-nav > li.active > a {
      background-color: #007bff !important;  /* Blue background for active tab */
      color: white !important;  /* White text for active tab */
    }
  "))),
  
  navbarPage(
    title = "",
    
    # Home Tab
    tabPanel(
      title = HTML('<i class="fas fa-home" style="margin-right: 5px;"></i> Home'),  # Use "fas" for Font Awesome 5
      fluidPage(
        fluidRow(
          column(12, align = "left",
                 h4(HTML("Welcome to <strong style='color: blue;'>ARL-unbiased Shewhart-type charts</strong> App!")),
                 br(),
                 
                 # Goal section
                 h4(HTML('<i class="fas fa-bullseye" style="margin-right: 10px;"></i> <strong>Goal</strong>')),  # Add Font Awesome goal icon
                 p("The goal of this application is to provide insights into ARL-unbiased Shewhart-type charts to effectively visualize and analyze control charts for i.i.d and autoregressive counts."),
                 br(),
                 
                 # Intended Audience section
                 h4(HTML('<i class="fas fa-users" style="margin-right: 10px;"></i> <strong>Intended Audience</strong>')),
                 p("This app is intended for data scientists, statisticians, and analysts who are interested in Statistical Quality Control. 
                 It can also be useful for Quality Control practitioners and draw their attention to the pivotal role of ARL-unbiased control charts in statistical process control (SPC). "),
                 p("For those interested in more technical details, supplementary documentation and materials can be provided."),
                 br(),
                 
                 # How to Use This App section
                 h4(HTML('<i class="fas fa-info-circle" style="margin-right: 10px;"></i> <strong>How to Use This App</strong>')),
                 p("1. Select the tab ''App'' and navigate through the tabs to explore different features."),
                 p("2. To produce an ARL-unbiased chart, select if you are considering independent or AR(1) counts."),
                 p("3. Then, choose the desired distribution and, making use of the sliders, the control parameters."),
                 p("4. A chart will appear in the main panel, with the control limits and randomization probabilities."),
                 p("5. If you desire, you can download the plot by clicking on ''Download Plot as PDF''.")
          )
        )
      )
    ),
    
    # Main App Tab
    tabPanel(
      title = "App",
      fluidPage(
        titlePanel("ARL-unbiased Shewhart-type charts"),
        
        tabsetPanel(
          # Tab for IID counts
          tabPanel("Independent Counts",
                   sidebarLayout(
                     sidebarPanel(
                       radioButtons("distribution_iid", "Select Distribution:",
                                    choices = c("Poisson", "binomial", "geometric")),
                       
                       # Conditional panels based on the distribution selection
                       conditionalPanel(
                         condition = "input.distribution_iid == 'Poisson'",
                         sliderInput("l0_poisson", HTML("&lambda;<sub>0</sub>"), min = 0, max = 50, value = 19, step = 0.1),
                         sliderInput("arl0_poisson", HTML("ARL<sub>0</sub>"), min = 100, max = 500, value = 370.4, step = 0.1)
                       ),
                       
                       conditionalPanel(
                         condition = "input.distribution_iid == 'binomial'",
                         sliderInput("n_binomial", "n", min = 1, max = 100, value = 22),
                         sliderInput("p0_binomial", HTML("p<sub>0</sub>"), min = 0, max = 1, value = 0.05, step = 0.001),
                         sliderInput("arl0_binomial", HTML("ARL<sub>0</sub>"), min = 100, max = 500, value = 370.4, step = 0.1)
                       ),
                       
                       conditionalPanel(
                         condition = "input.distribution_iid == 'geometric'",
                         sliderInput("p0_geometric", HTML("p<sub>0</sub>"), min = 0, max = 1, value = 0.05, step = 0.001),
                         sliderInput("arl0_geometric", HTML("ARL<sub>0</sub>"), min = 100, max = 500, value = 370.4, step = 0.1)
                       ),
                       
                       # Download button for the plot
                       downloadButton("downloadPlotIID", "Download Plot as PDF")
                     ),
                     mainPanel(
                       withSpinner(plotOutput("arlPlotIID"), type = 8, color = "blue", size = 0.5), # Circular spinner
                       #br(),
                       DTOutput("summaryTableIID"),  # The summary table (without title)
                       br()
                     )
                   )
          ),
          
          # Tab for Autocorrelated counts
          tabPanel("AR(1) Counts",
                   sidebarLayout(
                     sidebarPanel(
                       radioButtons("distribution_ar", "Select Distribution:",
                                    choices = c("Poisson", "binomial", "geometric")),
                       
                       # Conditional panels for autocorrelated counts
                       conditionalPanel(
                         condition = "input.distribution_ar == 'geometric'",
                         sliderInput("p0_geometric_ar", HTML("p<sub>0</sub>"), min = 0, max = 1, value = 0.3, step = 0.001),
                         sliderInput("rho0_geometric_ar", HTML("&rho;<sub>0</sub>"), min = 0, max = 1, value = 0.4, step = 0.01),
                         sliderInput("arl0_geometric_ar", HTML("ARL<sub>0</sub>"), min = 100, max = 500, value = 370.4, step = 0.1),
                         radioButtons("target_geometric", "Shifts in:", choices = c("p", "rho"))
                       ),
                       
                       conditionalPanel(
                         condition = "input.distribution_ar == 'Poisson'",
                         sliderInput("lambda0_poisson_ar", HTML("&lambda;<sub>0</sub>"), min = 0, max = 1.5, value = 0.4636, step = 0.001),
                         sliderInput("beta0_poisson_ar", HTML("&beta;<sub>0</sub>"), min = 0, max = 1, value = 0.81, step = 0.01),
                         sliderInput("arl0_poisson_ar", HTML("ARL<sub>0</sub>"), min = 100, max = 500, value = 500, step = 0.1),
                         radioButtons("target_poisson", "Select Target", choices = c("lambda", "beta"))
                       ),
                       
                       conditionalPanel(
                         condition = "input.distribution_ar == 'binomial'",
                         sliderInput("n_binomial_ar", "n", min = 1, max = 30, value = 15, step = 1),
                         sliderInput("p0_binomial_ar", HTML("p<sub>0</sub>"), min = 0, max = 1, value = 0.38, step = 0.001),
                         sliderInput("rho0_binomial_ar", HTML("&rho;<sub>0</sub>"), min = 0, max = 1, value = 0.97, step = 0.001),
                         sliderInput("arl0_binomial_ar", HTML("ARL<sub>0</sub>"), min = 100, max = 500, value = 370.4, step = 0.1),
                         radioButtons("target_binomial", "Shifts in:", choices = c("p", "rho"))
                       ),
                       
                       # Download button for the plot
                       downloadButton("downloadPlotAR", "Download Plot as PDF")
                     ),
                     mainPanel(
                       withSpinner(plotOutput("arlPlotAR"), type = 8, color = "blue", size = 0.5), # Circular spinner
                       #br(),
                       DTOutput("summaryTableAR"),  # The summary table (without title)
                       br()
                     )
                   )
          )
        )
      )
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  
  ## iid
  
  # Countdown timer for the plot
  countdown_time <- reactiveVal(0)
  
  observeEvent(input$l0_poisson, {
    countdown_time(5) # Assume 5 seconds for the countdown
    invalidateLater(1000, session)
  })
  
  observe({
    if (countdown_time() > 0) {
      isolate({
        countdown_time(countdown_time() - 1)
        invalidateLater(1000, session)
      })
    }
  })
  
  # Plot rendering logic
  output$arlPlotIID <- renderPlot({
    distribution_iid <- input$distribution_iid
    
    if (distribution_iid == "Poisson") {
      l0 <- input$l0_poisson
      arl0 <- input$arl0_poisson
      alpha <- 1 / arl0  # Set alpha as the inverse of arl0
      
      if (l0 == 0) {
        message <- "λ must be different than 0. "
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        results <- search.GandC_c(l0, alpha)
        cL <- results$cL
        cU <- results$cU
        gL <- results$gL
        gU <- results$gU
        
        l0_vals <- seq(1, 50, by = 0.001)
        arl_vals <- arl2c(cL, cU, gL, gU, l0_vals)
        
        plot(l0_vals, arl_vals, type = 'l', col = 'blue', xlab = expression(lambda), ylab = 'ARL')
        abline(v = l0, h = arl0, lty = 4, col = "grey")
        
        #legend("topright", legend = c(paste("LCL =", cL), paste("UCL =", cU), paste("γL =", round(gL, 4)), paste("γU =", round(gU, 4))), cex = 0.8, inset = c(0.025, 0.1))
      }
    } else if (distribution_iid == "binomial") {
      n <- input$n_binomial
      p0 <- input$p0_binomial
      arl0 <- input$arl0_binomial
      alpha <- 1 / arl0
      
      if (p0 == 0 | p0 == 1) {
        message <- "p0 must be different than 0 and 1."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        result <- search.GandC_np(n, p0, alpha)
        cL1 <- result$cL
        cU1 <- result$cU
        gL1 <- result$gL
        gU1 <- result$gU
        
        delta_vals <- seq(-p0, 1 - p0, by = 0.001)
        arl_vals <- sapply(delta_vals, function(delta) {
          binom.unb.arl(n, p0, arl0, delta)
        })
        
        plot(delta_vals + p0, arl_vals, type = 'l', col = 'blue', xlab = 'p', ylab = 'ARL')
        abline(v = p0, h = arl0, lty = 4, col = "grey")
        
        #legend("topright", legend = c(paste("LCL =", cL1), paste("UCL =", cU1), paste("γL =", round(gL1, 4)), paste("γU =", round(gU1, 4))), cex = 0.8, inset = c(0.025, 0.1))
      }
    } else if (distribution_iid == "geometric")
    {
      p0 <- input$p0_geometric
      arl0 <- input$arl0_geometric
      alpha <- 1 / arl0
      
      if (p0 == 0 | p0 == 1) {
        message <- "p0 must be different than 0 and 1."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        results <- search.GandC_g(p0, alpha)
        cL1 <- results$cL
        cU1 <- results$cU
        gL1 <- results$gL
        gU1 <- results$gU
        
        delta_vals <- seq(-p0, 1 - p0, by = 0.001)
        arl_vals <- sapply(delta_vals, function(delta) {
          geom.unb.arl(p0, arl0, p0 + delta)
        })
        
        plot(delta_vals + p0, arl_vals, type = 'l', col = 'blue', xlab = 'p', ylab = 'ARL')
        abline(v = p0, h = arl0, lty = 4, col = "grey")
        
        #legend("topright", legend = c(paste("LCL =", cL1), paste("UCL =", cU1), paste("γL =", round(gL1, 4)), paste("γU =", round(gU1, 4))), cex = 0.8, inset = c(0.025, 0.1))
      }
    }
  })
  
  ## ar(1)
  
  observeEvent(input$distribution_ar, {
    # Update the target radio buttons based on the selected distribution
    if (input$distribution_ar == "Poisson") {
      updateRadioButtons(session, "targetpoisson", choices = c("lambda", "beta"), selected = "lambda")
    } else if (input$distribution_ar == "binomial") {
      updateRadioButtons(session, "targetbinomial", choices = c("p", "rho"), selected = "p")
    } else if (input$distribution_ar == "geometric") {
      updateRadioButtons(session, "targetgeometric", choices = c("p", "rho"), selected = "p")
    }
  })
  
  output$arlPlotAR <- renderPlot({
    distribution_ar <- input$distribution_ar
    
    if (distribution_ar == "Poisson") {
      # Wait for the countdown to finish
      if (countdown_time() > 0) return(NULL)
      
      # Retrieve input values
      lambda0 <- input$lambda0_poisson_ar
      beta0 <- input$beta0_poisson_ar
      arl0 <- input$arl0_poisson_ar
      target <- input$target_poisson
      
      LUg <- tryCatch({
        inar1.get.UL2(arl0, lambda0, beta0, target)
      }, error = function(e) {
        NULL
      })
      
      if (is.null(LUg)) {
        message <- "λ must be different than 0 and β must be in (0,1). "
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        L  <- LUg$L
        gL <- round(LUg$gL, digits=6)
        U  <- LUg$U
        gU <- round(LUg$gU, digits=6)
        
        if (target == "beta") {
          # overall ARL curve for varying beta and fixed lambda
          beta   <- seq(0, 0.99, by=0.001)
          LL <- rep(NA, length(beta))
          for (i in 1:length(beta)) LL[i] <- inar1.arl2(L, U, gL, gU, lambda0, beta[i])
          
          plot(beta, LL, type="l", col = 'blue', xlab=expression(beta), ylab= "ARL")
          abline(h=arl0, v=beta0, lty=4, col="grey")
          
        } else if (target == "lambda") {
          # overall ARL curve for varying lambda and fixed beta
          lambda <- seq(0.001, 1.5, by=0.001)
          LL <- rep(NA, length(lambda))
          for (i in 1:length(lambda)) LL[i] <- inar1.arl2(L, U, gL, gU, lambda[i], beta0)
          
          plot(lambda, LL, type="l",col = 'blue', xlab=expression(lambda), ylab="ARL")
          abline(h=arl0, v=lambda0, lty=4, col="grey")
        }
        
        #legend("topright", legend = c(paste("LCL =", L), paste("UCL =", U), paste("γL =", round(gL, 4), "  "), paste("γU =", round(gU, 4))), cex = 0.8, inset = c(0.025, 0.1))
      }
    } else if (distribution_ar == "binomial") {
      # Wait for the countdown to finish
      if (countdown_time() > 0) return(NULL)
      
      # Retrieve input values
      n <- input$n_binomial_ar
      p0 <- input$p0_binomial_ar
      rho0 <- input$rho0_binomial_ar
      arl0 <- input$arl0_binomial_ar
      target <- input$target_binomial
      
      rho_range <- calculate_rho_range(input$p0_binomial_ar)
      
      LUg <- tryCatch({
        ar1.get.UL2(ARL0 = arl0, n, p0, rho0, target="p")
      }, error = function(e) {
        NULL
      })
      
      if (is.null(LUg)) {
        message <- "Unable to plot for this combination of values."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
        return()
      }
      
      L  <- LUg$L
      gL <- round(LUg$gL, digits=6)
      U  <- LUg$U
      gU <- round(LUg$gU, digits=6)
      
      if (target == "p") {
        # overall ARL curve for varying p and fixed rho
        p   <- seq(0.01, 0.99, by=0.001)
        LL <- rep(NA, length(p))
        for (i in 1:length(p)) LL[i] <- ar1.arl2(L, U, gL, gU, n, p[i], rho0)
        
        plot(p, LL, type="l", col = 'blue', xlab=expression(p), ylab= "ARL")
        abline(h=arl0, v=p0, lty=4, col="grey")
        
      } else if (target == "rho") {
        # overall ARL curve for varying rho and fixed p
        rho_seq <- seq(rho_range[1], 0.99, by=0.001)
        LL <- rep(NA, length(rho_seq))
        for (i in 1:length(rho_seq)) LL[i] <- ar1.arl2(L, U, gL, gU, n, p0, rho_seq[i])
        
        plot(rho_seq, LL, type="l", col = 'blue', xlab=expression(rho), ylab="ARL")
        abline(h=arl0, v=rho0, lty=4, col="grey")
      }
      
      #legend("topleft", legend = c(paste("LCL =", L), paste("UCL =", U), paste("γL =", round(gL, 4), "  "), paste("γU =", round(gU, 4))), cex = 0.8, inset = c(0.025, 0.05))
      
    } else if (distribution_ar == "geometric") {
      # Wait for the countdown to finish
      if (countdown_time() > 0) return(NULL)
      
      # Retrieve input values
      p0 <- input$p0_geometric_ar
      rho0 <- input$rho0_geometric_ar
      arl0 <- input$arl0_geometric_ar
      targetgeometric <- input$target_geometric
      
      LUg <- tryCatch({
        ginar1.get.UL2(ARL0 = arl0, lambda = p0, beta = rho0, target = "lambda")
      }, error = function(e) {
        NULL
      })
      
      if (is.null(LUg)) {
        message <- "Unable to plot for this combination of values."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
        return()
      }
      
      L  <- LUg$L
      gL <- round(LUg$gL, digits=6)
      U  <- LUg$U
      gU <- round(LUg$gU, digits=6)
      
      if (targetgeometric == "p") {
        # overall ARL curve for varying p and fixed rho
        p   <- seq(0.01, 0.99, by=0.001)
        LL <- rep(NA, length(p))
        for (i in 1:length(p)) LL[i] <- ginar1.arl2(L, U, gL, gU, p[i], rho0)
        
        plot(p, LL, type="l", col = 'blue', xlab=expression(p), ylab= "ARL")
        abline(h=arl0, v=p0, lty=4, col="grey")
        
      } else if (targetgeometric == "rho") {
        # overall ARL curve for varying rho and fixed p
        rho <- seq(0.01, 0.99, by=0.001)
        LL <- rep(NA, length(rho))
        for (i in 1:length(rho)) LL[i] <- ginar1.arl2(L, U, gL, gU, p0, rho[i])
        
        plot(rho, LL, type="l",col = 'blue', xlab=expression(rho), ylab="ARL")
        abline(h=arl0, v=rho0, lty=4, col="grey")
      }
      
      #legend("topleft", legend = c(paste("LCL =", L), paste("UCL =", U), paste("γL =", round(gL, 4), "  "), paste("γU =", round(gU, 4))), cex = 0.8, inset = c(0.025, 0.1))
    }
  })
  
  ## download iid
  
  # Download handler for the PDF plot
  output$downloadPlotIID <- downloadHandler(
    filename = function() {
      paste("arl_plot_iid_", input$distribution_iid, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 7)  # Set appropriate PDF dimensions
      par(mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins
      
      # Reproduce the plot based on the selected distribution
      distribution_iid <- input$distribution_iid
      
      if (distribution_iid == "Poisson") {
        l0 <- input$l0_poisson
        arl0 <- input$arl0_poisson
        alpha <- 1 / arl0  # Set alpha as the inverse of arl0
        
        if (l0 == 0) {
          message <- "λ must be different than 0. "
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, message, cex = 1.2)
        } else {
          results <- search.GandC_c(l0, alpha)
          cL <- results$cL
          cU <- results$cU
          gL <- results$gL
          gU <- results$gU
          
          l0_vals <- seq(1, 50, by = 0.001)
          arl_vals <- arl2c(cL, cU, gL, gU, l0_vals)
          
          plot(l0_vals, arl_vals, type = 'l', col = 'blue', xlab = expression(lambda), ylab = 'ARL')
          abline(v = l0, h = arl0, lty = 4, col = "grey")
          
          #legend("topright", legend = c(paste("LCL =", cL), paste("UCL =", cU), bquote(gamma[L] == .(round(gL, 4))), bquote(gamma[U] == .(round(gU, 4)))), cex = 0.8, inset = c(0.025, 0.1))
        }
        
      } else if (distribution_iid == "binomial") {
        n <- input$n_binomial
        p0 <- input$p0_binomial
        arl0 <- input$arl0_binomial
        alpha <- 1 / arl0
        
        if (p0 == 0 | p0 == 1) {
          message <- "p0 must be different than 0 and 1."
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, message, cex = 1.2)
        } else {
          result <- search.GandC_np(n, p0, alpha)
          cL1 <- result$cL
          cU1 <- result$cU
          gL1 <- result$gL
          gU1 <- result$gU
          
          delta_vals <- seq(-p0, 1 - p0, by = 0.001)
          arl_vals <- sapply(delta_vals, function(delta) {
            binom.unb.arl(n, p0, arl0, delta)
          })
          
          plot(delta_vals + p0, arl_vals, type = 'l', col = 'blue', xlab = 'p', ylab = 'ARL')
          abline(v = p0, h = arl0, lty = 4, col = "grey")
          
          #legend("topright", legend = c(paste("LCL =", cL1), paste("UCL =", cU1, "     "), bquote(gamma[L] == .(round(gL1, 4))), bquote(gamma[U] == .(round(gU1, 4)))), cex = 0.8, inset = c(0.025, 0.1))
        }
        
      } else if (distribution_iid == "geometric") {
        p0 <- input$p0_geometric
        arl0 <- input$arl0_geometric
        alpha <- 1 / arl0
        
        if (p0 == 0 | p0 == 1) {
          message <- "p0 must be different than 0 and 1."
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, message, cex = 1.2)
        } else {
          results <- search.GandC_g(p0, alpha)
          cL1 <- results$cL
          cU1 <- results$cU
          gL1 <- results$gL
          gU1 <- results$gU
          
          delta_vals <- seq(-p0, 1 - p0, by = 0.001)
          arl_vals <- sapply(delta_vals, function(delta) {
            geom.unb.arl(p0, arl0, p0 + delta)
          })
          
          plot(delta_vals + p0, arl_vals, type = 'l', col = 'blue', xlab = 'p', ylab = 'ARL')
          abline(v = p0, h = arl0, lty = 4, col = "grey")
          
          #legend("topright", legend = c(paste("LCL =", cL1), paste("UCL =", cU1), bquote(gamma[L] == .(round(gL1, 4))), bquote(gamma[U] == .(round(gU1, 4)))), cex = 0.8, inset = c(0.025, 0.1))
        }
      }
      dev.off()  # Close the PDF device
    })
  
  ## download ar(1)
  
  output$downloadPlotAR <- downloadHandler(
    filename = function() {
      paste("arl_plot_ar_", input$distribution_ar, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 7)  # Set appropriate PDF dimensions
      par(mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins
      
      # Reproduce the plot based on the selected distribution
      distribution_ar <- input$distribution_ar
      if (distribution_ar == "Poisson") {
        lambda0 <- input$lambda0_poisson_ar
        beta0 <- input$beta0_poisson_ar
        arl0 <- input$arl0_poisson_ar
        target <- input$target_poisson
        
        LUg <- tryCatch({
          inar1.get.UL2(arl0, lambda0, beta0, target)
        }, error = function(e) {
          NULL
        })
        
        if (is.null(LUg)) {
          message <- "λ must be different than 0 and β must be in (0,1). "
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, message, cex = 1.2)
        } else {
          L  <- LUg$L
          gL <- round(LUg$gL, digits=6)
          U  <- LUg$U
          gU <- round(LUg$gU, digits=6)
          
          if (target == "beta") {
            # overall ARL curve for varying beta and fixed lambda
            beta   <- seq(0, 0.99, by=0.001)
            LL <- rep(NA, length(beta))
            for (i in 1:length(beta)) LL[i] <- inar1.arl2(L, U, gL, gU, lambda0, beta[i])
            
            plot(beta, LL, type="l", col = 'blue', xlab=expression(beta), ylab= "ARL")
            abline(h=arl0, v=beta0, lty=4, col="grey")
            
          } else if (target == "lambda") {
            # overall ARL curve for varying lambda and fixed beta
            lambda <- seq(0.001, 1.5, by=0.001)
            LL <- rep(NA, length(lambda))
            for (i in 1:length(lambda)) LL[i] <- inar1.arl2(L, U, gL, gU, lambda[i], beta0)
            
            plot(lambda, LL, type="l",col = 'blue', xlab=expression(lambda), ylab="ARL")
            abline(h=arl0, v=lambda0, lty=4, col="grey")
          }
          
          #legend("topright", legend = c(paste("LCL =", L), paste("UCL =", U), bquote(gamma[L] == .(round(gL, 4))), bquote(gamma[U] == .(round(gU, 4)))), cex = 0.8, inset = c(0.025, 0.1))
        }
      } else if (distribution_ar == "binomial") {
        
        n <- input$n_binomial_ar
        p0 <- input$p0_binomial_ar
        rho0 <- input$rho0_binomial_ar
        arl0 <- input$arl0_binomial_ar
        target <- input$target_binomial
        
        rho_range <- calculate_rho_range(input$p0_binomial_ar)
        
        LUg <- tryCatch({
          ar1.get.UL2(ARL0 = arl0, n, p0, rho0, target="p")
        }, error = function(e) {
          NULL
        })
        
        if (is.null(LUg)) {
          message <- "Unable to plot for this combination of values."
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, message, cex = 1.2)
          return()
        }
        
        L  <- LUg$L
        gL <- round(LUg$gL, digits=6)
        U  <- LUg$U
        gU <- round(LUg$gU, digits=6)
        
        if (target == "p") {
          # overall ARL curve for varying p and fixed rho
          p   <- seq(0.01, 0.99, by=0.001)
          LL <- rep(NA, length(p))
          for (i in 1:length(p)) LL[i] <- ar1.arl2(L, U, gL, gU, n, p[i], rho0)
          
          plot(p, LL, type="l", col = 'blue', xlab=expression(p), ylab= "ARL")
          abline(h=arl0, v=p0, lty=4, col="grey")
          
        } else if (target == "rho") {
          # overall ARL curve for varying rho and fixed p
          rho_seq <- seq(rho_range[1], 0.99, by=0.001)
          LL <- rep(NA, length(rho_seq))
          for (i in 1:length(rho_seq)) LL[i] <- ar1.arl2(L, U, gL, gU, n, p0, rho_seq[i])
          
          plot(rho_seq, LL, type="l", col = 'blue', xlab=expression(rho), ylab="ARL")
          abline(h=arl0, v=rho0, lty=4, col="grey")
        }
        
        #legend("topleft", legend = c(paste("LCL =", L), paste("UCL =", U), bquote(gamma[L] == .(round(gL, 4))), bquote(gamma[U] == .(round(gU, 4)))), cex = 0.8, inset = c(0.025, 0.05))
        
      } else if (distribution_ar == "geometric") {
        
        p0 <- input$p0_geometric_ar
        rho0 <- input$rho0_geometric_ar
        arl0 <- input$arl0_geometric_ar
        targetgeometric <- input$target_geometric
        
        LUg <- tryCatch({
          ginar1.get.UL2(ARL0 = arl0, lambda = p0, beta = rho0, target = "lambda")
        }, error = function(e) {
          NULL
        })
        
        if (is.null(LUg)) {
          message <- "Unable to plot for this combination of values."
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, message, cex = 1.2)
          return()
        }
        
        L  <- LUg$L
        gL <- round(LUg$gL, digits=6)
        U  <- LUg$U
        gU <- round(LUg$gU, digits=6)
        
        if (targetgeometric == "p") {
          # overall ARL curve for varying p and fixed rho
          p   <- seq(0.01, 0.99, by=0.001)
          LL <- rep(NA, length(p))
          for (i in 1:length(p)) LL[i] <- ginar1.arl2(L, U, gL, gU, p[i], rho0)
          
          plot(p, LL, type="l", col = 'blue', xlab=expression(p), ylab= "ARL")
          abline(h=arl0, v=p0, lty=4, col="grey")
          
        } else if (targetgeometric == "rho") {
          # overall ARL curve for varying rho and fixed p
          rho <- seq(0.01, 0.99, by=0.001)
          LL <- rep(NA, length(rho))
          for (i in 1:length(rho)) LL[i] <- ginar1.arl2(L, U, gL, gU, p0, rho[i])
          
          plot(rho, LL, type="l",col = 'blue', xlab=expression(rho), ylab="ARL")
          abline(h=arl0, v=rho0, lty=4, col="grey")
        }
        
        #legend("topleft", legend = c(paste("LCL =", L), paste("UCL =", U), bquote(gamma[L] == .(round(gL, 4))), bquote(gamma[U] == .(round(gU, 4)))), cex = 0.8, inset = c(0.025, 0.1))
      }
      dev.off()  # Close the PDF device
    })
  
  ## table iid 
  
  output$summaryTableIID <- renderDT({
    
    distribution_iid <- input$distribution_iid
    if (distribution_iid == "Poisson") {
      l0 <- input$l0_poisson
      arl0 <- input$arl0_poisson
      alpha <- 1 / arl0  # Set alpha as the inverse of arl0
      
      if (l0 == 0) {
        message <- "λ must be different than 0. "
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        results <- search.GandC_c(l0, alpha)
        cL <- results$cL
        cU <- results$cU
        gL <- results$gL
        gU <- results$gU
        
        # Create a summary table in vertical format
        #summary_data <- data.frame(
        #  Parameter = c(htmltools::HTML("&lambda;<sub>0</sub>"), htmltools::HTML("ARL<sub>0</sub>"), "LCL", "UCL", htmltools::HTML("&gamma;<sub>L</sub>"), htmltools::HTML("&gamma;<sub>U</sub>")),
        #  Value = c(l0, arl0, cL, cU, round(gL, 4), round(gU,4))
        #)
        
        #datatable(summary_data, 
        #          rownames = FALSE,  # Remove row numbers
        #          escape = FALSE,    # Allow HTML tags to be rendered
        #          options = list(dom = 't',  # Only show the table
        #                         paging = FALSE,  # No pagination needed for small table
        #                         ordering = FALSE,  # Disable sorting
        #                         searching = FALSE))  # Disable search
        
        
        # Create a summary table in horizontal format
        summary_data <- data.frame(
          Lambda0 = l0,
          ARL0 = arl0,
          LCL = cL,
          UCL = cU,
          GammaL = round(gL, 4),
          GammaU = round(gU, 4)
        )
        
        datatable(summary_data, 
                  rownames = FALSE,  # Remove row numbers
                  escape = FALSE,    # Allow HTML tags (for subscript)
                  options = list(
                    dom = 't',       # Only show the table
                    paging = FALSE,  # No pagination needed for small table
                    ordering = FALSE,  # Disable sorting
                    searching = FALSE,  # Disable search
                    columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center all columns
                  ),
                  # Apply HTML for the column names (for subscripts)
                  colnames = c(
                    htmltools::HTML("&lambda;<sub>0</sub>"),
                    htmltools::HTML("ARL<sub>0</sub>"),
                    "LCL",
                    "UCL",
                    htmltools::HTML("&gamma;<sub>L</sub>"),
                    htmltools::HTML("&gamma;<sub>U</sub>")
                  )
        ) %>%
          formatStyle(columns = names(summary_data), textAlign = 'center')  # Center the text in all columns
      }
      
    } else if (distribution_iid == "binomial") {
      n <- input$n_binomial
      p0 <- input$p0_binomial
      arl0 <- input$arl0_binomial
      alpha <- 1 / arl0
      
      if (p0 == 0 | p0 == 1) {
        message <- "p0 must be different than 0 and 1."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        result <- search.GandC_np(n, p0, alpha)
        cL1 <- result$cL
        cU1 <- result$cU
        gL1 <- result$gL
        gU1 <- result$gU
        
        # Create a summary table in vertical format
        #summary_data <- data.frame(
        #  Parameter = c("n", htmltools::HTML("p<sub>0</sub>"), htmltools::HTML("ARL<sub>0</sub>"), "LCL", "UCL", htmltools::HTML("&gamma;<sub>L</sub>"), htmltools::HTML("&gamma;<sub>U</sub>")),
        #  Value = c(n, p0, arl0, cL1, cU1, round(gL1, 4), round(gU1,4))
        #)
        
        #datatable(summary_data, 
        #          rownames = FALSE,  # Remove row numbers
        #          escape = FALSE,    # Allow HTML tags to be rendered
        #          options = list(dom = 't',  # Only show the table
        #                         paging = FALSE,  # No pagination needed for small table
        #                         ordering = FALSE,  # Disable sorting
        #                         searching = FALSE))  # Disable search
        
        
        # Create a summary table in horizontal format
        summary_data <- data.frame(
          n = n,
          p0 = p0,
          ARL0 = arl0,
          LCL = cL1,
          UCL = cU1,
          GammaL = round(gL1, 4),
          GammaU = round(gU1, 4)
        )
        
        datatable(summary_data, 
                  rownames = FALSE,  # Remove row numbers
                  escape = FALSE,    # Allow HTML tags (for subscript)
                  options = list(
                    dom = 't',       # Only show the table
                    paging = FALSE,  # No pagination needed for small table
                    ordering = FALSE,  # Disable sorting
                    searching = FALSE,  # Disable search
                    columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center all columns
                  ),
                  # Apply HTML for the column names (for subscripts)
                  colnames = c(
                    "n",
                    htmltools::HTML("p<sub>0</sub>"),
                    htmltools::HTML("ARL<sub>0</sub>"),
                    "LCL",
                    "UCL",
                    htmltools::HTML("&gamma;<sub>L</sub>"),
                    htmltools::HTML("&gamma;<sub>U</sub>")
                  )
        ) %>%
          formatStyle(columns = names(summary_data), textAlign = 'center')  # Center the text in all columns
        
      }
      
    } else if (distribution_iid == "geometric") {
      p0 <- input$p0_geometric
      arl0 <- input$arl0_geometric
      alpha <- 1 / arl0
      
      if (p0 == 0 | p0 == 1) {
        message <- "p0 must be different than 0 and 1."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        results <- search.GandC_g(p0, alpha)
        cL1 <- results$cL
        cU1 <- results$cU
        gL1 <- results$gL
        gU1 <- results$gU
        
        # Create a summary table in vertical format
        #summary_data <- data.frame(
        #  Parameter = c(htmltools::HTML("p<sub>0</sub>"), htmltools::HTML("ARL<sub>0</sub>"), "LCL", "UCL", htmltools::HTML("&gamma;<sub>L</sub>"), htmltools::HTML("&gamma;<sub>U</sub>")),
        #  Value = c(p0, arl0, cL1, cU1, round(gL1,4), round(gU1,4))
        #)
        
        #datatable(summary_data, 
        #          rownames = FALSE,  # Remove row numbers
        #          escape = FALSE,    # Allow HTML tags to be rendered
        #          options = list(dom = 't',  # Only show the table
        #                         paging = FALSE,  # No pagination needed for small table
        #                         ordering = FALSE,  # Disable sorting
        #                         searching = FALSE))  # Disable search
        
        
        # Create a summary table in horizontal format
        summary_data <- data.frame(
          p0 = p0,
          ARL0 = arl0,
          LCL = cL1,
          UCL = cU1,
          GammaL = round(gL1, 4),
          GammaU = round(gU1, 4)
        )
        
        datatable(summary_data, 
                  rownames = FALSE,  # Remove row numbers
                  escape = FALSE,    # Allow HTML tags (for subscript)
                  options = list(
                    dom = 't',       # Only show the table
                    paging = FALSE,  # No pagination needed for small table
                    ordering = FALSE,  # Disable sorting
                    searching = FALSE,  # Disable search
                    columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center all columns
                  ),
                  # Apply HTML for the column names (for subscripts)
                  colnames = c(
                    htmltools::HTML("p<sub>0</sub>"),
                    htmltools::HTML("ARL<sub>0</sub>"),
                    "LCL",
                    "UCL",
                    htmltools::HTML("&gamma;<sub>L</sub>"),
                    htmltools::HTML("&gamma;<sub>U</sub>")
                  )
        ) %>%
          formatStyle(columns = names(summary_data), textAlign = 'center')  # Center the text in all columns
        
      }
    }
    
  })
  
  ## table ar(1)
  
  output$summaryTableAR <- renderDT({
    
    distribution_ar <- input$distribution_ar
    
    if (distribution_ar == "Poisson") {
      
      lambda0 <- input$lambda0_poisson_ar
      beta0 <- input$beta0_poisson_ar
      arl0 <- input$arl0_poisson_ar
      target <- input$target_poisson
      
      LUg <- tryCatch({
        inar1.get.UL2(arl0, lambda0, beta0, target)
      }, error = function(e) {
        NULL
      })
      
      if (is.null(LUg)) {
        message <- "λ must be different than 0 and β must be in (0,1). "
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
        L  <- LUg$L
        gL <- round(LUg$gL, digits=6)
        U  <- LUg$U
        gU <- round(LUg$gU, digits=6)
        
        # Create a summary table in horizontal format
        summary_data <- data.frame(
          Lambda0 = lambda0,
          beta0 = beta0,
          ARL0 = arl0,
          LCL = L,
          UCL = U,
          GammaL = round(gL, 4),
          GammaU = round(gU, 4)
        )
        
        datatable(summary_data, 
                  rownames = FALSE,  # Remove row numbers
                  escape = FALSE,    # Allow HTML tags (for subscript)
                  options = list(
                    dom = 't',       # Only show the table
                    paging = FALSE,  # No pagination needed for small table
                    ordering = FALSE,  # Disable sorting
                    searching = FALSE,  # Disable search
                    columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center all columns
                  ),
                  # Apply HTML for the column names (for subscripts)
                  colnames = c(
                    htmltools::HTML("&lambda;<sub>0</sub>"),
                    htmltools::HTML("&beta;<sub>0</sub>"),
                    htmltools::HTML("ARL<sub>0</sub>"),
                    "LCL",
                    "UCL",
                    htmltools::HTML("&gamma;<sub>L</sub>"),
                    htmltools::HTML("&gamma;<sub>U</sub>")
                  )
        ) %>%
          formatStyle(columns = names(summary_data), textAlign = 'center')  # Center the text in all columns
      }
      
    } else if (distribution_ar == "binomial") {
      # Retrieve input values
      n <- input$n_binomial_ar
      p0 <- input$p0_binomial_ar
      rho0 <- input$rho0_binomial_ar
      arl0 <- input$arl0_binomial_ar
      target <- input$target_binomial
    
      LUg <- tryCatch({
        ar1.get.UL2(ARL0 = arl0, n, p0, rho0, target="p")
      }, error = function(e) {
        NULL
      })
      
      if (is.null(LUg)) {
        message <- "Unable to plot for this combination of values."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
      
      L  <- LUg$L
      gL <- round(LUg$gL, digits=6)
      U  <- LUg$U
      gU <- round(LUg$gU, digits=6)
    
      
      # Create a summary table in horizontal format
      summary_data <- data.frame(
        n = n,
        p0 = p0,
        rho0 = rho0,
        ARL0 = arl0,
        LCL = L,
        UCL = U,
        GammaL = round(gL, 4),
        GammaU = round(gU, 4)
      )
      
      datatable(summary_data, 
                rownames = FALSE,  # Remove row numbers
                escape = FALSE,    # Allow HTML tags (for subscript)
                options = list(
                  dom = 't',       # Only show the table
                  paging = FALSE,  # No pagination needed for small table
                  ordering = FALSE,  # Disable sorting
                  searching = FALSE,  # Disable search
                  columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center all columns
                ),
                # Apply HTML for the column names (for subscripts)
                colnames = c(
                  "n",
                  htmltools::HTML("p<sub>0</sub>"),
                  htmltools::HTML("&rho;<sub>0</sub>"),
                  htmltools::HTML("ARL<sub>0</sub>"),
                  "LCL",
                  "UCL",
                  htmltools::HTML("&gamma;<sub>L</sub>"),
                  htmltools::HTML("&gamma;<sub>U</sub>")
                )
      ) %>%
        formatStyle(columns = names(summary_data), textAlign = 'center')  # Center the text in all columns
      }
      
    } else if (distribution_ar == "geometric") {
      # Retrieve input values
      p0 <- input$p0_geometric_ar
      rho0 <- input$rho0_geometric_ar
      arl0 <- input$arl0_geometric_ar
      target <- input$target_geometric
      
      LUg <- tryCatch({
        ginar1.get.UL2(ARL0 = arl0, lambda = p0, beta = rho0, target = "lambda")
      }, error = function(e) {
        NULL
      })
      
      if (is.null(LUg)) {
        message <- "Unable to plot for this combination of values."
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, message, cex = 1.2)
      } else {
      
      L  <- LUg$L
      gL <- round(LUg$gL, digits=6)
      U  <- LUg$U
      gU <- round(LUg$gU, digits=6)
      
      # Create a summary table in horizontal format
      summary_data <- data.frame(
        p0 = p0,
        rho0 = rho0,
        ARL0 = arl0,
        LCL = L,
        UCL = U,
        GammaL = round(gL, 4),
        GammaU = round(gU, 4)
      )
      
      datatable(summary_data, 
                rownames = FALSE,  # Remove row numbers
                escape = FALSE,    # Allow HTML tags (for subscript)
                options = list(
                  dom = 't',       # Only show the table
                  paging = FALSE,  # No pagination needed for small table
                  ordering = FALSE,  # Disable sorting
                  searching = FALSE,  # Disable search
                  columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center all columns
                ),
                # Apply HTML for the column names (for subscripts)
                colnames = c(
                  htmltools::HTML("p<sub>0</sub>"),
                  htmltools::HTML("&rho;<sub>0</sub>"),
                  htmltools::HTML("ARL<sub>0</sub>"),
                  "LCL",
                  "UCL",
                  htmltools::HTML("&gamma;<sub>L</sub>"),
                  htmltools::HTML("&gamma;<sub>U</sub>")
                )
      ) %>%
        formatStyle(columns = names(summary_data), textAlign = 'center')  # Center the text in all columns
      }
    }
  })
  
}

# Run the App
shinyApp(ui, server)








