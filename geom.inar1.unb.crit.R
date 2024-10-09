####
# ARL of a modified g-chart without randomization for GINAR(1) counts
# lambda and beta below are our p and rho, respectively
# Auxiliary functions
####

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

##########

geom.inar1.unb.crit <- function(p, rho, arl0){
  LUg <- ginar1.get.UL2(ARL0 = arl0, lambda = p, beta = rho, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4)
  L  <- LUg$L
  gL <- round(LUg$gL, digits=6)
  U  <- LUg$U
  gU <- round(LUg$gU, digits=6)
  result <- paste("LCL:", L, " ", "UCL:", U," ", "gL:", round(gL,4), " ", "gU:", round(gU,4))
  return(result)
}

# Example Table 1, MoraisWittenbergKnoth2023
arl0 <- 200
p0 <- 0.63
rho0 <- 0
geom.inar1.unb.crit(p0, rho0, arl0)



####### PARA TABELA TESE #######

# Define the geom.inar1.unb.crit function
geom.inar1.unb.crit <- function(p, rho, arl0){
  LUg <- ginar1.get.UL2(ARL0 = arl0, lambda = p, beta = rho, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4)
  L  <- LUg$L
  gL <- round(LUg$gL, digits=6)
  U  <- LUg$U
  gU <- round(LUg$gU, digits=6)
  result <- paste("(p0, rho0):", p, rho, "(LCL,UCL):", L, U, "(gL, gU):", round(gL,6), round(gU,6))
  return(result)
}

# Define the values for p0
p0_values <- c(0.1, 0.25, 0.5)

# Define the values for rho0
rho0_values <- c(0, 0.25, 0.5)

# Set the arl0 value
arl0 <- 370.4

# Initialize an empty list to store results
results <- list()

# Loop over all combinations of lambda0, beta0, and target

for (p in p0_values) {
  for (rho in rho0_values) {
    # Call the function for the current combination
    result <- geom.inar1.unb.crit(p, rho, arl0)
    # Append the result to the list
    results <- c(results, result)
  }
}


# Print all results
for (res in results) {
  print(res)
}

