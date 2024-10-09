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
  if (rho<max(-p/(1-p),-(1-p)/p) || rho>1) break
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

binom.ar1.unb.arl <- function(n, p0, rho0, arl0, p, rho){
  LUg <- ar1.get.UL2(arl0, n, p0, rho0, target="p", OUTPUT=FALSE, eps=1e-6, delta=1e-4)
  L  <- LUg$L
  gL <- round(LUg$gL, digits=6)
  U  <- LUg$U
  gU <- round(LUg$gU, digits=6)
  arl <- round(ar1.arl2(L, U, gL, gU, n, p, rho), digits=4)
  return(arl)
}

# Define the range and steps for p and rho
p_seq <- seq(0.01, 0.5, length.out = 30)  # Example range for p
rho_seq <- seq(0.01, 0.5, length.out = 30)    # Example range for rho

# Create a matrix to store L values
L_matrix <- matrix(0, nrow = length(p_seq), ncol = length(rho_seq))
U_matrix <- matrix(0, nrow = length(p_seq), ncol = length(rho_seq))
gL_matrix <- matrix(0, nrow = length(p_seq), ncol = length(rho_seq))
gU_matrix <- matrix(0, nrow = length(p_seq), ncol = length(rho_seq))


# Fill the L_matrix using nested loops
for (i in 1:length(p_seq)) {
  print(i)
  for (j in 1:length(rho_seq)) {
    n <- 15
    ARL0 <- 370.4
    p <- p_seq[i]
    rho <- rho_seq[j]
    
    # Call the function and extract L, U, gL and gU
    result <- ar1.get.UL2(ARL0, n, p, rho, target="p", OUTPUT=FALSE, eps=1e-6, delta=1e-4)
    L_matrix[i, j] <- result$L
    U_matrix[i, j] <- result$U
    gL_matrix[i, j] <- result$gL
    gU_matrix[i, j] <- result$gU
  }
}

# Write the matrices to a text file
write.table(L_matrix, file = "L_matrix_output.txt", row.names = FALSE, col.names = TRUE)
write.table(U_matrix, file = "U_matrix_output.txt", row.names = FALSE, col.names = TRUE)
write.table(gL_matrix, file = "gL_matrix_output.txt", row.names = FALSE, col.names = TRUE)
write.table(gU_matrix, file = "gU_matrix_output.txt", row.names = FALSE, col.names = TRUE)

## já tenho os valores guardados em .txt

L_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/2 Binomial (Done)/3D Plot/L_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
U_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/2 Binomial (Done)/3D Plot/U_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
gL_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/2 Binomial (Done)/3D Plot/gL_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
gU_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/2 Binomial (Done)/3D Plot/gU_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

L_matrix <- as.matrix(L_matrix)
U_matrix <- as.matrix(U_matrix)
gL_matrix <- as.matrix(gL_matrix)
gU_matrix <- as.matrix(gU_matrix)


# Create the 3D plot using persp
persp(p_seq, rho_seq, L_matrix, theta = -30, phi = 15, expand = 0.5, col = "lightgray",
      xlab = "p0", ylab = "ρ0", zlab = "LCL",
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0))
persp(p_seq, rho_seq, U_matrix, theta = -30, phi = 15, expand = 0.5, col = "lightgray",
      xlab = "p0", ylab = "ρ0", zlab = "UCL",
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0))
persp(p_seq, rho_seq, gL_matrix, theta = -30, phi = 10, expand = 0.5, col = "lightgray",
      xlab = "p0", ylab = "ρ0", zlab = "ɣL",
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0))
persp(p_seq, rho_seq, gU_matrix, theta = -30, phi = 10, expand = 0.5, col = "lightgray",
      xlab = "p0", ylab = "ρ0", zlab = "ɣU",
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0))

##########################################################

# ARL com p design

n <- 15
p0 <- 0.3
# rho0 <- varia

p_seq <- seq(0.1, 0.9, by=0.01)
rho_seq <- seq(0.1, 0.9, by=0.01)

# Create a matrix to store ARL values
arl_matrix <- matrix(0, nrow = length(rho_seq), ncol = length(p_seq))

# Fill the L_matrix using nested loops
for (i in 1:length(rho_seq)) {
  print(i)
  for (j in 1:length(p_seq)) {
    
    ARL0 <- 370.4
    
    arl_result <- binom.ar1.unb.arl(n, p0, rho0 = rho_seq[i], arl0 = ARL0, p = p_seq[j], rho = rho_seq[i])
    arl_matrix[i, j] <- arl_result
  }
}

write.table(arl_matrix, file = "arl_matrix_output.txt", row.names = FALSE, col.names = TRUE)

arl_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/2 Binomial (Done)/3D Plot/arl_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
arl_matrix <- as.matrix(arl_matrix)

# Create the 3D plot using persp
persp(rho_seq, p_seq, arl_matrix, theta = 70, phi = 10, expand = 0.4, col = "lightgray",
      ylab = "p", xlab = "ρ0", zlab = "ARL",
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0))

# with:
# ARL0 <- 370.4
# p0 <- 0.3
# n <- 50
# rho0 <- varia

#### CHATGPT 

# Your existing code to generate the 3D plot
persp_result <- persp(rho_seq, p_seq, arl_matrix, theta = 50, phi = 10, expand = 0.4, col = "lightgray",
                      ylab = "p", xlab = "ρ0", zlab = "ARL",
                      ticktype = "detailed", # Use detailed ticks to avoid overlay
                      cex.axis = 0.5, # Adjust axis label size if needed
                      cex.lab = 1,  # Adjust axis label size
                      mgp = c(2, 0.5, 0))

# Define specific values for x and y where the line should be drawn
p_0 <- 0.3
rho_0 <- 0.9

# Define z coordinates for the vertical line (from min to max ARL values)
z_coords <- c(min(arl_matrix), max(arl_matrix))

# Transform the x, y, and z coordinates to the 2D perspective of the plot
transformed <- trans3d(
  x = rep(rho_0, 2),       # x = 0.9 for both points
  y = rep(p_0, 2),     # y = 0.5 for both points
  z = z_coords,             # from min to max of ARL values
  pmat = persp_result       # transformation matrix from persp
)

# Draw the vertical line using the transformed coordinates
lines(transformed, col = "red", lwd = 2)







