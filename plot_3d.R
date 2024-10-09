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

pois.inar1.unb.arl <- function(lambda0, beta0, arl0, lambda, beta, target){
  LUg <- inar1.get.UL2(arl0, lambda0, beta0, target)
  L  <- LUg$L
  gL <- round(LUg$gL, digits=6)
  U  <- LUg$U
  gU <- round(LUg$gU, digits=6)
  arl <- round(inar1.arl2(L, U, gL, gU, lambda, beta), digits=4)
  return(arl)
}

### JA ESTA CERTO !!!!

#lambda0 <- 0.5
#beta0 <- 0
#arl0 <- 370.4
#beta <- beta0
#lambda <- 1.2 * lambda0
#pois.inar1.unb.arl(lambda0, beta0, arl0, lambda, beta, target = "lambda")

##########################################################

# LIMITS and GAMMAS

# Define the range and steps for lambda and beta
lambda_seq <- seq(0.5, 10, length.out = 30)  # Example range for lambda
beta_seq <- seq(0.1, 0.5, length.out = 30)    # Example range for beta

# Create a matrix to store L values
L_matrix <- matrix(0, nrow = length(lambda_seq), ncol = length(beta_seq))
U_matrix <- matrix(0, nrow = length(lambda_seq), ncol = length(beta_seq))
gL_matrix <- matrix(0, nrow = length(lambda_seq), ncol = length(beta_seq))
gU_matrix <- matrix(0, nrow = length(lambda_seq), ncol = length(beta_seq))

# Fill the L_matrix using nested loops
for (i in 1:length(lambda_seq)) {
  print(i)
  for (j in 1:length(beta_seq)) {
    ARL0 <- 370.4
    lambda0 <- lambda_seq[i]
    beta0 <- beta_seq[j]
    
    # Call the function and extract L, U, gL and gU
    result <- inar1.get.UL2(ARL0, lambda0, beta0, target="lambda", OUTPUT=FALSE)
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

L_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/1 Poisson (Done)/3D Plot/L_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
U_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/1 Poisson (Done)/3D Plot/U_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
gL_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/1 Poisson (Done)/3D Plot/gL_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
gU_matrix <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/1 Poisson (Done)/3D Plot/gU_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

L_matrix <- as.matrix(L_matrix)
U_matrix <- as.matrix(U_matrix)
gL_matrix <- as.matrix(gL_matrix)
gU_matrix <- as.matrix(gU_matrix)

# Create the 3D plot using persp
persp(lambda_seq, beta_seq, L_matrix, theta = 30, phi = 10, expand = 0.5, col = "lightgray",
      xlab = "λ0", ylab = "β0", zlab = "LCL", ticktype = "detailed",
      cex.axis = 0.5, cex.lab = 1, mgp = c(2, 0.5, 0))
persp(lambda_seq, beta_seq, U_matrix, theta = -30, phi = 10, expand = 0.5, col = "lightgray",
      xlab = "λ0", ylab = "β0", zlab = "UCL", ticktype = "detailed",
      cex.axis = 0.5, cex.lab = 1, mgp = c(2, 0.5, 0))
persp(lambda_seq, beta_seq, gL_matrix, theta = 20, phi = 10, expand = 0.5, col = "lightgray",
      xlab = "λ0", ylab = "β0", zlab = "ɣL", ticktype = "detailed",
      cex.axis = 0.5, cex.lab = 1, mgp = c(2, 0.5, 0))
persp(lambda_seq, beta_seq, gU_matrix, theta = -30, phi = 10, expand = 0.5, col = "lightgray",
      xlab = "λ0", ylab = "β0", zlab = "ɣU", ticktype = "detailed",
      cex.axis = 0.5, cex.lab = 1, mgp = c(2, 0.5, 0))

##########################################################

# ARL com lambda design

lambda0 <- 0.5
# beta0 <- varia

lambda_seq <- seq(0.1, 1.5, by=0.025)
beta_seq <- seq(0.1, 0.9, by=0.01)

# Create a matrix to store L values
arl_matrix <- matrix(0, nrow = length(beta_seq), ncol = length(lambda_seq))

# Fill the L_matrix using nested loops
for (i in 1:length(beta_seq)) {
  print(i)
  for (j in 1:length(lambda_seq)) { 

    ARL0 <- 370.4
    #lambda_seq <- lambda_seq[i]
    #beta_seq <- beta_seq[j]
    
    arl_result <- pois.inar1.unb.arl(lambda0, beta_seq[i], arl0 = ARL0, lambda_seq[j], beta_seq[i], target="lambda")
    #print(arl_result)
    arl_matrix[i, j] <- arl_result
  }
}

write.table(arl_matrix, file = "arl_matrix_output.txt", row.names = FALSE, col.names = TRUE)

data_arl <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/1 Poisson (Done)/3D Plot/arl_matrix_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
data_arl <- as.matrix(data_arl)

# Create the 3D plot using persp
persp(beta_seq, lambda_seq, arl_matrix, theta = 60, phi = 10, expand = 0.4, col = "lightgray",
      ylab = "λ", xlab = "β0", zlab = "ARL",
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0))

# with:
# ARL0 <- 370.4
# lambda0 <- 0.5
# beta0 <- varia

#### CHATGPT 

# Your existing code to generate the 3D plot
persp_result <- persp(
  beta_seq, lambda_seq, arl_matrix, theta = 50, phi = 10, expand = 0.4, col = "lightgray",
  ylab = "λ", xlab = "β0", zlab = "ARL",
  ticktype = "detailed",
  cex.axis = 0.5,
  cex.lab = 1,
  mgp = c(2, 0.5, 0)
)

# Define specific values for x and y where the line should be drawn
beta_0 <- 0.9
lambda_0 <- 0.5

# Define z coordinates for the vertical line (from min to max ARL values)
z_coords <- c(min(arl_matrix), max(arl_matrix))

# Transform the x, y, and z coordinates to the 2D perspective of the plot
transformed <- trans3d(
  x = rep(beta_0, 2),       # x = 0.9 for both points
  y = rep(lambda_0, 2),     # y = 0.5 for both points
  z = z_coords,             # from min to max of ARL values
  pmat = persp_result       # transformation matrix from persp
)

# Draw the vertical line using the transformed coordinates
lines(transformed, col = "red", lwd = 2)

#_#_#_#_#_#_#_#_ com plano em vez de linha

# Your existing code to generate the 3D plot
persp_result <- persp(
  beta_seq, lambda_seq, arl_matrix, theta = 120, phi = 10, expand = 0.4, col = "lightgray",
  ylab = "λ", xlab = "β0", zlab = "ARL",
  ticktype = "detailed",
  cex.axis = 0.5,
  cex.lab = 1,
  mgp = c(2, 0.5, 0)
)

# Define the specific y value for the plane (lambda)
lambda_0 <- 0.5

# Define the range for x (beta) and z (ARL) for the plane
beta_range <- range(beta_seq)  # full range of beta
arl_range <- range(arl_matrix)  # full range of ARL (z-axis)

# Create coordinates for the four corners of the plane
x_coords <- c(beta_range[1], beta_range[2], beta_range[2], beta_range[1])  # x-coordinates (beta)
y_coords <- rep(lambda_0, 4)  # fixed y-coordinate (lambda_0)
z_coords <- c(arl_range[1], arl_range[1], arl_range[2], arl_range[2])  # z-coordinates (ARL range)

# Transform the 3D coordinates to the 2D plotting perspective
transformed_plane <- trans3d(
  x = x_coords,
  y = y_coords,
  z = z_coords,
  pmat = persp_result  # transformation matrix from persp
)

# Draw the plane using the polygon function with transformed coordinates
polygon(transformed_plane, col = rgb(1, 0, 0, alpha = 0.3), border = "red")  # semi-transparent red plane



##########################################################

# ARL com beta design

#lambda0 <- varia
beta0 <- 0.5

lambda_seq <- seq(0.1, 1.5, by=0.05)
beta_seq <- seq(0.1, 0.9, by=0.01)

# Create a matrix to store L values
arl_matrix_beta <- matrix(0, nrow = length(lambda_seq), ncol = length(beta_seq))

# Fill the L_matrix using nested loops
for (i in 1:length(lambda_seq)) {
  print(i)
  for (j in 1:length(beta_seq)) { 
    
    ARL0 <- 370.4
    #lambda_seq <- lambda_seq[i]
    #beta_seq <- beta_seq[j]
    
    arl_result <- pois.inar1.unb.arl(lambda_seq[i], beta0, arl0 = ARL0, lambda_seq[i], beta_seq[j], target="beta")
    #print(arl_result)
    arl_matrix_beta[i, j] <- arl_result
  }
}

write.table(arl_matrix_beta, file = "arl_matrix_beta_output.txt", row.names = FALSE, col.names = TRUE)

data_beta_arl <- read.table("C:/Users/j-gue/OneDrive/Área de Trabalho/2º Ano - Mestrado Matemática Aplicada e Computação/2º Semestre/Dissertação/Codes R/ARL Unbiased/3 AR(1) Shewhart (Done)/1 Poisson (Done)/3D Plot/arl_matrix_beta_output.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
data_beta_arl <- as.matrix(data_beta_arl)

# Create the 3D plot using persp
persp(lambda_seq, beta_seq, data_beta_arl, theta = 50, phi = 10, expand = 0.4, col = "lightgray",
      xlab = "λ0", ylab = "β", zlab = "ARL",
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0))

# with:
# ARL0 <- 370.4
# lambda0 <- varia
# beta0 <- 0.5

#### CHATGPT 

# Your existing code to generate the 3D plot
persp_result <- persp(lambda_seq, beta_seq, data_beta_arl, theta = 140, phi = 20, expand = 0.4, col = "lightgray",
                      xlab = "λ0", ylab = "β", zlab = "ARL",
                      ticktype = "detailed", # Use detailed ticks to avoid overlay
                      cex.axis = 0.5, # Adjust axis label size if needed
                      cex.lab = 1,  # Adjust axis label size
                      mgp = c(2, 0.5, 0))

# Define specific values for x and y where the line should be drawn
beta_0 <- 0.5
lambda_0 <- 1.5

# Define z coordinates for the vertical line (from min to max ARL values)
z_coords <- c(min(arl_matrix_beta), max(arl_matrix_beta))

# Transform the x, y, and z coordinates to the 2D perspective of the plot
transformed <- trans3d(
  x = rep(lambda_0, 2),   
  y = rep(beta_0, 2),       
  z = z_coords,             
  pmat = persp_result       
)

# Draw the vertical line using the transformed coordinates
lines(transformed, col = "red", lwd = 2)



########################################
#### NAO INTERESSA DAQUI PARA BAIXO ####
########################################

# CALCULAR VALOR ARL0 para dif. valores de lambda0 e beta0

# Define the range and steps for lambda and beta
lambda_seq <- seq(0.01, 1.5, by=0.05)
beta_seq <- seq(0.01, 0.99, by=0.05)

# Create a matrix to store L values
arl0_matrix <- matrix(0, nrow = length(lambda_seq), ncol = length(beta_seq))

# Fill the L_matrix using nested loops
for (i in 1:length(lambda_seq)) {
  for (j in 1:length(beta_seq)) {
    print(j)
    ARL0 <- 370.4
    
    LUg <- inar1.get.UL2(ARL0, lambda_seq[i], beta_seq[j], target = "lambda")
    L  <- LUg$L
    gL <- round(LUg$gL, digits=6)
    U  <- LUg$U
    gU <- round(LUg$gU, digits=6)
    
    arl_result <- inar1.arl2(L, U, gL, gU, lambda_seq[i], beta_seq[j])
    arl0_matrix[i, j] <- arl_result
  }
}

# Create the 3D plot using persp
persp(lambda_seq, beta_seq, arl0_matrix, theta = 30, phi = 20, expand = 0.4, col = "lightgray",
      xlab = "λ0", ylab = "β0", zlab = "ARL")


persp(lambda_seq, beta_seq, arl0_matrix, theta = 30, phi = 10, expand = 0.4, col = "lightgray",
      xlab = "λ0", ylab = "β0", zlab = "ARL", 
      ticktype = "detailed", # Use detailed ticks to avoid overlay
      cex.axis = 0.5, # Adjust axis label size if needed
      cex.lab = 1,  # Adjust axis label size
      mgp = c(2, 0.5, 0)) # Adjusts margin lines for the axis titles


# Write the matrix to a text file
write.table(arl0_matrix, file = "arl0_matrix_output.txt", row.names = FALSE, col.names = TRUE)


