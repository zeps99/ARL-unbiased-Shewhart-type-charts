# Defining the value of the ARL unbiased g-chart

# function F etc.

# Function P, F and etc
P <- function(x, p0) {(1-p0)^x * p0} 

F <- function(x, p0) {1 - (1-p0)^(x+1)}

iF1 <- function(f, p0) {
  if (f == p0) {
    stop("error")
  } else {
    ceiling(log(1-f)/log(1-p0)) -1
  }
} 

iF2 <- function(f, p0) {
  if (f == p0) {
    stop("error")
  } else {
    ceiling(log(1-f)/log(1-p0)) -1
  }
}


# Function G and etc
G_ <- function(x, p0) {1 - (1-p0)^x * (1+p0*x)}

G <- Vectorize("G_", "x")

iG1_ <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( G(x, p0) < g  ) {
        while ( G(x, p0) < g ) x <- x + 1
      } else {
        while ( G(x, p0) >= g & x > 0 ) x <- x - 1
        if ( G(x, p0) < g ) x <- x + 1
      }
    }
  }
  x
}
iG1 <- Vectorize("iG1_", "g")

iG2_ <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( G(x, p0) <= g  ) {
        while ( G(x, p0) <= g ) x <- x + 1
      } else {
        while ( G(x, p0) > g & x > 0 ) x <- x - 1
        if ( G(x, p0) <= g ) x <- x + 1
      }
    }
  }
  x
}
iG2 <- Vectorize("iG2_", "g")

# Bounds

c.bounds <- function(alpha, p0) {
  Lmax1 <- iF2(alpha, p0)
  Lmax2 <- iG2(alpha, p0)
  Umin1 <- iF1(1-alpha, p0)
  Umin2 <- iG1(1-alpha, p0)
  Lmax <- min(Lmax1, Lmax2)
  Umin <- max(Umin1, Umin2)
  
  Lmin1 <- iF1(max(F(Umin-1,p0)-1+alpha,0), p0)
  Lmin2 <- iG1(max(G(Umin-1,p0)-1+alpha,0), p0)
  Umax1 <- iF2(min(F(Lmax,p0)+1-alpha,1), p0)
  Umax2 <- iG2(min(G(Lmax,p0)+1-alpha,1), p0)
  Lmin <- max(Lmin1, Lmin2)
  Umax <- min(Umax1, Umax2)
  
  data.frame(p0, Lmin, Lmax, Umin, Umax)
}


# Randomization Constants
gammas <- function(cL, cU, p0, alpha) {
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
search.GandC <- function(p0, alpha, OUTPUT=FALSE) {
  CB <- c.bounds(alpha, p0)
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
    
    GG <- gammas(cL, cU, p0, alpha)
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

#' Compute ARLs of unbiased g-charts
#' 
#' @description
#' Computation of the ARL unbiased g-chart at a given p0, ARL0 and p
#'
#' @param p0 target value of the fraction non-conforming
#' @param arl0 average run length in control or supremum of the ARL function
#' @param p fraction nonconforming
#'
#' @return it returns single value which resembles the ARL for a Geometric distribution
#'
#' @export
#'
#' @author José Guedes, Manuel C. Morais
#'
#' @references Morais, M.C. (2017). ARL-unbiased geometric and CCCG control charts. Sequential Analysis 36, 513-527. http://dx.doi.org/10.1080/07474946.2017.1394717
#'
#' @examples ## In-control process, ρ = 1
#' arl0 <- 370.4
#' p0 <-  0.001
#' p <- 1*p0
#' geom.unb.arl(p0,arl0,p)
#'
geom.unb.arl <- function(p0,arl0, p){
  results <- search.GandC(p0, (arl0)^(-1))
  cL4  <- results$cL
  cU4  <- results$cU
  gL4  <- results$gL
  gU4  <- results$gU

  # p <- p0*rho
  aL <- (1 - (1-l)^(cL4)) + gL4*((1-l)^cL4 * l) 
  aU <- (1-l)^(cU4+1) + gU4*((1-l)^cU4 * l)
  ARL  <- 1 / ( aL + aU )
  ARL
}
