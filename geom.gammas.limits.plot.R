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

#search.GandC(p0 = 0.00001, alpha = 1/ARL0, OUTPUT=FALSE)

##### for small p0s ####

#ARL0 <- 370.4
prob <- seq(0.0001, 0.005, by=0.000005)
ccc <- NULL


for ( i in 1:length(prob) ) {
  cat(paste(i, "\t", prob[i], "\n"))
  ccc <- rbind(ccc, search.GandC(prob[i], alpha = 1/ARL0, OUTPUT=FALSE))
}

dL  <- c(1, which( abs(diff(ccc$cL)) > 0 ) + 1) #da os indices onde ? verdadeiro, ou seja, onde h? altera??o do L 
gL  <- l2 <- NULL
for ( i in 1:(length(dL)-1) ) {
  gL <- c(gL, ccc$gL[dL[i]:(dL[i+1]-1)], NA) #da os gL e mete NA entre as altera??es de L
  l2 <- c(l2, prob[dL[i]:(dL[i+1]-1)], NA) #da os p e mete NA entre as altera??es de L
}
gL <- c(gL, ccc$gL[max(dL):nrow(ccc)]) #junta os ultimos gl depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 
l2 <- c(l2, prob[max(dL):nrow(ccc)]) #junta os ultimos p depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 

dU  <- c(1, which( abs(diff(ccc$cU)) > 0 ) + 1)
gU  <- l3 <- NULL
for ( i in 1:(length(dU)-1) ) {
  gU <- c(gU, ccc$gU[dU[i]:(dU[i+1]-1)], NA)
  l3 <- c(l3, prob[dU[i]:(dU[i+1]-1)], NA)
}
gU <- c(gU, ccc$gU[max(dU):nrow(ccc)])
l3 <- c(l3, prob[max(dU):nrow(ccc)])

par(mar=c(4,4,1,4))
plot(prob, ccc$cL, ylim=c(0, max(ccc$cL)*1.5), xlab=expression(p[0]), ylab="LCL", cex=.1)
par(new=T)
plot(l2, gL, type="l", col="black", axes=FALSE, bty="n", xlab="", ylab="")
abline(h=0:1, lty=4, col="grey")
abline(v=0.0001, lty=4, col="grey")
abline(v=0.005, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[L]), 4, 2)
graphics.off()

par(mar=c(4,4,1,4))
plot(prob, ccc$cU, ylim=c(0, max(ccc$cU)*1.5), xlab=expression(p[0]), ylab="UCL", cex=.1)
par(new=T)
plot(l3, gU, type="p", col="black", axes=FALSE, bty="n", xlab="", ylab="", pch=16, cex=0.2)
abline(h=0:1, lty=4, col="grey")
abline(v=0.0001, lty=4, col="grey")
abline(v=0.005, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[U]), 4, 2)
graphics.off()

## SÓ PARA CHECKAR ##
install.packages("zoom")
library("zoom")
zm()

####################################
####  NEW - expanding small p0  ####
####################################

install.packages("plotrix")
library(plotrix)
# Scatter plot with zoom
par(mar=c(4,4,1,2))
zoomInPlot(x = l3, y = gU, xlim= c(0,0.005), ylim = c(0, 1), 
           pch=16, cex=0.3, col = "black",
           xlab=expression(p[0]), ylab=expression(gamma[U]),
           rxlim = c(0.0042, 0.0048), 
           rylim = c(0.5, 0.9))


##### for "normal" p0s ####

ARL0 <- 370.4
prob <- seq(0.01, 0.5, by=0.0001)
ccc <- NULL


for ( i in 1:length(prob) ) {
  cat(paste(i, "\t", prob[i], "\n"))
  ccc <- rbind(ccc, search.GandC(prob[i], alpha = 1/ARL0, OUTPUT=FALSE))
}

dL  <- c(1, which( abs(diff(ccc$cL)) > 0 ) + 1) #da os indices onde ? verdadeiro, ou seja, onde h? altera??o do L 
gL  <- l2 <- NULL
for ( i in 1:(length(dL)-1) ) {
  gL <- c(gL, ccc$gL[dL[i]:(dL[i+1]-1)], NA) #da os gL e mete NA entre as altera??es de L
  l2 <- c(l2, prob[dL[i]:(dL[i+1]-1)], NA) #da os p e mete NA entre as altera??es de L
}
gL <- c(gL, ccc$gL[max(dL):nrow(ccc)]) #junta os ultimos gl depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 
l2 <- c(l2, prob[max(dL):nrow(ccc)]) #junta os ultimos p depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 

dU  <- c(1, which( abs(diff(ccc$cU)) > 0 ) + 1)
gU  <- l3 <- NULL
for ( i in 1:(length(dU)-1) ) {
  gU <- c(gU, ccc$gU[dU[i]:(dU[i+1]-1)], NA)
  l3 <- c(l3, prob[dU[i]:(dU[i+1]-1)], NA)
}
gU <- c(gU, ccc$gU[max(dU):nrow(ccc)])
l3 <- c(l3, prob[max(dU):nrow(ccc)])

par(mar=c(4,4,1,4))
plot(prob, ccc$cL, ylim=c(0, 1), xlab=expression(p[0]), ylab="LCL", cex=.1)
abline(h = 1, lty = 4, col="grey")
par(new=T)
plot(l2, gL, type="l", col="black", axes=FALSE, bty="n", xlab="", ylab="")
abline(h=0, lty=4, col="grey")
abline(v=0.01, lty=4, col="grey")
abline(v=0.5, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[L]), 4, 2)
graphics.off()

par(mar=c(4,4,1,4))
plot(prob, ccc$cU, ylim=c(0, max(ccc$cU)*1.5), xlab=expression(p[0]), ylab="UCL", cex=.1)
par(new=T)
plot(l3, gU, type="l", col="black", axes=F, bty="n", xlab="", ylab="")
abline(h=0:1, lty=4, col="grey")
abline(v=0.1, lty=4, col="grey")
abline(v=0.5, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[U]), 4, 2)
graphics.off()

