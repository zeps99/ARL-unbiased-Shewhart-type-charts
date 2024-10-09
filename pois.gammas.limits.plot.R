
# Paulino, Morais & Knoth (2016a) - Algorithm

# F functions

F <- function(x, l0) ppois(x, l0)

iF1 <- function(f, l0) qpois(f, l0)

iF2_ <- function(f, l0) {
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
}
iF2 <- Vectorize("iF2_", "f")


# G functions

G_ <- function(x, l0) sum( dpois(0:x, l0) * (0:x) )/l0

G <- Vectorize("G_", "x")

iG1_ <- function(g, l0) {
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
}
iG1 <- Vectorize("iG1_", "g")

iG2_ <- function(g, l0) {
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
}
iG2 <- Vectorize("iG2_", "g")


# control limits

c.bounds <- function(alpha, l0) {
  Lmax1 <- iF2(alpha, l0)
  Lmax2 <- iG2(alpha, l0)  
  Umin1 <- iF1(1-alpha, l0)
  Umin2 <- iG1(1-alpha, l0)  
  Lmax <- min(Lmax1, Lmax2)
  Umin <- max(Umin1, Umin2)
  
  Lmin1 <- iF1(max(F(Umin-1,l0)-1+alpha,0), l0)
  Lmin2 <- iG1(max(G(Umin-1,l0)-1+alpha,0), l0)
  Umax1 <- iF2(min(F(Lmax,l0)+1-alpha,1), l0)
  Umax2 <- iG2(min(G(Lmax,l0)+1-alpha,1), l0)
  Lmin <- max(Lmin1, Lmin2)
  Umax <- min(Umax1, Umax2)
  
  data.frame(l0, Lmin, Lmax, Umin, Umax)
}


# Randomization probabilities

gammas <- function(cL, cU, l0, alpha) {
  gU <- ( l0*( alpha - G(cL-1,l0) - 1+G(cU,l0) ) - cL*( alpha - F(cL-1,l0) - 1+F(cU,l0) ) ) / (dpois(cU,l0) * (cU-cL))
  gL <- ( cU*( alpha - F(cL-1,l0) - 1+F(cU,l0) ) - l0*( alpha - G(cL-1,l0) - 1+G(cU,l0) ) ) / (dpois(cL,l0) * (cU-cL))
  gg <- data.frame(gL, gU)
  return(gg)
}  


# Search procedure to obtain the control limits and the randomization probabilities

search.GandC <- function(l0, alpha, OUTPUT=FALSE) {
  CB <- c.bounds(alpha, l0)  
  results <- NULL
  for ( cL in CB$Lmin:CB$Lmax ) {
    for ( cU in CB$Umin:CB$Umax ) {
      GG <- gammas(cL, cU, l0, alpha)
      good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
      if ( OUTPUT ) cat(paste("cL =", cL, ", cU =", cU, ", gL =", GG[1], ", gU =", GG[2], "\n"))
      if ( good ) {
        results <-data.frame(l0, CB$Lmin, CB$Lmax, cL, CB$Umin, CB$Umax, cU, gL=GG[1], gU=GG[2])
        break
      }
    }
  } 
  return(results)
}


ARL0 <- 370.4
lambda <- seq(0.05, 20, by=0.01)
ccc <- NULL

for ( i in 1:length(lambda) ) {
  cat(paste(i, "\t", lambda[i], "\n"))
  ccc <- rbind(ccc, search.GandC(l0 = lambda[i], alpha = 1/ARL0, OUTPUT=FALSE))
}

dL  <- c(1, which( diff(ccc$cL) > 0 ) + 1) #da os indices onde ? verdadeiro, ou seja, onde h? altera??o do L 
gL  <- l2 <- NULL
for ( i in 1:(length(dL)-1) ) {
  gL <- c(gL, ccc$gL[dL[i]:(dL[i+1]-1)], NA) #da os gL e mete NA entre as altera??es de L
  l2 <- c(l2, lambda[dL[i]:(dL[i+1]-1)], NA) #da os p e mete NA entre as altera??es de L
}
gL <- c(gL, ccc$gL[max(dL):nrow(ccc)]) #junta os ultimos gl depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 
l2 <- c(l2, lambda[max(dL):nrow(ccc)]) #junta os ultimos p depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 

dU  <- c(1, which( diff(ccc$cU) > 0 ) + 1)
gU  <- l3 <- NULL
for ( i in 1:(length(dU)-1) ) {
  gU <- c(gU, ccc$gU[dU[i]:(dU[i+1]-1)], NA)
  l3 <- c(l3, lambda[dU[i]:(dU[i+1]-1)], NA)
}
gU <- c(gU, ccc$gU[max(dU):nrow(ccc)])
l3 <- c(l3, lambda[max(dU):nrow(ccc)])

par(mar=c(4,4,1,4))
plot(lambda, ccc$cL, ylim=c(0, max(ccc$cL)*1.5), xlab=expression(lambda[0]), ylab="LCL", cex=.1)
par(new=T)
plot(l2, gL, type="l", col="black", axes=F, bty="n", xlab="", ylab="")
abline(h=0:1, lty=4, col="grey")
abline(v=0, lty=4, col="grey")
abline(v=20, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[L]), 4, 3)
graphics.off()

par(mar=c(4,4,1,4))
plot(lambda, ccc$cU, ylim=c(0, max(ccc$cU)*1.5), xlab=expression(lambda[0]), ylab="UCL", cex=.1)
par(new=T)
plot(l3, gU, type="l", col="black", axes=F, bty="n", xlab="", ylab="")
abline(h=0:1, lty=4, col="grey")
abline(v=0, lty=4, col="grey")
abline(v=20, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[U]), 4, 3)
graphics.off()
