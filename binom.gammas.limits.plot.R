# Auxiliary functions (ARL-unbiased np-chart)
# - search procedure described by Morais (2016)

F <- function(x,n,p0) pbinom(x,n,p0)

P <- function(y,n,p0) dbinom(y,n,p0)

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

gammas <- function(cL,cU,n,p0,alpha) {
  a<-P(cL,n,p0)
  b<-P(cU,n,p0)
  c<-(cL*P(cL,n,p0))
  d<-(cU*P(cU,n,p0))
  e<-(alpha-1+F(cU,n,p0)-F(cL-1,n,p0))
  f<-(n*p0*(alpha-1+F(cU-1,n-1,p0)-F(cL-2,n-1,p0)))
  gL<-((d*e-b*f)/(a*d-b*c))
  gU<-((a*f-c*e)/(a*d-b*c))
  gg <- data.frame(gL, gU)
  gg
}

emitirsinal <- function(n,p0,cU,cL,gL,gU) {
  prob <- (1-F(cU,n,p0)+F(cL-1,n,p0)+gL*P(cL,n,p0)+gU*P(cU,n,p0))
  arl <- (1/prob)
  emissao <- data.frame(prob,arl)
  emissao
}

search.GandC <- function(n,p0,alpha) {
  results <- NULL
  cL <- initialcL(alpha,n,p0)
  q <- (1-alpha+F(cL-1,n,p0))
  cU <- qbinom(q,n,p0)
  GG <- gammas(cL,cU,n,p0,alpha)
  good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
  cL <- (cL+1)
  while (cL>0) {
    cL <- (cL-1)
    q <- (1-alpha+F(cL-1,n,p0))
    cU <- qbinom(q,n,p0)
    GG <- gammas(cL,cU,n,p0,alpha)
    good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
    if ( good ) {
      sinal <- emitirsinal(n,p0,cU,cL,GG[1],GG[2])
      results <- data.frame(n,p0, cL, cU, gL=GG[1], gU=GG[2],prob=sinal[1],ARL=sinal[2])
      break
    } else {
      cU<-(cU+1)
      GG <- gammas(cL,cU,n,p0,alpha)
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


ARL0 <- 370.4
prob <- seq(0.001, 0.5, by=0.0005)
ccc <- NULL


for ( i in 1:length(prob) ) {
  cat(paste(i, "\t", prob[i], "\n"))
  ccc <- rbind(ccc, search.GandC(n = 100, p0 = prob[i], alpha = 1/ARL0))
}

dL  <- c(1, which( diff(ccc$cL) > 0 ) + 1) #da os indices onde ? verdadeiro, ou seja, onde h? altera??o do L 
gL  <- l2 <- NULL
for ( i in 1:(length(dL)-1) ) {
  gL <- c(gL, ccc$gL[dL[i]:(dL[i+1]-1)], NA) #da os gL e mete NA entre as altera??es de L
  l2 <- c(l2, prob[dL[i]:(dL[i+1]-1)], NA) #da os p e mete NA entre as altera??es de L
}
gL <- c(gL, ccc$gL[max(dL):nrow(ccc)]) #junta os ultimos gl depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 
l2 <- c(l2, prob[max(dL):nrow(ccc)]) #junta os ultimos p depois da ultima altera??o do L, uma vez que n?o s?o adicionados anteriormente 

dU  <- c(1, which( diff(ccc$cU) > 0 ) + 1)
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
plot(l2, gL, type="l", col="black", axes=F, bty="n", xlab="", ylab="")
abline(h=0:1, lty=4, col="grey")
abline(v=0, lty=4, col="grey")
abline(v=0.5, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[L]), 4, 3)
graphics.off()

par(mar=c(4,4,1,4))
plot(prob, ccc$cU, ylim=c(0, max(ccc$cU)*1.5), xlab=expression(p[0]), ylab="UCL", cex=.1)
par(new=T)
plot(l3, gU, type="l", col="black", axes=F, bty="n", xlab="", ylab="")
abline(h=0:1, lty=4, col="grey")
abline(v=0, lty=4, col="grey")
abline(v=0.5, lty=4, col="grey")
axis(4, at=pretty(c(0,1)))
mtext(expression(gamma[U]), 4, 3)
graphics.off()
