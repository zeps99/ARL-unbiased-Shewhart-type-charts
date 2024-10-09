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

# ARL calculation

emitirsinal <- function(n,p0,cU,cL,gL,gU) {
  prob <- (1-F(cU,n,p0)+F(cL-1,n,p0)+gL*P(cL,n,p0)+gU*P(cU,n,p0))
  arl <- (1/prob)
  emissao <- data.frame(prob,arl)
  emissao
}


# arl function

arl2 <- function(n,p0,cU,cL,gL,gU) {
  prob <- (1-F(cU,n,p0)+F(cL-1,n,p0)+gL*P(cL,n,p0)+gU*P(cU,n,p0))
  arl <- (1/prob)
  return(arl)
}


#' Compute ARLs of unbiased np-charts
#'
#' @description
#' Computation of the Average Run Length (ARL) at a given fraction nonconforming p
#'
#'
#' @param n sample size
#' @param p0 target value of the fraction non-conforming
#' @param arl0 average run length in control or supremum of the ARL function
#' @param p fraction nonconforming
#'
#' @return it returns single value which resembles the ARL for a Binomial distribution
#'
#' @export
#'
#' @author José Guedes, Manuel C. Morais
#'
#' @references Morais, M.C. (2016). An ARL-unbiased np-chart. Economic Quality Control 31, 11-21. http://www.degruyter.com/view/j/eqc.2016.31.issue-1/eqc-2015-0013/eqc-2015-0013.xml
#' @references Morais, M.C., Wittenberg, P., and Cruz, C.J. (2022). The np-chart with 3-sigma limits and the ARL-unbiased np-chart revisited. Stochastics and Quality Control 37, 107-116. https://doi.org/10.1515/eqc-2022-0032
#'
#' @examples ## Out-of-control process, M. C. Morais et al., 2022, page 6
#' n <- 400
#' p0 <- 0.03
#' arl0 <- 370.4
#' p <- 1.1*p0
#' binom.unb.arl(n,p0, arl0, p)
#'
binom.unb.arl <- function(n,p0,arl0, p) {
  alpha <- 1/arl0

  result <- search.GandC(n,p0,alpha)

  cL1 <- result$cL
  cU1 <- result$cU
  gL1 <- result$gL
  gU1 <- result$gU

  # p = p0 + delta
  ARL <- arl2(n,p,cU1,cL1,gL1,gU1)
  return(ARL)
}

n <- 400
p0 <- 0.03
arl0 <- 370.4
p <- 1.1*p0
binom.unb.arl(n,p0, arl0, p)

## correndo um a um dá certo! quando meto no package deixa de dar
## e não entendo o porquê

