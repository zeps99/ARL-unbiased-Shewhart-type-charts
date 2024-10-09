# Defining the control limits and randomization probabilities of ARL unbiased np-chart

# Auxiliary functions

F <- function(x,n,p0) pbinom(x,n,p0)
P <- function(y,n,p0) dbinom(y,n,p0)

initialcL_ <- function(alpha,n,p0) {
  if ( alpha < 0 | alpha > 1 ) {
    stop("error")
  } else {
    if ( alpha == 1 ) {
      x <- Inf
    } else {
      inv<-qbinom(alpha,n,p0)
      if (F(inv,n,p0)>alpha) {
        x <- (inv)
      } else {
        x <- (inv+1)
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
  return(gg)
}

# ARL calculation
emitirsinal <- function(n,p0,cU,cL,gL,gU) {
  prob <- (1-F(cU,n,p0)+F(cL-1,n,p0)+gL*P(cL,n,p0)+gU*P(cU,n,p0))
  arl <- (1/prob)
  emissao <- data.frame(prob,arl)
  emissao
}

# Search procedure to obtain the control limits and the randomization probabilities
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
  return(results)
}

#' Randomization probabilities and control limits of the ARL unbiased np-chart
#'
#' @description
#' Computation of the randomization probabilities and control limits of the Average Run Length (ARL) unbiased np-chart at a given n, p0, ARL0
#'
#' @param n sample size
#' @param p0 target value of the fraction nonconforming
#' @param arl0 average run length in control or supremum of the ARL function
#'
#' @return it returns four values, respectively, the lower and upper control limits (LCL and UCL) the lower and upper randomization probabilities of a unbiased np-chart
#' @export
#'
#' @author José Guedes, Manuel C. Morais
#'
#' @references Morais, M.C. (2016). An ARL-unbiased np-chart. Economic Quality Control 31, 11-21. http://www.degruyter.com/view/j/eqc.2016.31.issue-1/eqc-2015-0013/eqc-2015-0013.xml
#' @references Morais, M.C., Wittenberg, P., and Cruz, C.J. (2022). The np-chart with 3-sigma limits and the ARL-unbiased np-chart revisited. Stochastics and Quality Control 37, 107-116. https://doi.org/10.1515/eqc-2022-0032
#'
#'
#' @examples
#' n <- 100
#' p0 <- 0.120
#' arl0 <- 370.4
#' binom.unb.crit(n,p0,arl0)
binom.unb.crit <- function(n, p0, arl0) {
  alpha <- 1/arl0
  result_vector <- search.GandC(n,p0,alpha)
  cL <- result_vector[3]
  cU <- result_vector[4]
  gL <- result_vector[5]
  gU <- result_vector[6]
  result <- paste("LCL:", cL, " ", "UCL:", cU," ", "gL:", round(gL,4), " ", "gU:", round(gU,4))
  return(result)
}

n <- 400
p0 <- 0.03
arl0 <- 370.4
p <- 1.1*p0
binom.unb.crit(n,p0, arl0)

####### PARA TABELA TESE #######

# Define the binom.unb.crit function
binom.unb.crit <- function(n, p0, arl0) {
  alpha <- 1 / arl0
  result_vector <- search.GandC(n, p0, alpha)
  cL <- result_vector[3]
  cU <- result_vector[4]
  gL <- result_vector[5]
  gU <- result_vector[6]
  result <- paste("n:", n, "p0:", p0, "(LCL, UCL):", cL, cU, "(gL,gU):", round(gL, 6), round(gU, 6))
  return(result)
}

# Define the values for n and p0
n_values <- c(25, 50, 100, 200, 500, 1000, 1500)
p0_values <- c(0.01, 0.05, 0.1, 0.25, 0.5)

# Set the arl0 value
arl0 <- 370.4

# Initialize an empty list to store results
results_list <- list()

# Nested loop to iterate over each combination of n and p0
for (n in n_values) {
  for (p0 in p0_values) {
    result <- binom.unb.crit(n, p0, arl0)
    # Store the result in the list
    results_list[[paste("n", n, "p0", p0, sep = "_")]] <- result
    # Print the result
    print(result)
  }
}
