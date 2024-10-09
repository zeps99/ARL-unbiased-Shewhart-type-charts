#######################################################################
### monthly counts storms Lakewood, Colorado (USA) (2014-2018)\2016 ### 
#######################################################################

library(qcc)

####
# ARL of a modified g-chart without randomization for GINAR(1) counts
# lambda and beta below are our p and rho, respectively
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

gammas <- function(arl0, p0, rho0) {
  ARL0 <- arl0
  lambda <- p0
  beta <- rho0
  LUg <- ginar1.get.UL2(ARL0, lambda, beta, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4)
  #L  <- LUg$L
  gL <- round(LUg$gL, digits=6)
  #U  <- LUg$U
  gU <- round(LUg$gU, digits=6)
  gg <- data.frame(gL, gU)
  return(gg)
}  


limites <- function(arl0, p0, rho0) {
  ARL0 <- arl0
  lambda <- p0
  beta <- rho0
  LUg <- ginar1.get.UL2(ARL0, lambda, beta, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4)
  L  <- LUg$L
  #gL <- round(LUg$gL, digits=6)
  U  <- LUg$U
  #gU <- round(LUg$gU, digits=6)
  ll <- data.frame(L, U)
  return(ll)
}

#arl0 <- 200
#p0 <- 0.63
#rho0 <- 0

#gammas(arl0, p0, rho0)
#limites(arl0, p0, rho0)


#-----------------------------------------------------------------------------#
#                                                                             #
#                     ADPATATION OF                                           #
#                     QUALITY CONTROL CHARTS IN R                             #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Luca Scrucca                                                   #
#              Department of Statistics                                       #
#              University of Perugia, ITALY                                   #
#              luca@stat.unipg.it                                             #
#                                                                             #
#-----------------------------------------------------------------------------#

#
#  Main function to create a 'qcc' object
#

qcc <- function(data, type = c("xbar", "R", "S", "xbar.one", "p", "np", "c", "u", "g"), sizes, center, std.dev, limits, data.name, labels, newdata, newsizes, newdata.name, newlabels, nsigmas = 3, confidence.level, rules = shewhart.rules, plot = TRUE, ...)
{
  call <- match.call()
  
  if (missing(data))
    stop("'data' argument is not specified")
  
  if(identical(type, eval(formals(qcc)$type)))
  { type <- as.character(type)[1]
  warning("chart 'type' not specified, assuming \"", type, "\"",
          immediate. = TRUE) }
  if(!exists(paste("stats.", type, sep = ""), mode="function") |
     !exists(paste("sd.", type, sep = ""), mode="function") |
     !exists(paste("limits.", type, sep = ""), mode="function"))
    stop(paste("invalid", type, "control chart. See help(qcc) "))
  
  if (missing(data.name)) 
    data.name <- deparse(substitute(data))
  data <- data.matrix(data)
  if (missing(sizes)) 
  { if (any(type==c("p", "np", "u")))
    stop(paste("sample 'sizes' must be given for a", type, "Chart"))
    else
      sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
  { if (length(sizes)==1)
    sizes <- rep(sizes, nrow(data))
  else if (length(sizes) != nrow(data))
    stop("sizes length doesn't match with data") }
  
  if (missing(labels))
  { if (is.null(rownames(data))) labels <- 1:nrow(data)
  else                         labels <- rownames(data) }
  
  stats <- paste("stats.", type, sep = "")
  if (!exists(stats, mode="function"))
    stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, sizes))
  statistics <- stats$statistics
  if (missing(center)) center <- stats$center
  
  sd <- paste("sd.", type, sep = "")
  if (!exists(sd, mode="function"))
    stop(paste("function", sd, "is not defined!"))
  missing.std.dev <- missing(std.dev)
  if (missing.std.dev)
  { std.dev <- NULL
  std.dev <- switch(type, 
                    "xbar" = { if(any(sizes > 25)) "RMSDF"
                      else                "UWAVE-R" },
                    "xbar.one" = "MR",
                    "R" = "UWAVE-R",
                    "S" = "UWAVE-SD",
                    NULL)
  std.dev <- do.call(sd, list(data, sizes, std.dev)) }
  else 
  { if (is.character(std.dev))
  { std.dev <- do.call(sd, list(data, sizes, std.dev)) }
    else
    { if (!is.numeric(std.dev))
      stop("if provided the argument 'std.dev' must be a method available or a numerical value. See help(qcc).")  }
  }
  
  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")
  
  object <- list(call = call, type = type, 
                 data.name = data.name, data = data, 
                 statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev)
  # check for new data provided and update object
  if (!missing(newdata))
  {   if (missing(newdata.name))
  {newdata.name <- deparse(substitute(newdata))}
    newdata <- data.matrix(newdata)
    if (missing(newsizes))
    { if (any(type==c("p", "np", "u")))
      stop(paste("sample sizes must be given for a", type, "Chart"))
      else
        newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) }
    else
    { if (length(newsizes)==1)
      newsizes <- rep(newsizes, nrow(newdata))
    else if (length(newsizes) != nrow(newdata))
      stop("newsizes length doesn't match with newdata") }
    stats <- paste("stats.", type, sep = "")
    if (!exists(stats, mode="function"))
      stop(paste("function", stats, "is not defined"))
    newstats <- do.call(stats, list(newdata, newsizes))$statistics
    if (missing(newlabels))
    { if (is.null(rownames(newdata)))
    { start <- length(statistics)
    newlabels <- seq(start+1, start+length(newstats)) }
      else
      { newlabels <- rownames(newdata) }
    }
    names(newstats) <- newlabels
    object$newstats <- newstats
    object$newdata  <- newdata
    object$newsizes <- newsizes
    object$newdata.name <- newdata.name
    statistics <- c(statistics, newstats)
    sizes <- c(sizes, newsizes)
  }
  
  conf <- nsigmas
  if (!missing(confidence.level))
    conf <- confidence.level
  if (conf >= 1)
  { object$nsigmas <- conf }
  else
    if (conf > 0 & conf < 1)
    { object$confidence.level <- conf } 
  
  # get control limits
  if (missing(limits))
  { limits <- paste("limits.", type, sep = "")
  if (!exists(limits, mode="function"))
    stop(paste("function", limits, "is not defined"))
  limits <- do.call(limits, list(center = center, std.dev = std.dev,
                                 sizes = sizes, conf = conf)) 
  }
  else 
  { if (!missing.std.dev)
    warning("'std.dev' is not used when limits is given")
    if (!is.numeric(limits))
      stop("'limits' must be a vector of length 2 or a 2-columns matrix")
    limits <- matrix(limits, ncol = 2)
    dimnames(limits) <- list(rep("",nrow(limits)), c("LCL ", "UCL"))
  }
  
  lcl <- limits[,1]
  ucl <- limits[,2]
  
  ##obtaining the randomization probabilities with arl0 = 200
  object$gaminhas <- gammas(arl0 = 200, p0 = 1/center, rho0 = 0.3274798) ############################									   
  
  object$limits <- limits
  if (is.function(rules)) violations <- rules(object)
  else                    violations <- NULL
  object$violations <- violations
  
  class(object) <- "qcc"
  if(plot) plot(object, ...) 
  return(object)
}

print.qcc <- function(x, ...) str(x,1)

summary.qcc <- function(object, digits =  getOption("digits"), ...)
{
  #object <- x   # Argh.  Really want to use 'object' anyway
  cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  data.name <- object$data.name
  type <- object$type
  cat(paste(type, "chart for", data.name, "\n"))
  statistics <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(statistics), digits = digits, ...)
  sizes <- object$sizes
  if(length(unique(sizes))==1)
    sizes <- sizes[1]
  if(length(sizes) == 1)
    cat("\nGroup sample size: ", format(sizes))
  else {
    cat("\nSummary of group sample sizes: ")
    tab <- table(sizes)
    print(matrix(c(as.numeric(names(tab)), tab), 
                 ncol = length(tab), byrow = TRUE, 
                 dimnames = list(c("  sizes", "  counts"),
                                 character(length(tab)))), 
          digits = digits, ...)
  }
  cat("\nNumber of groups: ", length(statistics))
  
  center <- object$center
  if(length(center) == 1)
  { cat("\nCenter of group statistics: ", format(center, digits = digits)) }
  else
  { out <- paste(format(center, digits = digits))
  out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)]      
  out <- paste0(paste(out, collapse = " "), " ...")
  cat("\nCenter of group statistics: ", out, sep = "")
  }
  
  sd <- object$std.dev
  if(length(sd) == 1)
  { cat("\nStandard deviation: ", format(sd, digits = digits), "\n") }
  else
  { out <- paste(format(sd, digits = digits))
  out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)]
  out <- paste0(paste(out, collapse = " "), " ...")
  cat("\nStandard deviation: ", out, "\n", sep = "")
  }
  
  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if (!is.null(newstats)) 
  { cat(paste("\nSummary of group statistics in ", 
              newdata.name, ":\n", sep = ""))
    print(summary(newstats), digits = digits, ...)
    newsizes <- object$newsizes
    if (length(unique(newsizes)) == 1)
      newsizes <- newsizes[1]
    if (length(newsizes) == 1)
      cat("\nGroup sample size: ", format(newsizes))
    else 
    { cat("\nSummary of group sample sizes:\n")
      new.tab <- table(newsizes)
      print(matrix(c(as.numeric(names(new.tab)), new.tab),
                   ncol = length(new.tab), byrow = TRUE, 
                   dimnames = list(c("  sizes", "  counts"),
                                   character(length(new.tab)))), 
            digits = digits, ...)
    }
    cat("\nNumber of groups: ", length(newstats), "\n")
  }
  
  limits <- object$limits
  if (!is.null(limits)) 
  { cat("\nControl limits:\n")
    .printShortMatrix(limits, digits = digits, ...) }
  
  invisible()
}


plot.qcc <- function(x, add.stats = TRUE, chart.all = TRUE, 
                     label.limits = c("LCL ", "UCL"),
                     title, xlab, ylab, ylim, axes.las = 0,
                     digits =  getOption("digits"),
                     restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "qcc")))
    stop("an object of class `qcc' is required")
  
  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- object$center
  stats <- object$statistics
  limits <- object$limits 
  lcl <- limits[,1]
  ucl <- limits[,2]
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  if(chart.all) 
  { statistics <- c(stats, newstats)
  indices <- 1:length(statistics) }
  else
  { if(is.null(newstats))
  { statistics <- stats
  indices <- 1:length(statistics) }
    else
    { statistics <- newstats 
    indices <- seq(length(stats)+1, length(stats)+length(newstats)) }
  }
  
  if (missing(title))
  { if (is.null(newstats))
    main.title <- paste(type, "Chart\nfor", data.name)
  else if (chart.all)
    main.title <- paste(type, "Chart\nfor", data.name, 
                        "and", newdata.name)
  else main.title <- paste(type, "Chart\nfor", newdata.name) 
  }
  else main.title <- paste(title)
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  mar <- pmax(oldpar$mar, c(4.1,4.1,3.1,2.1))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = if(add.stats) pmax(mar, c(7.6,0,0,0)) else mar)
  
  # plot Shewhart chart
  plot(indices, statistics, type="n",
       ylim = if(!missing(ylim)) ylim 
       else range(statistics, limits, center),
       ylab = if(missing(ylab)) "Group summary statistics" else ylab,
       xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  axis(1, at = indices, las = axes.las,
       labels = if(is.null(names(statistics))) 
         as.character(indices) else names(statistics))
  axis(2, las = axes.las)
  box()
  top.line <- par("mar")[3]-length(capture.output(cat(main.title)))
  top.line <- top.line - if(chart.all & (!is.null(newstats))) 0.1 else 0.5
  mtext(main.title, side = 3, line = top.line,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  
  lines(indices, statistics, type = "b", pch=20) 
  
  #if(length(center) == 1)                                    ###### alterado (sem CL no plot)
  #  abline(h = center)
  #else lines(indices, center[indices], type="s")
  
  if(length(lcl) == 1) 
  { abline(h = lcl, lty = 2)
    abline(h = ucl, lty = 2) }
  else 
  { lines(indices, lcl[indices], type="s", lty = 2)
    lines(indices, ucl[indices], type="s", lty = 2) }
  mtext(label.limits, side = 4, at = c(rev(lcl)[1], rev(ucl)[1]), 
        las = 1, line = 0.1, col = gray(0.3), cex = par("cex"))
  #mtext("CL", side = 4, at = rev(center)[1],                     ######## alterado (sem CL no plot)
  #      las = 1, line = 0.1, col = gray(0.3), cex = par("cex"))
  
  if(is.null(qcc.options("violating.runs")))
    stop(".qcc.options$violating.runs undefined. See help(qcc.options).")
  if(length(violations$violating.runs))
  { v <- violations$violating.runs
  if(!chart.all & !is.null(newstats))
  { v <- v - length(stats) 
  v <- v[v>0] }
  points(indices[v], statistics[v], 
         col = qcc.options("violating.runs")$col, 
         pch = qcc.options("violating.runs")$pch) 
  }
  
  if(is.null(qcc.options("beyond.limits")))
    stop(".qcc.options$beyond.limits undefined. See help(qcc.options).")
  if(length(violations$beyond.limits))
  { v <- violations$beyond.limits
  if(!chart.all & !is.null(newstats))
  { v <- v - length(stats) 
  v <- v[v>0] }
  points(indices[v], statistics[v], 
         col = qcc.options("beyond.limits")$col, 
         pch = qcc.options("beyond.limits")$pch) 
  }
  
  if(chart.all & (!is.null(newstats)))
  { len.obj.stats <- length(object$statistics)
  len.new.stats <- length(statistics) - len.obj.stats
  abline(v = len.obj.stats + 0.5, lty = 3)
  mtext(# paste("Calibration data in", data.name),
    "In-Control Data", cex = par("cex")*0.8,
    at = len.obj.stats/2, line = 0, adj = 0.5)
  mtext(# paste("New data in", object$newdata.name),  
    "Out-of-Control Data", cex = par("cex")*0.8, 
    at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
  }
  
  if(add.stats) 
  { 
    # computes the x margins of the figure region
    plt <- par()$plt; usr <- par()$usr
    px <- diff(usr[1:2])/diff(plt[1:2])
    xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
    at.col <- xfig[1] + diff(xfig[1:2])*c(0.10, 0.40, 0.65)
    top.line <- 4.5
    # write info at bottom
    mtext(paste("Number of groups = ", length(statistics), sep = ""), 
          side = 1, line = top.line, adj = 0, at = at.col[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    center <- object$center
    
    if(length(center) == 1)
    { mtext(paste("Center = ", signif(center[1], digits), sep = ""),
            side = 1, line = top.line+1, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    else 
    { mtext("Center is variable",
            side = 1, line = top.line+1, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    
    
    gaminhas <- gammas(arl0 = 200, p0 = 1/center, rho0 = 0.3274798) 
    
    
    mtext(paste("gL = ", signif(gaminhas[1,1], 6), sep = ""),
          side = 1, line = top.line+2, adj = 0, at = at.col[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("gU = ", signif(gaminhas[1,2], 6), sep = ""),
          side = 1, line = top.line+2, adj = 0, at = at.col[2],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    
    
    if(length(unique(lcl)) == 1)
    { mtext(paste("LCL = ", signif(lcl[1], digits), sep = ""), 
            side = 1, line = top.line, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    else 
    { mtext("LCL is variable", 
            side = 1, line = top.line, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    
    if(length(unique(ucl)) == 1)
    { mtext(paste("UCL = ", signif(ucl[1], digits), sep = ""),
            side = 1, line = top.line+1, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats")) 
    }
    else 
    { mtext("UCL is variable", 
            side = 1, line = top.line+1, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    
    if(!is.null(violations))
    { mtext(paste("Number beyond limits =",
                  length(unique(violations$beyond.limits))), 
            side = 1, line = top.line+1, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Randomization signals =",
                  length(unique(violations$violating.runs))), 
            side = 1, line = top.line+2, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
  }
  
  invisible()
}



#
#  Functions used to compute Shewhart charts statistics
#

.qcc.c4 <- function(n)
{ sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2)) }



# p Chart

#stats.p <- function(data, sizes)
#{
#  data <- as.vector(data)
#  sizes <- as.vector(sizes)
#  pbar <- sum(data)/sum(sizes)
#  list(statistics = data/sizes, center = pbar)
#}

#sd.p <- function(data, sizes, ...)
#{
#  data <- as.vector(data)
#  sizes <- as.vector(sizes)
#  pbar <- sum(data)/sum(sizes)
#  std.dev <- sqrt(pbar * (1 - pbar))
#  return(std.dev)
#}

#limits.p <- function(center, std.dev, sizes, conf)
#{ 
#  limits.np(center * sizes, std.dev, sizes, conf) / sizes
#}

# np Chart

#stats.np <- function(data, sizes)
#{
#  data <- as.vector(data)
#  sizes <- as.vector(sizes)
#  pbar <- sum(data)/sum(sizes)
#  center <- sizes * pbar
#  
#  if (length(unique(center)) == 1)
#    center <- center[1]
#  list(statistics = data, center = center)
#}

#sd.np <- function(data, sizes, ...)
#{
#  data <- as.vector(data)
#  sizes <- as.vector(sizes)
#  pbar <- sum(data)/sum(sizes)
#  std.dev <- sqrt(sizes * pbar * (1 - pbar))
#  if (length(unique(std.dev)) == 1)
#    std.dev <- std.dev[1]
#  return(std.dev)
#}

#limits.np <- function(center, std.dev, sizes, conf)
#{ 
#  sizes <- as.vector(sizes)
#  if (length(unique(sizes)) == 1)
#    sizes <- sizes[1]
#  pbar <- mean(center / sizes)
#  limites <- search.GandC(sizes,pbar,0.0027)
#  lcl<-limites[,3]
#  ucl<-limites[,4]
#  limits <- matrix(c(lcl, ucl), ncol = 2)
#  rownames(limits) <- rep("", length = nrow(limits))
#  colnames(limits) <- c("LCL", "UCL")
#  return(limits)
#}

# c Chart

#stats.c <- function(data, sizes)
#{
#  data <- as.vector(data)
#  sizes <- as.vector(sizes)
#  if (length(unique(sizes)) != 1)
#    stop("all sizes must be be equal for a c chart")
#  statistics <- data
#  center <- mean(statistics)
#  list(statistics = statistics, center = center)
#}

#sd.c <- function(data, sizes, ...)
#{
#  data <- as.vector(data)
#  std.dev <- sqrt(mean(data))
#  return(std.dev)
#}

#limits.c <- function(center, std.dev, sizes, conf)
#{
#  sizes <- as.vector(sizes)
#  if(length(unique(sizes))==1)
#    sizes <- sizes[1]
#  limites <- search.GandC(center, 0.0027)
#  lcl <- limites[,4]
#  ucl <- limites[,7]
#  
#  limits <- matrix(c(lcl, ucl), ncol = 2)
#  rownames(limits) <- rep("", length = nrow(limits))
#  colnames(limits) <- c("LCL", "UCL")
#  return(limits)
#}

# u Chart

#stats.u <- function(data, sizes)
#{
#  data <- as.vector(data)
#  sizes <- as.vector(sizes)
#  statistics <- data/sizes
#  center <- sum(sizes * statistics)/sum(sizes)
#  list(statistics = statistics, center = center)
#}

#sd.u <- function(data, sizes, ...)
#{
#  data <- as.vector(data)
#  sizes <- as.vector(sizes)
#  std.dev <- sqrt(sum(data)/sum(sizes))
#  return(std.dev)
#}

#limits.u <- function(center, std.dev, sizes, conf)
#{
#  sizes <- as.vector(sizes)
#  if (length(unique(sizes))==1)
#    sizes <- sizes[1]
#  limits.c(center * sizes, std.dev, sizes, conf) / sizes
#}

# g Chart   (BY ME)

stats.g <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) != 1)
    stop("all sizes must be be equal for a g chart")
  pbar <- 1/(mean(data)+1)
  statistics <- data
  center <- 1/pbar
  list(statistics = statistics, center = center)
}

sd.g <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- 1/(mean(data)+1)
  std.dev <- sqrt(1-pbar)/pbar
  if (length(unique(std.dev)) == 1)
    std.dev <- std.dev[1]
  return(std.dev)
}

limits.g <- function(center, std.dev, sizes, conf)
{
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) == 1)
    sizes <- sizes[1]
  pbar <- 1/center
  limites <- limites(arl0 = 200, p0 = pbar, rho0 = 0.3274798)
  lcl<-limites[,1]
  ucl<-limites[,2]
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}



#
# Functions used to signal points out of control 
#

shewhart.rules <- function(object, limits = object$limits, run.length = qcc.options("run.length"))
{
  # Return a list of cases beyond limits and violating runs
  bl <- beyond.limits(object, limits = limits)
  vr <- violating.runs(object, run.length = run.length)
  list(beyond.limits = bl, violating.runs = vr)
}


beyond.limits <- function(object, limits = object$limits)
{
  # Return cases beyond limits
  statistics <- c(object$statistics, object$newstats) 
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl, index.below.lcl))
}

violating.runs <- function(object, run.length = qcc.options("run.length"), limits=object$limits)
{
  # Return indices of points violating runs
  if(run.length == 0)
    return(numeric())
  center <- object$center
  
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.ucl <- seq(along = statistics)[statistics == ucl]
  index.lcl <- seq(along = statistics)[statistics == lcl]
  gaminhas <- gammas(arl0 = 200, p0 = 1/center, rho0 = 0.3274798) #alterado
  gL<-gaminhas[,1]
  gU <- gaminhas[,2]
  for (i in 1:length(index.ucl)){
    aleatorizacao<-rbinom(1,1,gU)
    if (aleatorizacao[1]==0){index.ucl<-index.ucl[-i]}
  }
  for (i in 1:length(index.lcl)){
    aleatorizacao<-rbinom(1,1,gL)
    if (aleatorizacao[1]==0){index.lcl<-index.lcl[-i]}
  }
  
  return(c(index.ucl, index.lcl))
}



###################################################################
###### ARL-unbiased geometric INAR(1) profile with REAL DATA ######
###################################################################

data1 <- c(0, 2, 0, 12, 9, 2, 6, 5, 8, 0, 0, 3, 0, 0, 0, 0, 4, 1, 0, 4, 0,
        15, 0, 0, 1, 0, 1, 3, 4, 0, 2, 0, 1, 0, 0, 0, 2, 0, 2, 2, 1, 4, 
        6, 2, 5, 1, 0, 0)

p0_yw <- 1/(1+mean(data1))
p0_yw #0.3076923

set.seed(123)
data2 <- rgeom(48, p0_yw)

# Create frequency tables for both datasets
df1 <- as.data.frame(table(data1))
df2 <- as.data.frame(table(data2))

# Rename the columns for merging
colnames(df1) <- c("Value", "Data Counts")
colnames(df2) <- c("Value", "geometric Fit")

# Merge the two datasets by 'Value'
df_combined <- merge(df1, df2, by = "Value", all = TRUE)

# Convert to long format for ggplot
df_long <- reshape2::melt(df_combined, id.vars = "Value", 
                          variable.name = "Dataset", value.name = "Count")

# Create the bar plot with overlayed bars
ggplot(df_long, aes(x = Value, y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_identity(), 
           color = "black", alpha = 0.7, width = 0.7) +  # Overlay bars with transparency
  theme_minimal() +
  scale_fill_manual(values = c("Data Counts" = "#2D2D2D", "geometric Fit" = "#BFBFBF")) +  # Darker gray and lighter gray
  theme(
    axis.title.x = element_blank(),   # Remove x-axis title
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_text(size = 12),    # Keep x-axis text
    axis.text.y = element_text(size = 12),    # Keep y-axis text
    plot.title = element_blank(),     # Remove plot title
    legend.position = c(0.87, 0.87),    # Position legend inside the plot
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_blank(),   # Remove legend title
    legend.text = element_text(size = 8),  # Smaller legend text
    panel.grid.major = element_blank(),    # Remove major grid lines
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.line.y = element_line(color = "black", size = 0.8), # Add the y-axis line
    axis.line.x = element_line(color = "black", size = 0.8)  # Add the x-axis line
  )

##################################
##### YULE-WALKER ESTIMATORS #####
##################################

p0_yw <- 1/(1+mean(data1))
p0_yw #0.3076923

### calculate rho_yw
diff <- data - p0_yw
products <- diff[-1] * diff[-length(diff)]
numerator <- sum(products)
denominator <- sum(diff^2)
rho0_yw <- numerator / denominator
rho0_yw #0.3274798

arl0 <- 200
#gammas(arl0, p0, rho0)
#limites(arl0, p0, rho0)

ginar1.get.UL2(ARL0 = arl0, lambda = p0_yw, beta = rho0_yw, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4)

p1 = 0.5*p0
dados  <- c(rgeom(n=30, prob = p1))
dados
obj <- qcc(data[39:48],type="g", sizes = 1, center = 1/p0, newdata = dados[1:30], newsizes = 1, 
           ylab="Number of conformities", xlab="Sample number", 
           title="ARL-unbiased GINAR(1)-chart")

##########################
##### CLS ESTIMATORS #####
##########################

calculate_rho_cls <- function(X) {
  n <- length(X)
  
  sum_Xt_Xt1 <- sum(X[2:n] * X[1:(n-1)])
  sum_Xt <- sum(X)
  sum_Xt1 <- sum(X[2:n])
  sum_Xt1_sq <- sum(X[2:n] ^ 2)
  sum_Xt1_sum_sq <- (sum(X[2:n]))^2 / n
  
  numerator <- sum_Xt_Xt1 - (sum_Xt * sum_Xt1 / n)
  denominator <- sum_Xt1_sq - sum_Xt1_sum_sq
  
  rho_hat_CLS <- numerator / denominator
  
  return(rho_hat_CLS)
}

rho0_cls <- calculate_rho_cls(data)
rho0_cls #0.08840864

calculate_p_cls <- function(X) {
  T <- length(X)
  rho_hat_CLS <- calculate_rho_cls(X)
  
  sum_Xt <- sum(X)
  sum_Xt1 <- sum(X[2:T])
  
  p_hat_CLS <- (sum_Xt - rho0_cls * sum_Xt1) / (T * (1 - rho0_cls))
  
  return(p_hat_CLS)
}

p0_cls <- calculate_p_cls(data)
p0_cls # 2.25

### OI???? algo está mal com estes estimadores


