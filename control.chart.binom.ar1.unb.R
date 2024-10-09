library(qcc)

# ARL of an np-chart (without randomization) falsely assuming independent
# output when dealing with AR(1) binomial counts 

ar1.arl <- function(L, U, n, p, rho) {
  # Ensure L and U are single values # ADDED
  L <- if (length(L) > 1) L[1] else L
  U <- if (length(U) > 1) U[1] else U
  
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
  # Ensure L and U are single values
  L <- if (length(L) > 1) L[1] else L
  U <- if (length(U) > 1) U[1] else U
  
  if (rho<max(-p/(1-p),-(1-p)/p) || rho>1) return(NA)   ## alterei de break para return(NA)
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
  # Ensure n and p are single value # ADDED
  n <- if (length(n) > 1) n[1] else n
  p <- if (length(p) > 1) p[1] else p
  
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

gammas <- function(arl0, n, p0, rho0) {
  ARL0 <- arl0
  p <- p0
  rho <- rho0
  LUg <- ar1.get.UL2(ARL0, n, p, rho, target="p")
  #L  <- LUg$L
  gL <- round(LUg$gL, digits=6)
  #U  <- LUg$U
  gU <- round(LUg$gU, digits=6)
  gg <- data.frame(gL, gU)
  return(gg)
}  


limites <- function(arl0, n, p0, rho0) {
  ARL0 <- arl0
  p <- p0
  rho <- rho0
  LUg <- ar1.get.UL2(ARL0, n, p, rho, target="p")
  L  <- LUg$L
  #gL <- round(LUg$gL, digits=6)
  U  <- LUg$U
  #gU <- round(LUg$gU, digits=6)
  ll <- data.frame(L, U)
  return(ll)
}

#-----------------------------------------------------------------------------#
#                                                                             #
#                   ADAPTATION OF QUALITY CONTROL CHARTS IN R                 #
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
  
  ##obtaining the randomization probabilities with arl0 = 370.4
  object$gaminhas <- gammas(arl0 = 370.4, n = 15, p0 = 0.38, rho0 = 0.97) ############################									   
  
  
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
  #object <- x
  cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  data.name <- object$data.name
  type <- object$type
  cat(paste("chart for", type, "\n"))
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
  
  if(length(center) == 1)
    abline(h = center, lty="dashed") ## Modified
  else lines(indices, center[indices], type="s")
  
  if(length(lcl) == 1) 
  { abline(h = lcl, lty = "solid")## Modified
    abline(h = ucl, lty = "solid") }## Modified
  else 
  { lines(indices, lcl[indices], type="s", lty = 2)
    lines(indices, ucl[indices], type="s", lty = 2) }
  mtext(label.limits, side = 4, at = c(rev(lcl)[1], rev(ucl)[1]), 
        las = 1, line = 0.1, col = gray(0.3), cex = par("cex"))
  mtext("CL", side = 4, at = rev(center)[1], 
        las = 1, line = 0.1, col = gray(0.3), cex = par("cex"))
  
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
  
  ## ADDED
  mtext(# paste("Calibration data in", data.name),
    "In-control data", cex = par("cex")*0.8,
    at = len.obj.stats/2, line = 0, adj = 0.5)
  mtext(# paste("New data in", object$newdata.name),  
    "Out-of-control data", cex = par("cex")*0.8, 
    at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
  ## End of ADDED
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
    sizes <- object$sizes
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
    
    lcl <- limits[,1]
    ucl <- limits[,2]
  
    ##obtaining the randomization probabilities with arl0 = 370.4
    gaminhas <- gammas(arl0 = 370.4, n = 15, p0 = 0.38, rho0 = 0.97) ############################									   
    
    mtext(paste("gL = ", signif(gaminhas[1,1], digits), sep = ""),
          side = 1, line = top.line+2, adj = 0, at = at.col[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("gU = ", signif(gaminhas[1,2], digits), sep = ""),
          side = 1, line = top.line+2, adj = 0, at = at.col[2],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    ## End of ADDED
    
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
    
    ## ADDED? 
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
  ## End of ADDED?
  
  invisible()
}

#
#  Functions used to compute Shewhart charts statistics
#

.qcc.c4 <- function(n)
{ sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2)) }


# np Chart
stats.np <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  center <- sizes * pbar
  ## ADDED
  center <- 5.7
  pbar <- 0.38
  ## End of ADDED
  if (length(unique(center)) == 1)
    center <- center[1]
  list(statistics = data, center = center)
}

sd.np <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(sizes * pbar * (1 - pbar))
  if (length(unique(std.dev)) == 1)
    std.dev <- std.dev[1]
  return(std.dev)
}

limits.np <- function(center, std.dev, sizes, conf)
{ 
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) == 1)
    sizes <- sizes[1]
  pbar <- center / sizes
  ## ADDED
  pbar <- 0.38
  limites <- limites(arl0 = 370.4, n = 15, p0 = pbar, rho0 = 0.97)
  lcl<-limites[,1]
  ucl<-limites[,2]
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  ## End of ADDED
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
  sizes <- object$sizes
  pbar <- center / sizes 
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.ucl <- seq(along = statistics)[statistics == ucl]
  index.lcl <- seq(along = statistics)[statistics == lcl]
  gaminhas <- gammas(arl0 = 370.4, n = 15, p0 = 0.38, rho0 = 0.97)
  gL <- gaminhas[,1]
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

n = 15
p0 = 0.38 
#p1 = 0.18 
ARL0 = 370.4
rho = 0.97

ar1.get.UL2(ARL0, n, p = p0, rho, target="p", OUTPUT=FALSE, eps=1e-6, delta=1e-4)
  
L =2
U = 10
gL = 0.038700 
gU = 0.084907

ar1.arl(L, U, n, p = p0, rho)

gaminhas <- gammas(arl0 = 370.4, n = 15, p0 = 0.38, rho0 = 0.97)
gL <- gaminhas[,1]
gL
gU <- gaminhas[,2]
gU


##### PLOT CONTROL-CHART #################
dados  <- c(rbinom(n=10, size=15, .38), rbinom(n=30, size=15, .18))
dados
sum(dados)

pbar <- 215 * 0.38
pbar

qcc(dados[1:10], type="np", sizes=15, newdata=dados[11:30],newsizes=15, ylab="Nonconforming items", xlab="Group", title="ARL binomial AR(1)-chart")





# Target values, etc.
n    <- 15
p0   <- .38
p1   <- .18
rho0 <- .97
ARL0 <- 370.4

beta0  <- p0*(1-rho0)
alpha0 <- beta0 + rho0
beta1  <- p1*(1-rho0)
alpha1 <- beta1 + rho0

# Control limits and randomization probabilities
L  <- 2
U  <- 10
gL <- .03870007
gU <- .08490665


##### PLOT/code adapted from PaulinoMoraisKnoth2019 #################

# Stuff
nobs       <- 10+29
veclogins  <- rep(NA, nobs+1)
vecalarms1 <- rep(NA, nobs+1)
vecalarms2 <- rep(NA, nobs+1)
vecalarms3 <- rep(NA, nobs+1)
vecalarms  <- rep(NA, nobs+1)

# Simulating data
set.seed(5961)
Z      <- rbinom(1, n, p0)
alarm1 <- ( Z < L  ) | ( Z > U )
alarm2 <- ( Z==L & rbinom(1, 1, gL) )
alarm3 <- ( Z==U & rbinom(1, 1, gU) )
alarm  <- alarm1 | alarm2 | alarm3
vecalarms1[1] <- alarm1
vecalarms2[1] <- alarm2
vecalarms3[1] <- alarm3
vecalarms[1]  <- alarm

ZZ <- Z
i  <- 1
while ( i <= nobs) {
  i <- i+1
  if( i<=10 ) {
    Z  <- rbinom(1, Z, alpha0) + rbinom(1, n-Z, beta0)
  }
  else {
    Z <- rbinom(1, Z, alpha1) + rbinom(1, n-Z, beta1)
  }
  ZZ <- c(ZZ, Z)
  alarm1 <- ( Z < L  ) | ( Z > U )
  alarm2 <- ( Z==L & rbinom(1, 1, gL) )
  alarm3 <- ( Z==U & rbinom(1, 1, gU) )
  alarm <- alarm1 | alarm2 | alarm3
  veclogins[i]  <- Z
  vecalarms1[i] <- alarm1
  vecalarms2[i] <- alarm2
  vecalarms3[i] <- alarm3
  vecalarms[i]  <- alarm
}


cat("nbeds \n", vecbeds, "\n\n")

cat("alarms due to obs. beyond the control limits \n", vecalarms1, "\n\n")

cat("alarms due to obs. equal to the lower control limit and randomization \n", vecalarms2, "\n\n")

cat("alarms due to obs. equal to the upper control limit and randomization \n", vecalarms3, "\n\n")

cat("alarms \n", vecalarms)


plot(ZZ, type="o", pch=19, cex=.8, xlab="t", ylab=expression(X[t]), ylim=c(0,10))
abline(h=c(L,U), col="red", lty=2)
abline(h=n*p0, lty=4, col="forestgreen")
abline(v=10.5, lty=3, col="blue")
points(ZZ, type="o", pch=19, cex=.8)
if ( any(vecalarms3) ) points(which(vecalarms3), rep(U, sum(vecalarms3)), col="red", pch=19, cex=1.2)
