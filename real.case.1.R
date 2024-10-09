##########################################################################
### Zeger1998 - monthly number of cases of poliomyelitis (1970 - 1983) ### 
##########################################################################

library(qcc)

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

gammas <- function(arl0, lambda0, beta0, target="lambda") {
  ARL0 <- arl0
  lambda <- lambda0
  beta <- beta0
  LUg <- inar1.get.UL2(ARL0, lambda, beta, target = "lambda")
  #L  <- LUg$L
  gL <- round(LUg$gL, digits=6)
  #U  <- LUg$U
  gU <- round(LUg$gU, digits=6)
  gg <- data.frame(gL, gU)
  return(gg)
}  

limites <- function(arl0, lambda0, beta0, target="lambda") {
  ARL0 <- arl0
  lambda <- lambda0
  beta <- beta0
  LUg <- inar1.get.UL2(ARL0, lambda, beta, target = "lambda")
  L  <- LUg$L
  #gL <- round(LUg$gL, digits=6)
  U  <- LUg$U
  #gU <- round(LUg$gU, digits=6)
  ll <- data.frame(L, U)
  return(ll)
}  


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
  
  ##obtaining the randomization probabilities with arl0 500
  object$gaminhas <- gammas(arl0 = 500, lambda0 = center, beta0 = 0.06416097, target="lambda") ############################									   
  
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
  
  #if(length(center) == 1)
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
  #mtext("CL", side = 4, at = rev(center)[1], 
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
    
    
    gaminhas <- gammas(arl0 = 500, lambda0 = center, beta0 = 0.06416097, target="lambda") 
    
    
    mtext(paste("gL= ", signif(gaminhas[1,1], 6), sep = ""),
          side = 1, line = top.line+2, adj = 0, at = at.col[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("gU= ", signif(gaminhas[1,2], 6), sep = ""),
          
          
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


# c Chart

stats.c <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) != 1)
    stop("all sizes must be be equal for a c chart")
  statistics <- data
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

sd.c <- function(data, sizes, ...)
{
  data <- as.vector(data)
  std.dev <- sqrt(mean(data))
  return(std.dev)
}

limits.c <- function(center, std.dev, sizes, conf)
{
  sizes <- as.vector(sizes)
  if(length(unique(sizes))==1)
    sizes <- sizes[1]
  limites <- limites(arl0 = 500, lambda0 = center, beta0 = 0.06416097, target="lambda")
  lcl <- limites[1,1]
  ucl <- limites[1,2]
  
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
  
  lambda0 <- center
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.ucl <- seq(along = statistics)[statistics == ucl]
  index.lcl <- seq(along = statistics)[statistics == lcl]
  gaminhas <- gammas(arl0 = 500, lambda0, beta0 = 0.06416097, target="lambda")
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

#################################################################
###### ARL-unbiased Poisson INAR(1) profile with REAL DATA ######
#################################################################


library(readxl)

data <- read_excel("data_poisson.xlsx")
data_mat <- as.matrix(data)
head(data_mat)
length(data_mat)

data_vec <- c(data_mat)
length(data_vec)

# 14 years - assume 14 years in-control
mu1 <- mean(data_vec)
print(mu1)
var(data_vec) #### muito diferente de lambda

set.seed(123)
poisson_data <- rpois(168, lambda = mu1)
hist(data_vec, breaks = 10, freq = TRUE, main = "Histogram with Poisson Fit", xlab = "Data", col = "transparent", border = rgb(0, 0, 1, 0.5), xlim = c(0,10))
hist(poisson_data, breaks = 5, freq = TRUE, col = "transparent", border = rgb(1, 0, 0, 0.5), add = TRUE)
legend("topright", legend = c("Data", "Poisson Fit"), fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), cex = 0.8)

### +- follows a Poisson distribution (preciso de checkar)

### calculate beta -> ver weisstestik2011, pagina 4 (Yule-Walker estimate)

diff <- data_vec - mu1
products <- diff[-1] * diff[-length(diff)]
numerator <- sum(products)
denominator <- sum(diff^2)

beta1 <- numerator / denominator
beta1 #0.06416097

lambda1 <- mu1*(1-beta1)
lambda1

inar1.get.UL2(ARL0 = 500, lambda = lambda1, beta = beta1, target="lambda", OUTPUT=FALSE, eps=1e-6, delta=1e-4)

shiftinmu1 <- 4
shiftinlambda1 <- shiftinmu1*(1-beta1)

newdata1 = rpois(n=20, lambda = lambda1 + shiftinlambda1)
newdata1

obj <- qcc(data_vec[158:168], type="c", newdata = newdata1, center=lambda1, ylab="Number of nonconformities", xlab="Sample number", title="ARL-unbiased Poisson INAR(1)-chart")

