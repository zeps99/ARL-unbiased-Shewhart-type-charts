library(qcc)

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

gammas(2, 8122, 0.001, 0.0027)

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
  
  ##obtaining the randomization probabilities with alpha=0.0027
  object$gaminhas <- gammas(lcl, ucl,1/center,0.0027) ############################									   
  
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
  plot(indices, statistics, type="n", log="y",
       ylim = if(!missing(ylim)) ylim 
       else range(statistics, limits, center),
       ylab = if(missing(ylab)) "Group summary statistics" else ylab,
       xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], exp(exp(par("usr")[4])),       ######ALTERADO (y max deixei esse)
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
  
  #if(length(center) == 1)                                           #ALTERADO
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
  #mtext("CL", side = 4, at = rev(center)[1],                        #ALTERADO
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
    "In-Control Data", cex = par("cex")*0.8,               #####ALTERADO
    at = len.obj.stats/2, line = 0, adj = 0.5)
  mtext(# paste("New data in", object$newdata.name),  
    "Out-of-Control Data", cex = par("cex")*0.8,           #####ALTERADO
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
    
    
    gaminhas <- gammas(cL=lcl,cU=ucl,1/center,0.0027)  
    
    
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
  limites <- search.GandC(pbar,0.0027)
  lcl<-limites[,4]
  ucl<-limites[,7]
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
  gaminhas <- gammas(cL=lcl,cU=ucl,1/center,0.0027) #alterado
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


#### ARL-unbiased g-chart with simulated data ####
#target_p = 0.001
#LCL = 2
#UCL = 8122
#gL = 0.406312
#gU = 0.224264


p0 = 0.001
p = 0.1*p0

dados  <- c(rgeom(n=10, prob = p0), rgeom(n=30, prob = p))
obj <- qcc(dados[1:10],type="g", sizes = 1, center = 1/p0, newdata = dados[11:40], newsizes = 1, 
           ylab="Number of conformities", xlab="Sample number", 
           title="ARL-unbiased g-chart (log scale)")

dados

#when a nonconformity = 0, the plot stops working. log(0) -> -inf
