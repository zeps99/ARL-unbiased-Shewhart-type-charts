library(qcc)

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


#-----------------------------------------------------------------------------#
#                                                                             #
#                   ADAPTATION OF QUALITY CONTROL CHARTS IN R                             #
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
  
  ## ADDED
  ## Obtaining the randomization probabilities with alpha=0.0027
  object$gaminhas <- gammas(lcl, ucl,sizes,center/sizes,0.0027) ############################
  
  
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
    
    
    ## ADDED    
    ##### Randomization probabilities
    gaminhas <- gammas(cL=lcl,cU=ucl,n=sizes,p0=center/sizes,0.0027) 
    
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
  center <- 3
  pbar <- 0.03
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
  pbar <- 0.03
  limites <- search.GandC(sizes,pbar,0.0027)
  lcl<-limites[,3]
  ucl<-limites[,4]
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
  gaminhas <- gammas(cL=lcl,cU=ucl,n=sizes,p0=pbar,0.0027)
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


##### PLOT  NP-CHART #################
dados  <- c(rbinom(n=10, size=100, .03), rbinom(n=30, size=100, .05))
dados  <- cbind(1:40, dados)
dados  <- data.frame(dados)
grupos <- qcc.groups(dados[, 2], dados[, 1]) 

qcc(dados[1:10,2], type="np", sizes=100, newdata=dados[11:40,2],newsizes=100, ylab="Nonconforming items", xlab="Group", title="ARL-unbiased np-chart")

