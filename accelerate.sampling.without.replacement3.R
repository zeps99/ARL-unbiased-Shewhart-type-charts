#  Libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(bench)
library(patchwork)
library(gridExtra)

#########################
###### OLD VERSION ######
#########################

# Function P, F and etc
P1 <- function(y, p0) dgeom(y, p0)

F1 <- function(x, p0) pgeom(x, p0)

iF11 <- function(f, p0) qgeom(f, p0)

iF2_1 <- function(f, p0) {
  if ( f < 0 | f > 1 ) {
    stop("error")
  } else {
    if ( f == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(f, p0)
      if ( F1(x, p0) <= f  ) {
        while ( F1(x, p0) <= f ) x <- x + 1
      } else {
        while ( F1(x, p0) > f & x > 0 ) x <- x - 1
        if ( F1(x, p0) <= f ) x <- x + 1
      }
    }
  }
  x
}
iF21 <- Vectorize("iF2_1", "f")

# Function G and etc
G_1 <- function(x, p0) sum( dgeom(0:x, p0) * (0:x) )/(1/p0 - 1)

G1 <- Vectorize("G_1", "x")

iG1_1 <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( G1(x, p0) < g  ) {
        while ( G1(x, p0) < g ) x <- x + 1
      } else {
        while ( G1(x, p0) >= g & x > 0 ) x <- x - 1
        if ( G1(x, p0) < g ) x <- x + 1
      }
    }
  }
  x
}
iG11 <- Vectorize("iG1_1", "g")

iG2_1 <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( G1(x, p0) <= g  ) {
        while ( G1(x, p0) <= g ) x <- x + 1
      } else {
        while ( G1(x, p0) > g & x > 0 ) x <- x - 1
        if ( G1(x, p0) <= g ) x <- x + 1
      }
    }
  }
  x
}
iG21 <- Vectorize("iG2_1", "g")

# Bounds
c.bounds1 <- function(alpha, p0) {
  Lmax1 <- iF21(alpha, p0)
  Lmax2 <- iG21(alpha, p0)
  Umin1 <- iF11(1-alpha, p0)
  Umin2 <- iG11(1-alpha, p0)
  Lmax <- min(Lmax1, Lmax2)
  Umin <- max(Umin1, Umin2)
  
  Lmin1 <- iF11(max(F1(Umin-1,p0)-1+alpha,0), p0)
  Lmin2 <- iG11(max(G1(Umin-1,p0)-1+alpha,0), p0)
  Umax1 <- iF21(min(F1(Lmax,p0)+1-alpha,1), p0)
  Umax2 <- iG21(min(G1(Lmax,p0)+1-alpha,1), p0)
  Lmin <- max(Lmin1, Lmin2)
  Umax <- min(Umax1, Umax2)
  
  data.frame(p0, Lmin, Lmax, Umin, Umax)
}

# Randomization Constants

gammas1 <- function(cL, cU, l0, alpha) {
  gU <- ( (1/l0 - 1)*( alpha - G1(cL-1,l0) - 1+G1(cU,l0) ) - cL*( alpha - F1(cL-1,l0) - 1+F1(cU,l0) ) ) / (dgeom(cU,l0) * (cU-cL))
  gL <- ( cU*( alpha - F1(cL-1,l0) - 1+F1(cU,l0) ) - (1/l0 - 1)*( alpha - G1(cL-1,l0) - 1+G1(cU,l0) ) ) / (dgeom(cL,l0) * (cU-cL))
  gg <- data.frame(gL, gU)
  gg
}

# Identify Admissible Constants
search.GandC1 <- function(l0, alpha, OUTPUT=FALSE) {
  CB <- c.bounds1(alpha, l0)
  results <- NULL
  for ( cL in CB$Lmin:CB$Lmax ) {
    for ( cU in CB$Umin:CB$Umax ) {
      GG <- gammas1(cL, cU, l0, alpha)
      good <- ( 0 <= GG[1] & GG[1] <= 1 & 0 <= GG[2] & GG[2] <= 1 )
      if ( OUTPUT ) cat(paste("cL =", cL, ", cU =", cU, ", gL =", GG[1], ", gU =", GG[2], "\n"))
      if ( good ) {
        results <- rbind(results, data.frame(l0, CB$Lmin, CB$Lmax, cL, CB$Umin, CB$Umax, cU, gL=GG[1], gU=GG[2]))
        break
      }
    }
  }
  results
}

#########################
###### NEW VERSION ######
#########################

# Function P, F and etc
P2 <- function(x, p0) {(1-p0)^x * p0} 

F2 <- function(x, p0) {1 - (1-p0)^(x+1)}

iF12 <- function(f, p0) {
  if (f == p0) {
    stop("error")
  } else {
    ceiling(log(1-f)/log(1-p0)) -1
    }
} 

iF22 <- function(f, p0) {
  if (f == p0) {
    stop("error")
  } else {
    ceiling(log(1-f)/log(1-p0)) -1
  }
}


# Function G and etc
G_2 <- function(x, p0) {1 - (1-p0)^x * (1+p0*x)}

G2 <- Vectorize("G_2", "x")

iG1_2 <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( G2(x, p0) < g  ) {
        while ( G2(x, p0) < g ) x <- x + 1
      } else {
        while ( G2(x, p0) >= g & x > 0 ) x <- x - 1
        if ( G2(x, p0) < g ) x <- x + 1
      }
    }
  }
  x
}
iG12 <- Vectorize("iG1_2", "g")

iG2_2 <- function(g, p0) {
  if ( g < 0 | g > 1 ) {
    stop("error")
  } else {
    if ( g == 1 ) {
      x <- Inf
    } else {
      x <- qgeom(g, p0)
      if ( G2(x, p0) <= g  ) {
        while ( G2(x, p0) <= g ) x <- x + 1
      } else {
        while ( G2(x, p0) > g & x > 0 ) x <- x - 1
        if ( G2(x, p0) <= g ) x <- x + 1
      }
    }
  }
  x
}
iG22 <- Vectorize("iG2_2", "g")

# Bounds

c.bounds2 <- function(alpha, p0) {
  Lmax1 <- iF22(alpha, p0)
  Lmax2 <- iG22(alpha, p0)
  Umin1 <- iF12(1-alpha, p0)
  Umin2 <- iG12(1-alpha, p0)
  Lmax <- min(Lmax1, Lmax2)
  Umin <- max(Umin1, Umin2)
  
  Lmin1 <- iF12(max(F2(Umin-1,p0)-1+alpha,0), p0)
  Lmin2 <- iG12(max(G2(Umin-1,p0)-1+alpha,0), p0)
  Umax1 <- iF22(min(F2(Lmax,p0)+1-alpha,1), p0)
  Umax2 <- iG22(min(G2(Lmax,p0)+1-alpha,1), p0)
  Lmin <- max(Lmin1, Lmin2)
  Umax <- min(Umax1, Umax2)
  
  data.frame(p0, Lmin, Lmax, Umin, Umax)
}


# Randomization Constants
gammas2 <- function(cL, cU, p0, alpha) {
  a <- (1-p0)^cL * p0
  b <- (1-p0)^cU * p0
  c <- (cL*(1-p0)^cL * p0)
  d <- (cU*(1-p0)^cU * p0)
  e <- (alpha-1+(1-p0)^(cL)-(1-p0)^(1+cU))
  f <- ((alpha-1)*((1/p0)-1) + (1/p0)*((1-p0)^(cL)*(1+(-1+cL)*p0) + (1-p0)^(cU)*(-1+p0)*(1+p0*cU)))
  gL <- ((d*e-b*f)/(a*d-b*c))
  gU <- ((a*f-c*e)/(a*d-b*c))
  gg <- c(gL, gU)
  gg
}

# Identify Admissible Constants
search.GandC_sampling <- function(p0, alpha, OUTPUT=FALSE) {
  CB <- c.bounds2(alpha, p0)
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
    
    GG <- gammas2(cL, cU, p0, alpha)
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


p0 <- 0.001
alpha <- 0.0027
search.GandC1(p0,alpha)
search.GandC_sampling(p0,alpha)

#################################################
####### Calculations Times and Jitter Plot ######
#################################################

# Setup
min_iter <- 1000
iter <- 1000

# res1 Contains The Times For Each Method
res1 <- mark(
  iterations=iter, min_iterations=min_iter,
  "search1"=as.numeric(search.GandC1(p0, alpha)[,c("gL", "gU")]),
  "search2"=as.numeric(search.GandC_sampling(p0, alpha,OUTPUT=FALSE)[,c("gL", "gU")])
)


# Extracting The Times For Each Method
times <- c(res1$time[1], res1$time[2])

# Creating Vector With The Method Names
methods <- c("PMK", "PMK-r")

data <- data.frame(Method = rep(methods, each = length(res1$time[[1]])),
                   Time = unlist(times))

# Adding Numeric Variable For x-Axis Position For Each Method
data$MethodNum <- as.numeric(factor(data$Method))

# Adding Jitter For Better Visualization
data$jitter <- jitter(data$MethodNum)

filtered_data <- data %>% filter(Time <= 6)

# Scatter Plot With Rectangles And Points For Each Method
ggplot(filtered_data, aes(x = Time, y = jitter)) +
  geom_point(position = position_jitter(height = 0.2), size = 1, fill = "lightgrey", shape = 21, color = "grey", stroke = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 6)) +
  scale_y_continuous(breaks = 1:length(unique(filtered_data$MethodNum)), labels = unique(filtered_data$Method)) +
  labs(x = "Time (seconds)",
       y = "Method") +
  theme_minimal()

# Calculate The Mean Time For Each One
#meantime1 <- res1$total_time[1]/res1$n_itr[1]
#meantime2 <- res1$total_time[2]/res1$n_itr[2]

#print(meantime1)
#print(meantime2)

# Median time for each one
mediantime1 <- res1$median[1]
mediantime2 <- res1$median[2]

print(mediantime1)
print(mediantime2)

######################################################
###### Boxplots And Histogram (different plots) ######
######################################################

data1 <- subset(filtered_data, MethodNum == 1)
data2 <- subset(filtered_data, MethodNum == 2)

time1 <- as.numeric(data1$Time)
time2 <- as.numeric(data2$Time)

# Layout to split the screen
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), heights = c(0.5, 2))

# Boxplot 1
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(time2, horizontal=TRUE, ylim=c(min(time2), max(time2)), xaxt="n", col="lightgray", frame=F)
title(main = "PMK-r", adj = 0.5, line = -0.3)

# Boxplot 2
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(time1, horizontal=TRUE, ylim=c(min(time1), max(time1)), xaxt="n", col="darkgray", frame=F)
title(main = "PMK", adj = 0.5, line = -0.3)

# Histogram 1
par(mar=c(4, 3.1, 1.1, 2.1))
hist(time2, breaks=40, col="lightgray", border=T, main="", xlab="Time (s)", xlim=c(min(time2), max(time2)))

# Histogram 2
par(mar=c(4, 3.1, 1.1, 2.1))
hist(time1, breaks=40, col="darkgray", border=T, main="", xlab="Times (s)", xlim=c(min(time1), max(time1)))

################################################
###### Boxplots And Histogram (same plot) ######
################################################

# Transform the data to logarithmic scale
log_time1 <- log10(time1)
log_time2 <- log10(time2)

# Determine the range of the combined data
combined_range <- range(c(log_time1, log_time2))

# Define the number of breaks
num_breaks <- 100
common_breaks <- seq(combined_range[1], combined_range[2], length.out = num_breaks + 1)

# Layout to split the screen
layout(matrix(c(1, 2, 3), nrow = 3, byrow = TRUE), heights = c(0.5, 0.5, 4))

# Boxplot 1
par(mar=c(0, 4.1, 0.7, 2.1))  # Adjusted margins
boxplot(log_time1, horizontal=TRUE, ylim=c(floor(min(log_time2)), max(log_time1)), xaxt="n", col=rgb(0.8, 0.8, 0, 0.5), frame=F)
title(main = "Search.GandC", adj = 0.8, line = -0.2)

# Boxplot 2
par(mar=c(0, 4.1, 0.5, 2.1))  # Adjusted margins
boxplot(log_time2, horizontal=TRUE, ylim=c(floor(min(log_time2)), max(log_time1)), xaxt="n", col=rgb(0.2, 0.8, 0.5, 0.5), frame=F)
title(main = "Search.GandC_random", adj = 0.1, line = -0.2)

# Combined Histogram
par(mar=c(5, 4.1, 0, 0.5))  # Adjusted margins for the histogram
hist(log_time1, breaks=common_breaks, col=rgb(0.8, 0.8, 0, 0.5), border=T, main="", xlab="Time (Log scale)", xlim= c(floor(min(log_time2)), max(log_time1)), freq=FALSE)
hist(log_time2, breaks=common_breaks, col=rgb(0.2, 0.8, 0.5, 0.5), border=T, add=TRUE, freq=FALSE)

######################
###### Boxplots ######
######################

par(mar = c(4.5, 1, 3.5, 1))
plot(0, type = "n", xlim = c(floor(min(log_time2)), max(log_time1)), ylim = c(0.5, 2.5), ylab = "", xlab = "Time (Log scale)", axes = FALSE)

boxplot(log_time1, horizontal = TRUE, at = 1, add = TRUE, col = rgb(0.8, 0.8, 0, 0.5), frame = FALSE)
boxplot(log_time2, horizontal = TRUE, at = 2, add = TRUE, col = rgb(0.2, 0.8, 0.5, 0.5), frame = FALSE)
title(main = "Comparison of Search.GandC procedures", line = 2)
legend("topright", legend = c("Search.GandC", "Search.GandC_random"), fill = c(rgb(0.8, 0.8, 0, 0.5), rgb(0.2, 0.8, 0.5, 0.5)), border = "transparent", bg = "transparent", cex = 0.6)

