################################################################################
#' @export
logit <- function(p) log(p/(1-p))
#' @export
expit <- function(x) 1/(1+exp(-x))
################################################################################

################################################################################
#' Histogram from count data
#' 
#' @description Plotting a histogram from a set of counts in bins rather
#' than a vector of observations
#' 
#' @export
hist.from.counts <- function(counts, mids=NULL, breaks=NULL, col="light grey", 
                             xlim=NULL, ...) {
  # Define bin breaks and counts
  #breaks <- c(0, 10, 20, 30, 40, 50) # The boundaries of your bins
  #counts <- c(10, 25, 40, 15, 5) # Counts for each bin
  
  if(!is.null(breaks)) {
    # Calculate bin midpoints for potential labeling or other uses
    bin_mids <- (breaks[-length(breaks)] + breaks[-1]) / 2
    if(is.null(xlim)) xlim <- range(breaks)
  } else if(!is.null(mids)) {
    bin_mids <- mids
    dd <- diff(mids)[1]
    breaks <- c(bin_mids[1]-dd/2,bin_mids+dd/2)
    if(is.null(xlim)) xlim <- range(mids) + dd*c(-1,1)
  } else {
    bin_mids <- 1:length(counts)
    breaks <- c(bin_mids[1]-dd/2,bin_mids+dd/2)
    if(is.null(xlim)) xlim <- c(-0.5,length(counts)+0.5)
  }
  
  # Create an empty plot to set up the axes
  plot(NA, xlim = xlim, ylim = c(0, max(counts) * 1.1),
       axes=FALSE, ...)
  axis(1); axis(2)
  
  # Add rectangles for each bin
  for (i in 1:length(counts)) {
    rect(xleft = breaks[i], ybottom = 0, 
         xright = breaks[i+1], ytop = counts[i], col=col, 
         ...)
  }
}
################################################################################


#############################################################################
count.distinct <- function(x) length(unique(x))
vector.eqmatrix <- function(x) {
  matrix(as.numeric(outer(x,x,"==")),nrow=length(x))
}
rowmax <- function(xmat) apply(xmat,1,max)
colmax <- function(xmat) apply(xmat,2,max)
extend.range <- function(y, p=0) {
  # extend the range of y by proportion p
  r <- range(y)
  d <- p*diff(r)
  r <- r + d*c(-1,1)
  return(r)
}
#############################################################################
#' @export
rbeta.t <- function(n, shape1, shape2, ncp=0,
                    eps=.Machine$double.neg.eps) {
  # Prevent beta random variates from returning 0 or 1 values
  x <- rbeta(n, shape1, shape2, ncp)
  x <- pmin(1-eps,pmax(eps,x))
  return(x)
}
#############################################################################
#' @export
rcat <- function(n,prob) {
  # draw n times from a categorical variable with probability prob
  apply(rmultinom(n,1,prob),2,which.max)
}
#############################################################################
