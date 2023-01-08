################################################################################
logit <- function(p) log(p/(1-p))
expit <- function(x) 1/(1+exp(-x))
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
rbeta.t <- function(n, shape1, shape2, ncp=0,
                    eps=.Machine$double.neg.eps) {
  # Prevent beta random variates from returning 0 or 1 values
  x <- rbeta(n, shape1, shape2, ncp)
  x <- pmin(1-eps,pmax(eps,x))
  return(x)
}
#############################################################################
rcat <- function(n,prob) {
  # draw n times from a categorical variable with probability prob
  apply(rmultinom(n,1,prob),2,which.max)
}
#############################################################################
