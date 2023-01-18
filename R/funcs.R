#############################################################################
#' Make the DPP Stick Breaking Weights
#' 
#' @param vvec Vector of draws from Beta(1,alpha)
#' 
#' @description Make the DPP stick breaking weights from a set of 
#' draws from Beta(1,alpha)
#' @export
makew <- function(vvec) {
  # make the weight vector
  kmax <- length(vvec)
  wvec <- c(vvec[-kmax],1)*c(1,cumprod(1-vvec[-kmax]))
  return(wvec)
}

#' Make the DPP Beta draws
#' 
#' @param wvec Vector of DPP weights
#' 
#' @description Make the DPP Beta draws from a set of unscaled stick breaking weights
#' 
#' @export
makev <- function(uvec, use.Cpp=TRUE) {
  # make the vvec vector from the unscaled weights uvec
  # note: sum(uvec)==1
  if(use.Cpp) {
    return(makev_c(uvec))
  }
  kmax <- length(uvec)
  # direct version is numerically unstable
  ##vvec <- wvec/c(1, 1-cumsum(wvec[-kmax]))
  vvec[1] <- uvec[1]
  if(kmax>1) {
    cp <- 1
    for(k in 2:(kmax-1)) {
      vvec[k] <- uvec[k]/(cp-uvec[k-1])
      cp <- cp*(1-vvec[k-1])
    }
    vvec[kmax] <- 0.5
  }
  return(vvec)
}

#' Make the DPP Beta draws
#' 
#' @param wvec Vector of DPP weights
#' 
#' @description Make the unscaled weights from DPP Beta draws
#'  
#' @export
makeu <- function(vvec, use.Cpp=TRUE) {
  # make the uvec vector from the DPP beta draws vvec
  # note: sum(uvec)==1
  kmax <- length(uvec)
  uvec <- vvec*c(1,cumprod(1-vvec[1-kmax]))
  uvec[kmax] <- prod(1-vvec[-kmax])
  return(uvec)
}
################################################################################
# Gaussian process

gaussp_varfunc <- function(t, t0=0, sigma=1) {
  # integrated variance function
  (sigma^2)*(t-t0)
}

#' Generate from a Gaussian process
#' 
#' @param tvec Vector of locations
#' @param sigma Standard deviation
#' 
#' @export
rgaussp <- function(tvec, sigma=1) {
  # sample from a Gaussian process at locations tvec
  n <- length(tvec)
  if(any(diff(tvec)<0)) stop("tvec must be non-decreasing")
  # increments
  zvec <- c(0,rnorm(n-1)*sqrt(gaussp_varfunc(tvec[-1],tvec[-n],sigma)))
  # accumulate
  yvec <- cumsum(zvec)
  return(yvec)
}

#' Generate from a Gaussian Bridge
#' 
#' @param tvec Vector of locations at which to simulate
#' @param y0 Value of the process at first tvec
#' @param y1 Value of the process at last tvec (unconstrained if NULL)
#' @param sigma Standard deviation
#' 
#' @export
rgaussp_bridge <- function(tvec, y0=0, y1=NULL, sigma=1) {
  # sample from a Gaussian process from t0 to t1 but given starting
  # and ending values y0 and y1
  n <- length(tvec)
  if(any(diff(tvec)<0)) stop("tvec must be non-decreasing")
  if(is.null(y1)) {
    # no final value specified - just simulate forwards
    yvec <- y0 + rgaussp(tvec, sigma)
  } else {
    # otherwise interleave
    below <- c(NA,rep(1,n-2),NA)
    above <- c(NA,rep(n,n-2),NA)
    ivec <- 1:n
    isamp <- sample(2:(n-1))
    yvec <- c(y0,rep(NA,n-2),y1)
    for(i in isamp) {
      ib <- below[i]
      ia <- above[i]
      dt <- tvec[ia]-tvec[ib]
      dt0 <- tvec[i]-tvec[ib]
      dt1 <- tvec[ia]-tvec[i]
      wb <- dt1/dt
      wa <- dt0/dt
      yvec[i] <- rnorm(1, wb*yvec[ib]+wa*yvec[ia], sigma*sqrt(wa*wb*dt))

      below[i] <- NA
      above[i] <- NA
      kdx <- (ivec<i) & !is.na(above) & (above>i)
      above[kdx] <- pmax(ivec,i)[kdx]
      kdx <- (ivec>i) & !is.na(below) & (below<i)
      below[kdx] <- pmin(ivec,i)[kdx]
    }
  }
  return(yvec)
}

################################################################################
# Sample from a gamma process (alpha(),beta) - constant beta

gammp.alphafunc <- function(t1,t2,alpha=1) {
  # alpha() measure
  alpha*(t2-t1)
}

#' Generate from a Gamma Process
#' 
#' @param tvec Vector of locations at which to simulate
#' @param alpha Scaling of the alpha() measure [constant]
#' @param beta Weighting scale
#' 
#' @export
rgammp <- function(tvec, alpha=1, beta=1) {
  # sample from a Gamma process at locations tvec
  n <- length(tvec)
  if(any(diff(tvec)<0)) stop("tvec must be non-decreasing")
  # increments
  zvec <- c(0,rgamma(n-1,
                     gammp.alphafunc(tvec[-n],tvec[-1],alpha),
                     beta))
  # accumulate
  yvec <- cumsum(zvec)
  return(yvec)
}

#' Generate from a Gamma Bridge
#' 
#' @param tvec Vector of locations at which to simulate
#' @param y0 Value of the process at first tvec
#' @param y1 Value of the process at last tvec (unconstrained if NULL)
#' @param alpha Scaling of the alpha() measure [constant]
#' @param beta Weighting scale
#' 
#' @export
rgammp_bridge <- function(tvec, y0=0, y1=NULL, alpha=1, beta=1) {
  # sample from a Gaussian process from t0 to t1 but given starting
  # and ending values y0 and y1
  n <- length(tvec)
  if(any(diff(tvec)<0)) stop("tvec must be non-decreasing")
  if(is.null(y1)) {
    # no final value specified - just simulate forwards
    yvec <- y0 + rgammp(tvec, alpha, beta)
  } else {
    # otherwise interleave
    below <- c(NA,rep(1,n-2),NA)
    above <- c(NA,rep(n,n-2),NA)
    ivec <- 1:n
    isamp <- sample(2:(n-1))
    yvec <- c(y0,rep(NA,n-2),y1)
    for(i in isamp) {
      ib <- below[i]
      ia <- above[i]
      yvec[i] <- ( yvec[ib]
                   +
                   (yvec[ia]-yvec[ib])*
                     rbeta.t(1,gammp.alphafunc(tvec[ib],tvec[i],alpha),
                               gammp.alphafunc(tvec[i],tvec[ia],alpha)) )
      below[i] <- NA
      above[i] <- NA
      kdx <- (ivec<i) & !is.na(above) & (above>i)
      above[kdx] <- pmax(ivec,i)[kdx]
      kdx <- (ivec>i) & !is.na(below) & (below<i)
      below[kdx] <- pmin(ivec,i)[kdx]
    }
  }
  return(yvec)
}


################################################################################
#' Generate from a Gamma process using the RK Stick breaking construction
#' 
#' @param alpha Scaling of the alpha() measure [constant]
#' @param beta Weighting scale
#' @param g0 Base measure (see Details)
#' @param ... Arguments to g0
#' @param kmax Truncation value
#' 
#' @details \code{g0(n, ...)} is a function that generates \code{n} values from a random
#' distribution with parameters \code{...} .   For example \code{g0=rexp} (the default) generates
#' from an Exponential distribution.
#' 
#' @export
rgammp_sbrk <- function(alpha=1, beta=1, gamma=1, g0=rexp, ..., kmax=100) {
  nc <- 0
  cvec <- numeric()
  while(nc<kmax) {
    cval <- rpois(1,gamma)
    nc <- nc + cval
    cvec <- c(cvec,cval)
  }
  m <- length(cvec)
  cvec[m] <- kmax-sum(cvec[-m])
  thetavec <- g0(kmax, ...)
  evec <- rexp(kmax, beta)
  dvec <- rep(1:m, times=cvec)
  ttvec <- rgamma(kmax,dvec,alpha)
  wvec <- evec*exp(-ttvec)
  retval <- list(kmax=kmax, cvec=cvec, dvec=dvec, evec=evec, ttvec=ttvec,
                 wvec=wvec, thetavec=thetavec, odx=order(thetavec))
  return(retval)
}
################################################################################
#' Generate from a Gamma process using the Dirichlet Process version of
#' the Stick breaking construction
#' 
#' @param alpha Scaling of the alpha() measure [constant]
#' @param beta Weighting scale
#' @param g0 Base measure (see Details)
#' @param ... Arguments to g0
#' @param kmax Truncation value
#' 
#' @details \code{g0(n, ...)} is a function that generates \code{n} values from a random
#' distribution with parameters \code{...} .   For example \code{g0=rexp} (the default) generates
#' from an Exponential distribution.
#' 
#' @export
rgammp_sbdp <- function(alpha=1, beta=1, g0=rexp, ..., kmax=100) {
  gamma <- rgamma(1,alpha,beta)
  vvec <- c(rbeta.t(kmax-1,1,alpha),1)
  thetavec <- g0(kmax, ...)
  cvec <- cumprod(1-vvec)
  uvec <- vvec*c(1,cvec[-kmax])
  uvec[kmax] <- 1-sum(uvec[-kmax])
  wvec <- gamma*uvec
  retval <- list(kmax=kmax, gamma=gamma, vvec=vvec, uvec=uvec,
                 wvec=wvec, thetavec=thetavec)
  return(retval)
}


#' Plot a draw from the gamma process generated by Stick Breaking (RK)
#' 
#' @export
drawrgrk <- function(alpha=1, beta=1, gamma=1, g0=rexp, ..., kmax=100) {
  rg <- rgammp_sbrk(alpha=alpha, beta=beta, gamma=gamma, g0=g0, ..., kmax=kmax)
  nw <- sum(cumsum(rg$wvec)/sum(rg$wvec)<0.99)
  plot(rg$thetavec, rg$wvec,
       main=bquote("(SB)" ~ alpha==.(alpha) ~","~ beta==.(beta) ~","~ gamma==.(gamma) ~","~ n[w]==.(nw) ~","~ S[w]==.(sum(rg$wvec))),
       xlab=bquote(theta), ylab=bquote(w))
  invisible()
}

#' Plot a draw from the gamma process generated by Stick Breaking (Dirichlet Process)
#' 
#' @export
drawrgdp <- function(alpha=1, beta=1, g0=rexp, ..., kmax=100) {
  rg <- rgammp_sbdp(alpha=alpha, beta=beta, g0=g0, ..., kmax=kmax)
  nw <- sum(cumsum(rg$wvec)/sum(rg$wvec)<0.99)
  plot(rg$thetavec, rg$wvec,
       main=bquote("(DP)" ~ alpha==.(alpha) ~","~ beta==.(beta) ~","~ n[w]==.(nw) ~","~ S[w]==.(sum(rg$wvec))),
       xlab=bquote(theta), ylab=bquote(w))
  invisible()
}
########################################################################
