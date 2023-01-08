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
                     rbeta(1,gammp.alphafunc(tvec[ib],tvec[i],alpha),
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
  vvec <- c(rbeta(kmax-1,1,alpha),1)
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
################################################################################
#' Simulate from a hazard rate function 
#' 
#' @param n Number of failure times to simulate
#' @param model.list Model specification list
#' @param model IFR or DFR (increasing or decreasing failure rates)
#' 
#' @description Simulate from a hazard rate function generated by 
#' a set of discrete locations and weights and integrated
#' 
#' @export
rfail <- function(n, model.list, tau=Inf) {
  # simulate failure times, given a set of weights and locations
  # of the underlying process, and a simulation model
  # Censor at time tau
  if(model.list$model=="IFR") {
    tvec <- rfail.ifr(n, model.list$thetavec, model.list$wvec, model.list$lambda0)
  } else if(model.list$model=="DFR") {
    tvec <- rfail.dfr(n, model.list$thetavec, model.list$wvec, model.list$lambda0)
  } else if(model.list$model=="LWB") {
    tvec <- rfail.gen(n, model.list)
  } else if(model.list$model=="SBT") {
    tvec <- rfail.gen(n, model.list)
  } else if(model.list$model=="MBT") {
    tvec <- rfail.mbt(n, model.list$thetavec1, model.list$wvec1, 
                        model.list$thetavec2, model.list$wvec2, 
                        model.list$pival)
  } else if(model.list$model=="LCV") {
    tvec <- rfail.lcv(n, model.list$thetavec, model.list$wvec, model.list$lambda0, model.list$w0)
  } else {
    stop(paste0("Model ",model.list$model," not implemented"))
  }
  tvec <- pmin(tvec, tau)
  return(tvec)
}

#' @export
rfail.gen <- function(n, model.list) {
   # Generate failure times where the integrated hazard is piecewise linear
   if(model.list$model%in%c("IFR","DFR","LWB","SBT")) {
     # Integrated hazard is piecewise linear
     if(model.list$model %in% c("IFR","DFR")) {
        tknot <- sort(unique(model.list$thetavec))
     } else if(model.list$model=="LWB") {
        tknot <- sort(unique(c(model.list$thetavec, model.list$a)))
     } else if(model.list$model=="SBT") {
        tknot <- sort(unique(c(model.list$thetavec1, model.list$thetavec2)))
     } else {
        stop(paste0("Model ",model.list$model," not implemented"))
     }
     kmax <- length(tknot)
     kmaxp2 <- kmax+2
     tknot <- c(0,tknot,1.1*tknot[kmax])
     svec <- int.lambda.func(tknot, model.list)
     nluvec <- -log(runif(n,0,1))
     k1vec <- apply(outer(nluvec, svec, function(nlu,s) nlu>s),1,sum)
     k1vec <- pmin(k1vec,kmaxp2-1)
     mvec <- diff(tknot)/diff(svec)
     tvec <- tknot[k1vec] + mvec[k1vec]*(nluvec-svec[k1vec])
   } else {
     stop(paste0("Model ",model.list$model," not implemented"))
   }
   return(tvec)
}

#' Simulate from an IFR model
#' 
#' @param n Number of failure times to simulate
#' @param thetavec Set of support locations
#' @param wvec Set of associated weights
#' @param lambda0 Offset to be added to the hazard rate
#' 
#' @export
rfail.ifr <- function(n, thetavec, wvec, lambda0=0) {
  # simulate failure times from an IFR model, 
  # given a set of weights and locations
  odx <- order(thetavec)
  othetavec <- thetavec[odx]
  owvec <- wvec[odx]
  kmax <- length(othetavec)
  othetavec <- c(0,othetavec,1.1*othetavec[kmax])
  owvec <- c(0,owvec,0)
  kmaxp2 <- kmax+2
  nluvec <- -log(runif(n,0,1))
  
  s1vec <- cumsum(owvec)
  s2vec <- cumsum(owvec*othetavec)
  svec <- (lambda0+s1vec[-kmaxp2])*othetavec[-1] - s2vec[-kmaxp2]
  k1vec <- apply(outer(nluvec, svec, function(nlu,s) nlu>s),1,sum)
  k1vec <- pmin(k1vec,kmaxp2-1)
  k2vec <- k1vec+1
  tvec <- (nluvec + s2vec[k2vec])/(lambda0 + s1vec[k2vec])
    
  # check 
  #int.lambda.vec <- int.lambda.func(othetavec, thetavec, wvec, lambda0, "IFR")
  #uvec <- exp(-nluvec)
  #all( tvec>othetavec[kmaxp2] | exp(-int.lambda.vec[k2vec+1])<= uvec & uvec <= exp(-int.lambda.vec[k1vec+1]) )
  #all( tvec>othetavec[kmaxp2] | (othetavec[k2vec] <= tvec & tvec <= othetavec[k2vec+1]) )
  #all( nluvec>svec[kmaxp2-1] | svec[k1vec]<= nluvec & nluvec <= svec[k2vec] )
    
  return(tvec)
}

#' Simulate from an DFR model
#' 
#' @param n Number of failure times to simulate
#' @param thetavec Set of support locations
#' @param wvec Set of associated weights
#' @param lambda0 Offset to be added to the hazard rate
#' 
#' @export
rfail.dfr <- function(n, thetavec, wvec, lambda0=0) {
  # simulate failure times from a DFR model, 
  # given a set of weights and locations
  odx <- order(thetavec)
  othetavec <- thetavec[odx]
  owvec <- wvec[odx]
  kmax <- length(othetavec)
  othetavec <- c(0,othetavec,1.1*othetavec[kmax])
  owvec <- c(0,owvec,0)
  kmaxp2 <- kmax+2
  nluvec <- -log(runif(n,0,1))
  
  s1vec <- cumsum(owvec)
  s2vec <- cumsum(owvec*othetavec)
  s1max <- s1vec[kmaxp2]
  s2max <- s2vec[kmaxp2]
  svec <- (lambda0+s1max-s1vec)*othetavec + s2vec
  k1vec <- apply(outer(nluvec, svec, function(nlu,s) nlu>s),1,sum)
  k1vec <- pmin(k1vec,kmaxp2-1)
  k2vec <- k1vec+1
  tvec <- (nluvec-s2vec[k1vec])/(lambda0+s1max-s1vec[k1vec])
  
  # check 
  #all(k2vec==kmaxp2 | svec[k1vec]<= nluvec & nluvec <= svec[k2vec])
  #int.lambda.vec <- int.lambda.func(othetavec, thetavec, wvec, lambda0, "DFR")
  #range(svec-int.lambda.vec)
  #uvec <- exp(-nluvec)
  #all(k2vec==kmaxp2 | svec[k1vec]<= nluvec & nluvec <= svec[k2vec])
  #all(k2vec==kmaxp2 | othetavec[k1vec] <= tvec & tvec <= othetavec[k2vec])
  #k2vec==kmaxp2 | othetavec[k1vec] <= tvec & tvec <= othetavec[k2vec]
  #all(k2vec==kmaxp2 | (exp(-int.lambda.vec[k2vec])<= uvec & uvec <= exp(-int.lambda.vec[k1vec])))
  return(tvec)
}

#' @export
rfail.mbt <- function(n, thetavec1, wvec1, thetavec2, wvec2, pival=0.5) {
  # simulate failure times from an MBT model 
  
  # Numbers selected from DFR and IFR components respectively
  n1 <- rbinom(1, n, pival)
  n2 <- n-n1
  tvec1 <- rfail.dfr(n1, thetavec1, wvec1, 0)
  tvec2 <- rfail.ifr(n2, thetavec2, wvec2, 0)
  # Combine and shuffle
  tvec <- sample(c(tvec1,tvec2))
  return(tvec)
}  

#' @export
rfail.lcv <- function(n, thetavec, wvec, lambda0=1, w0=-1) {
  # simulate failure times from an LCV model
  odx <- order(thetavec) 
  othetavec <- thetavec[odx]  # theta*_k (1...K)
  owvec <- wvec[odx]          # w*_k     (1...K)
  kmax <- length(othetavec)   # K
  othetavec <- c(0,othetavec,1.1*othetavec[kmax]) # 0, thetavec, 1.1*last
  owvec <- c(0,owvec,0)                           # 0, wvec,     0
  kmaxp2 <- kmax+2            # K+2
  
  s1vec <- cumsum(owvec)
  s01vec <- w0 + s1vec
  s2vec <- cumsum(owvec*othetavec)
  ccvec <- lambda0*exp(-s2vec)/s01vec
  ccvec1 <- ccvec[-kmaxp2]
  s01vec1 <- s01vec[-kmaxp2]
  csvec1 <- ccvec1*(exp(s01vec1*othetavec[-1])-exp(s01vec1*othetavec[-kmaxp2]))
  ssvec <- c(0,cumsum(csvec1))
  
  nluvec <- -log(runif(n,0,1))
  k1vec <- apply(outer(nluvec, ssvec, function(nlu,ilamval) ilamval<=nlu),1,sum)
  k1vec <- pmax(1,k1vec)
  
  tvec <- (1/s01vec[k1vec])*log(exp(s01vec[k1vec]*othetavec[k1vec])
                                + (nluvec-ssvec[k1vec])/ccvec[k1vec])
  return(tvec)
}  

#' Old version of failure simulations
#' 
#' @export
rfail.old <- function(n, thetavec, wvec, lambda0=0, model="IFR") {
  # simulate failure times, given a set of weights and locations
  # of the underlying process, and a simulation model
  odx <- order(thetavec)
  othetavec <- thetavec[odx]
  owvec <- wvec[odx]
  kmax <- length(othetavec)
  othetavec <- c(0,othetavec,1.1*othetavec[kmax])
  owvec <- c(0,owvec,0)
  kmaxp2 <- kmax+2
  nluvec <- -log(runif(n,0,1))
  if(model=="IFR") {
    s1vec <- cumsum(owvec)
    s2vec <- cumsum(owvec*othetavec)
    svec <- (lambda0+s1vec)*othetavec - s2vec
    k1vec <- apply(outer(nluvec, svec, function(nlu,s) nlu>s),1,sum)
    k1vec <- pmin(k1vec,kmaxp2-1)
    k2vec <- k1vec+1
    tvec <- (othetavec[k1vec]
             + (othetavec[k2vec]-othetavec[k1vec])
             /(svec[k2vec]-svec[k1vec]) * (nluvec-svec[k1vec])
    )
  } else if(model=="DFR") {
    s1vec <- cumsum(owvec)
    s2vec <- cumsum(owvec*othetavec)
    s1max <- s1vec[kmaxp2]
    s2max <- s2vec[kmaxp2]
    svec <- (lambda0+s1max-s1vec)*othetavec + s2vec
    k1vec <- apply(outer(nluvec, svec, function(nlu,s) nlu>s),1,sum)
    k1vec <- pmin(k1vec,kmaxp2-1)
    k2vec <- k1vec+1
    tvec <- (othetavec[k1vec]
             + (othetavec[k2vec]-othetavec[k1vec])
             /(svec[k2vec]-svec[k1vec]) * (nluvec-svec[k1vec])
    )
  } else {
    stop(paste0("Model ",model," not implemented"))
  }
  return(tvec)
}

###########################################################################################################
#' Hazard rate function 
#' 
#' @param tvec Locations at which to evaluate the function
#' @param model.list Model specification
#' @param use.Cpp Use cpp to evaluate the function?
#' 
#' @description Hazard rate function generated by 
#' a set of discrete locations and weights and integrated
#' 
#' @export
lambda.func <- function(tvec, model.list, use.Cpp=TRUE) {
  if(model.list$model=="IFR") {
    return(lambda.func.ifr(tvec, model.list$thetavec, model.list$wvec, model.list$lambda0, use.Cpp=use.Cpp))
  } else if(model.list$model=="DFR") {
    return(lambda.func.dfr(tvec, model.list$thetavec, model.list$wvec, model.list$lambda0, use.Cpp=use.Cpp))
  } else if(model.list$model=="LWB") {
    return(lambda.func.lwb(tvec, model.list$thetavec, model.list$wvec, model.list$a, use.Cpp=use.Cpp))
  } else if(model.list$model=="SBT") {
    return(lambda.func.sbt(tvec, 
                           model.list$thetavec1, model.list$wvec1, 
                           model.list$thetavec2, model.list$wvec2, 
                           model.list$lambda0, use.Cpp=use.Cpp))
  } else if(model.list$model=="MBT") {
    return(lambda.func.mbt(tvec, 
                           model.list$thetavec1, model.list$wvec1, 
                           model.list$thetavec2, model.list$wvec2, 
                           model.list$pival, use.Cpp=use.Cpp))
  } else if(model.list$model=="LCV") {
    return(lambda.func.lcv(tvec, 
                           model.list$thetavec, model.list$wvec, 
                           model.list$lambda0, model.list$w0, 
                           use.Cpp=use.Cpp))
  } else {
    stop(paste0("Model ",model.list$model," not implemented"))
  }
}

#' Hazard rate function - IFR case
#' 
#' @param tvec Locations at which to evaluate the function
#' @param thetavec Set of support locations
#' @param wvec Set of associated weights
#' @param lambda0 Offset to be added to the hazard rate
#' @param use.Cpp Use C++ function?
#' 
#' @export
lambda.func.ifr <- function(tvec, thetavec, wvec, lambda0=0, use.Cpp=TRUE) {
  # hazard rate function - IFR case
  if(use.Cpp) {
    return(lambda_func_ifr_c(tvec, thetavec, wvec, lambda0))
  }
  odx <- order(thetavec)
  othetavec <- thetavec[odx]
  owvec <- wvec[odx]
  kmax <- length(othetavec)
  othetavec <- c(0,othetavec,1.1*othetavec[kmax])
  owvec <- c(0,owvec,0)
  kmaxp2 <- kmax+2
  k1vec <- apply(outer(tvec, othetavec, function(t,theta) theta<=t),1,sum)
  k1vec <- pmax(1,pmin(k1vec,kmaxp2-1))
  
  s1vec <- cumsum(owvec)
  lambda.vec <- lambda0 + s1vec[k1vec]
  
  return(lambda.vec)
}

#' @export
lambda.func.dfr <- function(tvec, thetavec, wvec, lambda0=0, use.Cpp=TRUE) {
  # hazard rate function
  if(use.Cpp) {
    return(lambda_func_dfr_c(tvec, thetavec, wvec, lambda0))
  }
  odx <- order(thetavec)
  othetavec <- thetavec[odx]
  owvec <- wvec[odx]
  kmax <- length(othetavec)
  othetavec <- c(0,othetavec,1.1*othetavec[kmax])
  owvec <- c(0,owvec,0)
  kmaxp2 <- kmax+2
  k1vec <- apply(outer(tvec, othetavec, function(t,theta) theta<=t),1,sum)
  k1vec <- pmax(1,pmin(k1vec,kmaxp2-1))
  
  s1vec <- cumsum(owvec)
  lambda.vec <- lambda0 + s1vec[kmaxp2] - s1vec[k1vec]
  
  return(lambda.vec)
}

#' @export
lambda.func.lwb <- function(tvec, thetavec, wvec, a, use.Cpp=TRUE) {
  # hazard rate function - Lo+Weng Bathtub
  if(use.Cpp) {
    return(lambda_func_lwb_c(tvec, thetavec, wvec, a))
  }
  lambda.vec <- as.vector(outer(tvec, thetavec, function(t,theta) theta<abs(t-a))%*%wvec)

  return(lambda.vec)
}

#' @export
lambda.func.sbt <- function(tvec, thetavec1, wvec1, thetavec2, wvec2, lambda0=0, use.Cpp=TRUE) {
  # hazard rate function - Superposition Bathtub
  lambda.vec <- ( lambda0 
                  + lambda.func.dfr(tvec, thetavec1, wvec1, 0, use.Cpp=use.Cpp)
                  + lambda.func.ifr(tvec, thetavec2, wvec2, 0, use.Cpp=use.Cpp) )
  return(lambda.vec)
}

#' @export
lambda.func.mbt <- function(tvec, thetavec1, wvec1, thetavec2, wvec2, pival, use.Cpp=TRUE) {
  # hazard rate function - Mixture Bathtub
  h1vec <- lambda.func.dfr(tvec, thetavec1, wvec1, 0, use.Cpp=use.Cpp)
  h2vec <- lambda.func.ifr(tvec, thetavec2, wvec2, 0, use.Cpp=use.Cpp)
  ih1vec <- int.lambda.func.dfr(tvec, thetavec1, wvec1, 0, use.Cpp=use.Cpp)
  ih2vec <- int.lambda.func.ifr(tvec, thetavec2, wvec2, 0, use.Cpp=use.Cpp)
  fs1vec <- exp(-ih1vec)
  fs2vec <- exp(-ih2vec)
  lambda.vec <- (pival*h1vec*fs1vec + (1-pival)*h2vec*fs2vec)/(pival*fs1vec+(1-pival)*fs2vec)
  return(lambda.vec)
}

#' @export
lambda.func.lcv <- function(tvec, thetavec, wvec, lambda0=1, w0=0, use.Cpp=TRUE) {
  # hazard rate function - Log Convex
  lambda.vec <- lambda0*exp(int.lambda.func.ifr(tvec, thetavec, wvec, w0, use.Cpp=use.Cpp))
  return(lambda.vec)
}

#' Integrated hazard rate function 
#' 
#' @param tvec Locations at which to evaluate the function
#' @param model.list Model specification
#' @param use.Cpp Use cpp to evaluate the function?
#' 
#' @description Hazard rate function generated by 
#' a set of discrete locations and weights and integrated
#' 
#' @export
int.lambda.func <- function(tvec, model.list, use.Cpp=TRUE) {
  if(model.list$model=="IFR") {
    return(int.lambda.func.ifr(tvec, model.list$thetavec, model.list$wvec, model.list$lambda0, use.Cpp=use.Cpp))
  } else if(model.list$model=="DFR") {
    return(int.lambda.func.dfr(tvec, model.list$thetavec, model.list$wvec, model.list$lambda0, use.Cpp=use.Cpp))
  } else if(model.list$model=="LWB") {
    return(int.lambda.func.lwb(tvec, model.list$thetavec, model.list$wvec, model.list$a, use.Cpp=use.Cpp))
  } else if(model.list$model=="SBT") {
    return(int.lambda.func.sbt(tvec, 
                               model.list$thetavec1, model.list$wvec1, 
                               model.list$thetavec2, model.list$wvec2, 
                               model.list$lambda0, use.Cpp=use.Cpp))
  } else if(model.list$model=="MBT") {
    return(int.lambda.func.mbt(tvec, 
                               model.list$thetavec1, model.list$wvec1, 
                               model.list$thetavec2, model.list$wvec2, 
                               model.list$pival, use.Cpp=use.Cpp))
  } else if(model.list$model=="LCV") {
    return(int.lambda.func.lcv(tvec, 
                               model.list$thetavec, model.list$wvec, 
                               model.list$lambda0, model.list$w0, 
                               use.Cpp=use.Cpp))
  } else {
    stop(paste0("Model ",model.list$model," not implemented"))
  }
}

#' Integrated hazard rate function - IFR case
#' 
#' @param tvec Locations at which to evaluate the function
#' @param thetavec Set of support locations
#' @param wvec Set of associated weights
#' @param lambda0 Offset to be added to the hazard rate
#' @param use.Cpp Use cpp to evaluate the function?
#' 
#' @description Integrated hazard rate function generated by 
#' a set of discrete locations and weights and integrated
#' 
#' @export
int.lambda.func.ifr <- function(tvec, thetavec, wvec, lambda0=0, use.Cpp=TRUE) {
  # integrated hazard rate function - IFR case
  if(use.Cpp) {
    return(int_lambda_func_ifr_c(tvec, thetavec, wvec, lambda0))
  }
  odx <- order(thetavec)
  othetavec <- thetavec[odx]
  owvec <- wvec[odx]
  kmax <- length(othetavec)
  othetavec <- c(0,othetavec,1.1*othetavec[kmax])
  owvec <- c(0,owvec,0)
  kmaxp2 <- kmax+2
  k1vec <- apply(outer(tvec, othetavec, function(t,theta) theta<=t),1,sum)
  k1vec <- pmax(1,pmin(k1vec,kmaxp2-1))
  k2vec <- k1vec+1
  
  s1vec <- cumsum(owvec)
  s2vec <- cumsum(owvec*othetavec)
  i1vec <- (lambda0 + s1vec[k1vec])*othetavec[k1vec] - s2vec[k1vec]
  i2vec <- (lambda0 + s1vec[k2vec])*othetavec[k2vec] - s2vec[k2vec]
  int.lambda.vec <- (i1vec
                      +(i2vec-i1vec)/(othetavec[k2vec]-othetavec[k1vec])
                        *(tvec-othetavec[k1vec])
                    )

  return(int.lambda.vec)
}

#' Integrated hazard rate function - DFR case
#' 
#' @param tvec Locations at which to evaluate the function
#' @param thetavec Set of support locations
#' @param wvec Set of associated weights
#' @param lambda0 Offset to be added to the hazard rate
#' @param use.Cpp Use cpp to evaluate the function?
#' 
#' @description Integrated hazard rate function generated by 
#' a set of discrete locations and weights and integrated
#' 
#' @export
int.lambda.func.dfr <- function(tvec, thetavec, wvec, lambda0=0, use.Cpp=TRUE) {
  # integrated hazard rate function - DFR case
  if(use.Cpp) {
    return(int_lambda_func_dfr_c(tvec, thetavec, wvec, lambda0))
  }
  odx <- order(thetavec)
  othetavec <- thetavec[odx]
  owvec <- wvec[odx]
  kmax <- length(othetavec)
  othetavec <- c(0,othetavec,1.1*othetavec[kmax])
  owvec <- c(0,owvec,0)
  kmaxp2 <- kmax+2
  k1vec <- apply(outer(tvec, othetavec, function(t,theta) theta<=t),1,sum)
  k1vec <- pmax(1,pmin(k1vec,kmaxp2-1))
  k2vec <- k1vec+1

  s1vec <- cumsum(owvec)
  s2vec <- cumsum(owvec*othetavec)
  i1vec <- ( (lambda0 + s1vec[kmaxp2] - s1vec[k1vec])*othetavec[k1vec]
               + s2vec[k1vec] )
  i2vec <- ( (lambda0 + s1vec[kmaxp2] - s1vec[k2vec])*othetavec[k2vec]
               + s2vec[k2vec] )
  int.lambda.vec <- (i1vec
                       +(i2vec-i1vec)/(othetavec[k2vec]-othetavec[k1vec])
                       *(tvec-othetavec[k1vec])
                     )

  return(int.lambda.vec)
}

#' @export
int.lambda.func.lwb <- function(tvec, thetavec, wvec, a, use.Cpp=TRUE) {
  # integrated hazard rate function - Lo-Weng Bathtub
  if(use.Cpp) {
    return(int_lambda_func_lwb_c(tvec, thetavec, wvec, a))
  }
  int.lambda.vec <- as.vector(outer(tvec, thetavec, 
                    function(t,theta) {
                      ifelse(t<=a, (theta<a)*pmin(t,a-theta),
                                   pmax(0,a-theta)+pmax(0,t-a-theta))
                    })%*%wvec)
  return(int.lambda.vec)
}

#' @export
int.lambda.func.sbt <- function(tvec, thetavec1, wvec1, thetavec2, wvec2, lambda0=0, use.Cpp=TRUE) {
  # integrated hazard rate function - Superposition Bathtub
  int.lambda.vec <- ( lambda0*tvec 
                      + int.lambda.func.dfr(tvec, thetavec1, wvec1, 0, use.Cpp=use.Cpp)
                      + int.lambda.func.ifr(tvec, thetavec2, wvec2, 0, use.Cpp=use.Cpp) )
  return(int.lambda.vec)
}

#' @export
int.lambda.func.mbt <- function(tvec, thetavec1, wvec1, thetavec2, wvec2, pival, use.Cpp=TRUE) {
  # integrated hazard rate function - Mixture Bathtub
  ih1vec <- int.lambda.func.dfr(tvec, thetavec1, wvec1, 0, use.Cpp=use.Cpp)
  ih2vec <- int.lambda.func.ifr(tvec, thetavec2, wvec1, 0, use.Cpp=use.Cpp)
  fs1vec <- exp(-ih1vec)
  fs2vec <- exp(-ih2vec)
  int.lambda.vec <- -log(pival*fs1vec+(1-pival)*fs2vec)
  return(int.lambda.vec)
}

#' @export
int.lambda.func.lcv <- function(tvec, thetavec, wvec, lambda0=1, w0=0, use.Cpp=TRUE) {
  # integrated hazard rate function - Log Convex

  if(use.Cpp) {
    return(int_lambda_func_lcv_c(tvec, thetavec, wvec, lambda0, w0))
  }
  odx <- order(thetavec) 
  othetavec <- thetavec[odx]  # theta*_k (1...K)
  owvec <- wvec[odx]          # w*_k     (1...K)
  kmax <- length(othetavec)   # K
  othetavec <- c(0,othetavec,1.1*othetavec[kmax]) # 0, thetavec, 1.1*last
  owvec <- c(0,owvec,0)                           # 0, wvec,     0
  kmaxp2 <- kmax+2            # K+2
  k1vec <- apply(outer(tvec, othetavec, function(t,theta) theta<=t),1,sum)
  k1vec <- pmax(1,k1vec)

  s1vec <- cumsum(owvec)
  s01vec <- w0 + s1vec
  s2vec <- cumsum(owvec*othetavec)
  ccvec <- lambda0*exp(-s2vec)/s01vec
  ccvec1 <- ccvec[-kmaxp2]
  s01vec1 <- s01vec[-kmaxp2]
  csvec1 <- ccvec1*(exp(s01vec1*othetavec[-1])-exp(s01vec1*othetavec[-kmaxp2]))
  ssvec <- c(0,cumsum(csvec1))
  
  int.lambda.vec <- ( ssvec[k1vec] 
                      + ccvec[k1vec]*( exp(s01vec[k1vec]*tvec)
                                      -exp(s01vec[k1vec]*othetavec[k1vec])) )

  return(int.lambda.vec)
}



#' Evaluate both the hazard and integrated hazard functions at the same locations
#' 
#' @param tvec Locations at which to evaluate the functions
#' @param model.list Model specification
#' @param use.Cpp Use cpp to evaluate the functions?
#' 
#' @description Hazard rate and Integrated hazard rate function generated by 
#' a set of discrete locations and weights and integrated
#' 
#' @export
both.lambda.funcs <- function(tvec, model.list, use.Cpp=TRUE) {
  lambda.vec <- lambda.func(tvec=tvec, model.list=model.list, use.Cpp=use.Cpp) 
  int.lambda.vec <- int.lambda.func(tvec=tvec, model.list=model.list, use.Cpp=use.Cpp) 
  return(list(lambda.vec=lambda.vec,int.lambda.vec=int.lambda.vec))
}


llike.fail <- function(tvec, model.list, tau=Inf, use.Cpp=TRUE) {
  # log likelihood of observations tvec (censored at tau)
  n <- length(tvec)
  censored <- tvec>=tau
  n0 <- n-sum(censored)
  tvec <- pmin(tvec,tau)

  lambda.vec <- lambda.func(tvec, model.list, use.Cpp=use.Cpp)
  int.lambda.vec <- int.lambda.func(tvec, model.list, use.Cpp=use.Cpp)

  retval <- sum(log(lambda.vec[!censored])) - sum(int.lambda.vec)
  return(retval)
}
################################################################################
#' Update the current state
#' 
#' @export
update.state <- function(state, datlist, fpar, ppar, model) {
  # update the state
  old.state <- state
  
  if(model%in%c("IFR","DFR")) { ##!!== not completed
    # update alpha ##!!==
    if(ppar$update["alpha"]) {
      alpha.old <- state$alpha
      a1star <- fpar$a1+fpar$kmax  ## Check: should be a1+kmax-1 ##!!==
      a2star <- ( fpar$a2-log(state$beta)-log(state$gamma)
                  -log(state$uvec[fpar$kmax]) )
      # log normal random walk proposal
      alpha.new <- exp( rnorm(1,log(alpha.old),ppar$sd.log.alpha))
      log.r <- ( lgamma(alpha.old)-lgamma(alpha.new)
                 + a1star*log(alpha.new/alpha.old)
                 - (alpha.new-alpha.old)*a2star )
      # gamma(a1star,a2star) proposal
      #alpha.new <- rgamma(1,a1star,a2star)
      #log.r <- lgamma(alpha.old)-lgamma(alpha.new)
      if(ppar$verbose) {
        cat(sprintf("alpha: %g->%g: logr=%g\n",
                    alpha.old, alpha.new, log.r))
      }
      if(runif(1)<exp(log.r)) {
        # accepted
        state$alpha <- alpha.new
        state$accepted["alpha"] <- 1
      } else {
        # reject
        state$alpha <- alpha.old
        state$accepted["alpha"] <- 0
      }
    }
    
    # update beta [OK]
    if(ppar$update["beta"]) {
      # Gibbs update
      b1star <- fpar$b1+state$alpha
      b2star <- fpar$b2+state$gamma
      state$beta <- rgamma(1, b1star, b2star)
      state$accepted["beta"] <- 1
    }
    
    # update phi [OK]
    if(ppar$update["phi"]) {
      # Gibbs update
      f1star <- fpar$f1+fpar$kmax
      f2star <- fpar$f2+sum(state$thetavec)
      state$phi <- rgamma(1, f1star, f2star)
      state$accepted["phi"] <- 1
    }
    
    # update lambda0 ##!!==
    if(ppar$update["lambda0"]) {
      state$lambda0 <- state$lambda0
    }
    
    # update gamma
    if(ppar$update["gamma"]) {
      cc <- sum(int.lambda.func(datlist$tvec, 
                                model.list=list(model=model, 
                                                thetavec=state$thetavec, 
                                                wvec=state$uvec, # NB - unscaled weights used here
                                                lambda0=state$lambda0),
                                fpar$use.Cpp))
      g1star <- state$alpha+datlist$n0
      g2star <- state$beta+cc
      gamma.new <- rgamma(1, g1star, g2star)
      if(state$lambda0==0) {
        # Gibbs update [OK]
        state$gamma <- gamma.new
        state$accepted["gamma"] <- 1
      } else {
        # MH update ##!!==
        llike.old <- llikef(state, datlist, fpar, model)
        lprior.old <- lpriorf(state, fpar, model)
        gamma.old <- state$gamma
        state$gamma <- gamma.new
        llike.new <- llikef(state, datlist, fpar, model)
        lprior.new <- lpriorf(state, fpar, model)
        lq.new <- dgamma(gamma.new, g1star, g2star, log=TRUE)
        lq.old <- dgamma(gamma.old, g1star, g2star, log=TRUE)
        log.r <- ( (llike.new+lprior.new+lq.old)
                   -(llike.old+lprior.old+lq.new) )
        #cat(sprintf("Gamma: %f->%f: r=%f\n",gamma.old, gamma.new, exp(log.r)))
        #cat(sprintf("q(.|a,b) = %f %f\n",
        #            state$alpha+datlist$n0, state$beta+cc))
        #cat(sprintf("  (llike,prior,lq): (%.3f %.3f %.3f)-> (%.3f %.3f %.3f)\n",
        #            llike.old, lprior.old, lq.new,
        #            llike.new, lprior.new, lq.old))
        if(runif(1)<exp(log.r)) {
          # accepted
          state$accepted["gamma"] <- 1
        } else {
          # reject
          state$gamma <- gamma.old
          state$accepted["gamma"] <- 0
        }
      }
    }
    
    # update vvec ##!!==
    if(ppar$update["vvec"]) {
      # logistic Normal proposal
      if(ppar$ksweep) {
        # sweep update
        ksamplevec <- 1:fpar$kmax
      } else {
        # random support point update
        ksamplevec <- sample(fpar$kmax, ppar$ksim,
                             prob=state$uvec, replace=TRUE)
      }
      naccepted <- 0
      llike.old <- llikef(state, datlist, fpar, model)
      if(ppar$verbose) cat("vvvec:")
      for(k in ksamplevec) {
        v.old <- state$vvec[k]
        v.new <- expit( rnorm(1,logit(v.old),ppar$sd.logit.v) )
        state$vvec[k] <- v.new
        state <- update.vw(state, fpar)
        llike.new <- llikef(state, datlist, fpar, model)
        
        log.r <- (state$alpha)*log((1-v.new)/(1-v.old))
        if(ppar$ksweep) {
          log.r <- log.r + log(v.new/v.old)
        } else {
          log.r <- log.r + 2*log(v.new/v.old)  
        }
        log.r <- log.r + llike.new - llike.old
        if(ppar$verbose) cat(k)
        if(runif(1)<exp(log.r)) {
          # accepted
          naccepted <- naccepted + 1
          llike.old <- llike.new
          if(ppar$verbose) cat("+")
        } else {
          # reject
          state$vvec[k] <- v.old
          state <- update.vw(state, fpar)
          if(ppar$verbose) cat("-")
        }
      }
      if(ppar$verbose) cat("\n")
      state$accepted["vvec"] <- naccepted
      state$llike <- llike.old
    }
    
    # update thetavec ##!!==
    if(ppar$update["thetavec"]) {
      # log Normal proposal
      if(ppar$ksweep) {
        # sweep update
        ksamplevec <- 1:fpar$kmax
      } else {
        # random support point update
        ksamplevec <- sample(fpar$kmax, ppar$ksim,
                             prob=state$uvec, replace=TRUE)
      }
      naccepted <- 0
      llike.old <- llikef(state, datlist, fpar, model)
      for(k in ksamplevec) {
        theta.old <- state$thetavec[k]
        theta.new <- exp( rnorm(1,log(theta.old),ppar$sd.log.theta) )
        state$thetavec[k] <- theta.new
        llike.new <- llikef(state, datlist, fpar, model)
        log.r <- log(theta.new/theta.old) - state$phi*(theta.new-theta.old)
        log.r <- log.r + llike.new - llike.old
        if(runif(1)<exp(log.r)) {
          # accepted
          naccepted <- naccepted + 1
          llike.old <- llike.new
        } else {
          # reject
          state$thetavec[k] <- theta.old
        }
      }
      state$accepted["thetavec"] <- naccepted
      state$llike <- llike.old
    }
  } else {
    stop("Specified model has not been implemented")
  }
  
  # update other constants
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  
  return(state)
}
########################################################################
