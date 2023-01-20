# General Fragments for updating

#' eta V1 (IFR/DFR)
#'
#' @export
update.eta.v1 <- function(state, datlist, fpar, ppar, model,
                          nm=c(eta="eta",
                               gamma="gamma",
                               thetavec="thetavec",
                               vvec="vvec",
                               alpha="alpha",
                               beta="beta",
                               phi="phi",
                               uvec="uvec",
                               wvec="wvec",
                               lambda0="lambda0",
                               thetaswap="thetaswap")) {

  state.old <- state
  eta.old <- state.old[[nm["eta"]]]
  eta.new <- exp( rnorm(1, log(eta.old), ppar$sd.log.eta) )
  state[[nm["eta"]]] <- eta.new
  state[[nm["lambda0"]]] <- state[[nm["gamma"]]]*state[[nm["eta"]]]
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- log.r + log(eta.new/eta.old) - fpar$nu*(eta.new-eta.old)
  if(ppar$verbose) {
    cat(sprintf("eta: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                eta.old, eta.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
    cat(sprintf("eta: %g->%g: logr=%g\n",
                eta.old, eta.new, log.r))
    browser()
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["eta"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["eta"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

#' gamma V1 (Gibbs) (IFR/DFR)
#'
#' @export
update.gamma.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(eta="eta",
                                 gamma="gamma",
                                 thetavec="thetavec",
                                 vvec="vvec",
                                 alpha="alpha",
                                 beta="beta",
                                 phi="phi",
                                 uvec="uvec",
                                 wvec="wvec",
                                 lambda0="lambda0",
                                 thetaswap="thetaswap")) {

  llike.old <- state$llike
  gamma.old <- state$gamma
  temp.model.list <- list(model=model, kmax=kmax, 
                          lambda0=state[[nm["eta"]]], ## Note: unscaled
                          thetavec=state[[nm["thetavec"]]],
                          wvec=state[[nm["uvec"]]])   ## Note: unscaled
  cc <- sum(int.lambda.func(datlist$tvec, model.list=temp.model.list))
  g1star <- state[[nm["alpha"]]]+datlist$n0
  g2star <- state[[nm["beta"]]] + cc
  gamma.new <- rgamma(1, g1star, g2star)
  state[[nm["gamma"]]] <- gamma.new
  state[[nm["lambda0"]]] <- state[[nm["gamma"]]]*state[[nm["eta"]]]
  state[[nm["wvec"]]] <- state[[nm["gamma"]]]*state[[nm["uvec"]]]
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  state$accepted[nm["gamma"]] <- 1
  if(ppar$verbose) {
    cat(sprintf("gamma: (%g;%g) %g->%g: logr=Gibbs\n",
                llike.old, state$llike,
                gamma.old, gamma.new))
    cat("+")
  }
  
  return(state)  
}

#' thetavec V1 (IFR/DFR)
#'
#' @export
update.thetavec.v1 <- function(state, datlist, fpar, ppar, model,
                               nm=c(eta="eta",
                                    gamma="gamma",
                                    thetavec="thetavec",
                                    vvec="vvec",
                                    alpha="alpha",
                                    beta="beta",
                                    phi="phi",
                                    uvec="uvec",
                                    wvec="wvec",
                                    lambda0="lambda0",
                                    thetaswap="thetaswap")) {

  sweep <- as.logical(rbinom(1,1,ppar$psweep))
  if(ppar$verbose && sweep) cat("(sweep):")
  if(sweep || ppar$ksweep) {
    ksamplevec <- 1:fpar$kmax
  } else {
    ksamplevec <- sample(fpar$kmax, ppar$ksim,
                         prob=state$uvec, replace=TRUE)
  }
  naccepted <- 0
  for(k in ksamplevec) {
    state.old <- state
    theta.old <- state[[nm["thetavec"]]][k]
    theta.new <- exp( rnorm(1,log(theta.old),ppar$sd.log.theta) )
    state[[nm["thetavec"]]][k] <- theta.new
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    log.r <- state$llike - state.old$llike
    log.r <- (log.r - state[[nm["phi"]]]*(theta.new-theta.old)
              + log(theta.new/theta.old))
    if(ppar$verbose) {
      cat(sprintf("thetavec[%d]: (%g;%g) %g->%g: logr=%g\n",
                  k, state.old$llike, state$llike,
                  theta.old, theta.new, log.r))
    }
    if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
      cat(sprintf("thetavec[%d]: %g->%g: logr=%g\n",
                  k, theta.old, theta.new, log.r))
      browser()
    }
    if(runif(1)<exp(log.r)) {
      # accept
      naccepted <- naccepted + 1
      if(ppar$verbose) cat("+")
    } else {
      # reject
      state <- state.old
      if(ppar$verbose) cat("-")
    }
  }
  state$accepted[nm["thetavec"]] <- naccepted/length(ksamplevec)

  return(state)
}

#' vvec V1 (IFR/DFR)
#'
#' @export
update.vvec.v1 <- function(state, datlist, fpar, ppar, model,
                           nm=c(eta="eta",
                                gamma="gamma",
                                thetavec="thetavec",
                                vvec="vvec",
                                alpha="alpha",
                                beta="beta",
                                phi="phi",
                                uvec="uvec",
                                wvec="wvec",
                                lambda0="lambda0",
                                thetaswap="thetaswap")) {
    
  sweep <- as.logical(rbinom(1,1,ppar$psweep))
  if(ppar$verbose && sweep) cat("(sweep):")
  if(sweep || ppar$ksweep) {
    ksamplevec <- 1:(fpar$kmax-1)
  } else {
    ksamplevec <- sample(fpar$kmax-1, ppar$ksim,
                         prob=state[[nm["uvec"]]][-fpar$kmax], replace=TRUE)
  }
  naccepted <- 0
  for(k in ksamplevec) {
    state.old <- state
    ukmax.old <- state[[nm["uvec"]]][fpar$kmax]
    v.old <- state[[nm["vvec"]]][k]
    v.new <- expit( rnorm(1,logit(v.old),ppar$sd.logit.v) )
    state[[nm["vvec"]]][k] <- v.new
    cp <- cumprod(1-state[[nm["vvec"]]][-fpar$kmax])
    state[[nm["uvec"]]] <- state[[nm["vvec"]]]*c(1,cp)
    state[[nm["uvec"]]][fpar$kmax] <- cp[fpar$kmax-1]
    #state[[nm["uvec"]]][fpar$kmax] <- 1-sum(state[[nm["uvec"]]][-fpar$kmax])
    state[[nm["wvec"]]] <- state[[nm["gamma"]]] * state[[nm["uvec"]]]
    ukmax.new <- state[[nm["uvec"]]][fpar$kmax]
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    log.r <- state$llike - state.old$llike
    if(sweep || ppar$ksweep) {
      # updating in sequence
      log.r <- (log.r + state[[nm["alpha"]]]*log((1-v.new)/(1-v.old)) 
                + log(v.new/v.old)
      )
    } else {
      # updating at random
      log.r <- (log.r + state[[nm["alpha"]]]*log((1-v.new)/(1-v.old)) 
                + 2*log(v.new/v.old)
                + log((1-ukmax.old)/(1-ukmax.new))
      )
    }
    if(ppar$verbose) {
      cat(sprintf("vvec[%d]: (%g;%g) %g->%g: logr=%g\n",
                  k, state.old$llike, state$llike,
                  v.old, v.new, log.r))
    }
    if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
      cat(sprintf("vvec[%d]: %g->%g: logr=%g\n",
                  k, v.old, v.new, log.r))
      browser()
    }
    if(runif(1)<exp(log.r)) {
      # accept
      naccepted <- naccepted + 1
      if(ppar$verbose) cat("+")
    } else {
      # reject
      state <- state.old
      if(ppar$verbose) cat("-")
    }
  }
  state$accepted[nm["vvec"]] <- naccepted/length(ksamplevec)
  
  return(state)
}

#' alpha V1 (IFR/DFR)
#'
#' @export
update.alpha.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(eta="eta",
                                 gamma="gamma",
                                 thetavec="thetavec",
                                 vvec="vvec",
                                 alpha="alpha",
                                 beta="beta",
                                 phi="phi",
                                 uvec="uvec",
                                 wvec="wvec",
                                 lambda0="lambda0",
                                 thetaswap="thetaswap")) {

  state.old <- state
  alpha.old <- state.old[[nm["alpha"]]]
  alpha.new <- exp( rnorm(1, log(alpha.old), ppar$sd.log.alpha) )
  state[[nm["alpha"]]] <- alpha.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.ukmax <- log(state[[nm["uvec"]]][fpar$kmax]) 
  #log.ukmax - sum(log(1-state[[nm["vvec"]]][-fpar$kmax])) 
  a1star <- fpar$a1 + fpar$kmax - 1
  a2star <- fpar$a2 - log.ukmax - log(state[[nm["beta"]]]*state[[nm["gamma"]]])
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + a1star*log(alpha.new/alpha.old) 
             - a2star*(alpha.new-alpha.old) 
             + lgamma(alpha.old) - lgamma(alpha.new) )
  if(ppar$verbose) {
    cat(sprintf("alpha: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                alpha.old, alpha.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
    cat(sprintf("alpha: %g->%g: logr=%g\n",
                alpha.old, alpha.new, log.r))
    browser()
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["alpha"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["alpha"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

#' beta V1 (IFR/DFR)
#'
#' @export
update.beta.v1 <- function(state, datlist, fpar, ppar, model,
                           nm=c(eta="eta",
                                gamma="gamma",
                                thetavec="thetavec",
                                vvec="vvec",
                                alpha="alpha",
                                beta="beta",
                                phi="phi",
                                uvec="uvec",
                                wvec="wvec",
                                lambda0="lambda0",
                                thetaswap="thetaswap")) {

  beta.old <- state[[nm["beta"]]]
  llike.old <- state$llike
  b1star <- fpar$b1 + state[[nm["alpha"]]]
  b2star <- fpar$b2 + state[[nm["gamma"]]]
  beta.new <- rgamma(1, b1star, b2star)
  state[[nm["beta"]]] <- beta.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  state$accepted[nm["beta"]] <- 1
  if(ppar$verbose) {
    cat(sprintf("beta: (%g;%g) %g->%g: logr=Gibbs\n",
                llike.old, state$llike,
                beta.old, beta.new))
    cat("+")
  }
  
  return(state)
}

#' phi V1 (IFR/DFR)
#'
#' @export
update.phi.v1 <- function(state, datlist, fpar, ppar, model,
                          nm=c(eta="eta",
                               gamma="gamma",
                               thetavec="thetavec",
                               vvec="vvec",
                               alpha="alpha",
                               beta="beta",
                               phi="phi",
                               uvec="uvec",
                               wvec="wvec",
                               lambda0="lambda0",
                               thetaswap="thetaswap",
                               thetaswap="thetaswap")) {
  
  phi.old <- state[[nm["phi"]]]
  llike.old <- state$llike
  f1star <- fpar$f1 + fpar$kmax
  f2star <- fpar$f2 + sum(state[[nm["thetavec"]]])
  phi.new <- rgamma(1, f1star, f2star)
  state[[nm["phi"]]] <- phi.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  state$accepted[nm["phi"]] <- 1
  if(ppar$verbose) {
    cat(sprintf("phi: (%g;%g) %g->%g: logr=Gibbs\n",
                llike.old, state$llike,
                phi.old, phi.new))
    cat("+")
  }
  
  return(state)
}

#' wvec V1 (IFR/DFR)
#'
#' @export
update.wvec.v1 <- function(state, datlist, fpar, ppar, model,
                           nm=c(eta="eta",
                                gamma="gamma",
                                thetavec="thetavec",
                                vvec="vvec",
                                alpha="alpha",
                                beta="beta",
                                phi="phi",
                                uvec="uvec",
                                wvec="wvec",
                                lambda0="lambda0",
                                thetaswap="thetaswap",
                                thetaswap="thetaswap")) {

  # Do not update wvec[kmax]
  #sweep <- as.logical(rbinom(1,1,ppar$psweep))
  sweep <- FALSE # do not do a global update of every wvec[k]
  if(ppar$verbose && sweep) cat("(sweep):")
  if(sweep || ppar$ksweep) {
    ksamplevec <- 1:(fpar$kmax-1)
  } else {
    ksamplevec <- sample(fpar$kmax-1, ppar$ksim,
                         prob=state[[nm["uvec"]]][-fpar$kmax], replace=TRUE)
  }
  naccepted <- 0
  for(k in ksamplevec) {
    state.old <- state
    vvec.old <- state[[nm["vvec"]]]
    gamma.old <- state[[nm["gamma"]]]
    wkmax.old <- state[[nm["wvec"]]][fpar$kmax]
    w.old <- state[[nm["wvec"]]][k]
    w.new <- exp( rnorm(1,log(w.old),ppar$sd.log.w) )
    state[[nm["wvec"]]][k] <- w.new
    gamma.new <- sum(state[[nm["wvec"]]])
    uvec.new <- state[[nm["wvec"]]]/gamma.new
    uvec.new <- pmax(.Machine$double.neg.eps, uvec.new) # needed to stabilise
    uvec.new <- uvec.new/sum(uvec.new)
    
    # v values only change for indices up to and including k
    vvec.new <- vvec.old
    vvec.new[1:k] <- wvec[1:k]/(c(1, 1-cumsum(wvec[-kmax]))[1:k])
    state[[nm["gamma"]]] <- gamma.new
    state[[nm["vvec"]]] <- vvec.new
    state[[nm["uvec"]]] <- uvec.new
    state[[nm["wvec"]]] <- state[[nm["gamma"]]]*state[[nm["uvec"]]]
    state[[nm["lambda0"]]] <- state[[nm["gamma"]]]*state[[nm["eta"]]]
    wkmax <- state[[nm["wvec"]]][fpar$kmax]
    
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    
    log.r <- state$llike - state.old$llike
    log.r <- (log.r - state[[nm["beta"]]]*(gamma.new-gamma.old)
              + sum(log(vvec.new[1:k]/vvec.old[1:k])))
    
    if(sweep || ppar$ksweep) {
      # updating sequentially
      log.r <- log.r
    } else {
      # updating at random
      log.r <- (log.r + log((w.new/gamma.new)/(w.old/gamma.old)))
    }
    
    if(ppar$verbose) {
      cat(sprintf("wvec[%d]: (%g;%g) %g->%g: logr=%g\n",
                  k, state.old$llike, state$llike,
                  w.old, w.new, log.r))
    }
    if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
      cat(sprintf("wvec[%d]: %g->%g: logr=%g\n",
                  k, w.old, w.new, log.r))
      browser()
    }
    if(runif(1)<exp(log.r)) {
      # accept
      naccepted <- naccepted + 1
      if(ppar$verbose) cat("+")
    } else {
      # reject
      state <- state.old
      if(ppar$verbose) cat("-")
    }
  }
  state$accepted[nm["wvec"]] <- naccepted/length(ksamplevec)
  
  return(state) 
}
    
#' thetaswap V1 (IFR/DFR)
#'
#' @export
update.thetaswap.v1 <- function(state, datlist, fpar, ppar, model,
                                nm=c(eta="eta",
                                     gamma="gamma",
                                     thetavec="thetavec",
                                     vvec="vvec",
                                     alpha="alpha",
                                     beta="beta",
                                     phi="phi",
                                     uvec="uvec",
                                     wvec="wvec",
                                     lambda0="lambda0",
                                     thetaswap="thetaswap",
                                     thetaswap="thetaswap")) {

  ksamplevec1 <- sample(fpar$kmax, ppar$kswap,
                        prob=state[[nm["uvec"]]], replace=TRUE)
  ksamplevec2 <- sample(fpar$kmax, ppar$kswap, replace=TRUE)
  naccepted <- 0
  for(i in 1:ppar$kswap) {
    k1 <- ksamplevec1[i]
    k2 <- ksamplevec2[i]
    state.old <- state
    state[[nm["thetavec"]]][c(k1,k2)] <- state.old[[nm["thetavec"]]][c(k2,k1)]
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    log.r <- state$llike - state.old$llike
    log.r <- (log.r + log(state[[nm["wvec"]]][k2]/state[[nm["wvec"]]][k1]))
    if(ppar$verbose) {
      cat(sprintf("thetaswap[%d,%d]: (%g;%g) %g->%g: logr=%g\n",
                  k1, k2, state.old$llike, state$llike,
                  state.old[[nm["thetavec"]]][k1], state.old[[nm["thetavec"]]][k2], log.r))
    }
    if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
      cat(sprintf("thetaswap[%d,%d]: %g->%g: logr=%g\n",
                  k1, k2, state.old[[nm["thetavec"]]][k1], state.old[[nm["thetavec"]]][k2], log.r))
      browser()
    }
    if(runif(1)<exp(log.r)) {
      # accept
      naccepted <- naccepted + 1
      if(ppar$verbose) cat("+")
    } else {
      # reject
      state <- state.old
      if(ppar$verbose) cat("-")
    }
  }
  state$accepted[nm["thetaswap"]] <- naccepted/ppar$kswap
  
  return(state)
}
