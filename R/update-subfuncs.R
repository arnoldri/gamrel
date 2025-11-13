# General Fragments for updating

#' lambda0 V0 (Gibbs) (CON)
#'
#' @export
update.lambda0.v0 <- function(state, datlist, fpar, ppar, model,
                              nm=c(lambda0="lambda0")) {
  last.seed <- .Random.seed
  llike.old <- state$llike
  temp.model.list <- list(model=model, lambda0=1) 
  cc <- sum(datlist$tvec)
  
  lambda0.old <- state[[nm["lambda0"]]]
  s1star <- fpar$s1 + datlist$nobs
  s2star <- fpar$s2 + cc
  lambda0.new <- rgamma(1, s1star, s2star)
  state[[nm["lambda0"]]] <- lambda0.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  if(ppar$verbose) {
    cat(sprintf("a: (%g;%g) %g->%g: logr=Gibbs\n",
                llike.old, state$llike,
                lambda0.old, lambda0.new))
  }
  state$accepted[nm["lambda0"]] <- 1
  
  return(state)
}

#' alpha V1 (IFR/DFR)
#'
#' @export
update.alpha.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(lambda0="lambda0",
                                 gamma="gamma",
                                 thetavec="thetavec",
                                 vvec="vvec",
                                 alpha="alpha",
                                 beta="beta",
                                 nu="nu",
                                 phi="phi",
                                 uvec="uvec",
                                 wvec="wvec",
                                 lambda0="lambda0")) {
  last.seed <- .Random.seed
  state.old <- state
  alpha.old <- state.old[[nm["alpha"]]]
  alpha.new <- exp( rnorm(1, log(alpha.old), ppar$sd.log.alpha) )
  state[[nm["alpha"]]] <- alpha.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.ukmax <- log(state[[nm["uvec"]]][fpar$kmax]) 
  #log.ukmax - sum(log(1-state[[nm["vvec"]]][-fpar$kmax])) 
  a1star <- fpar$a1 + fpar$kmax - 1
  a2star <- fpar$a2 - log.ukmax - log(state[[nm["gamma"]]]*state[[nm["beta"]]])
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + a1star*log(alpha.new/alpha.old) 
             - a2star*(alpha.new-alpha.old) 
             + lgamma(alpha.old) - lgamma(alpha.new) )
  if(ppar$verbose) {
    cat(sprintf("alpha: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                alpha.old, alpha.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("alpha: %g->%g: logr=%g\n",
                alpha.old, alpha.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["alpha"]
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata"))
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
                           nm=c(lambda0="lambda0",
                                gamma="gamma",
                                thetavec="thetavec",
                                vvec="vvec",
                                alpha="alpha",
                                beta="beta",
                                nu="nu",
                                phi="phi",
                                uvec="uvec",
                                wvec="wvec",
                                lambda0="lambda0")) {
  last.seed <- .Random.seed
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

#' nu V1 (IFR/DFR)
#'
#' @export
update.nu.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(lambda0="lambda0",
                                 gamma="gamma",
                                 thetavec="thetavec",
                                 vvec="vvec",
                                 alpha="alpha",
                                 beta="beta",
                                 nu="nu",
                                 phi="phi",
                                 uvec="uvec",
                                 wvec="wvec",
                                 lambda0="lambda0")) {
  last.seed <- .Random.seed
  state.old <- state
  nu.old <- state.old[[nm["nu"]]]
  nu.new <- exp( rnorm(1, log(nu.old), ppar$sd.log.nu) )
  state[[nm["nu"]]] <- nu.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  g1star <- fpar$g1
  g2star <- fpar$g2 - fpar$kmax*log(state[[nm["phi"]]]) - sum(log(state[[nm["thetavec"]]]))
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + g1star*log(nu.new/nu.old) 
             - g2star*(nu.new-nu.old) 
             + fpar$kmax*(lgamma(nu.old) - lgamma(nu.new)) )
  if(ppar$verbose) {
    cat(sprintf("nu: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                nu.old, nu.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("nu: %g->%g: logr=%g\n",
                nu.old, nu.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["nu"]
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata"))
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["nu"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["nu"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}


#' phi V1 (IFR/DFR)
#'
#' @export
update.phi.v1 <- function(state, datlist, fpar, ppar, model,
                          nm=c(lambda0="lambda0",
                               gamma="gamma",
                               thetavec="thetavec",
                               vvec="vvec",
                               alpha="alpha",
                               beta="beta",
                               nu="nu",
                               phi="phi",
                               uvec="uvec",
                               wvec="wvec",
                               lambda0="lambda0",
                               f1="f1", f2="f2")) {
  last.seed <- .Random.seed
  phi.old <- state[[nm["phi"]]]
  llike.old <- state$llike
  f1star <- fpar[[nm["f1"]]] + fpar$kmax*state[[nm["nu"]]]
  f2star <- fpar[[nm["f2"]]] + sum(state[[nm["thetavec"]]])
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

#' lambda0 V1 (IFR/DFR)
#'
#' @export
update.lambda0.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(lambda0="lambda0",
                                 gamma="gamma",
                                 thetavec="thetavec",
                                 vvec="vvec",
                                 alpha="alpha",
                                 beta="beta",
                                 nu="nu",
                                 phi="phi",
                                 uvec="uvec",
                                 wvec="wvec",
                                 lambda0="lambda0")) {
  last.seed <- .Random.seed
  state.old <- state
  lambda0.old <- state.old[[nm["lambda0"]]]
  lambda0.new <- exp( rnorm(1, log(lambda0.old), ppar$sd.log.lambda0) )
  state[[nm["lambda0"]]] <- lambda0.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  s1star <- fpar$s1
  s2star <- fpar$s2
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + s1star*log(lambda0.new/lambda0.old) 
             - s2star*(lambda0.new-lambda0.old)  )
  if(ppar$verbose) {
    cat(sprintf("lambda0: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                lambda0.old, lambda0.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("lambda0: %g->%g: logr=%g\n",
                lambda0.old, lambda0.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["lambda0"]
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata"))
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["lambda0"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["lambda0"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

#' thetavec V1 (IFR/DFR)
#'
#' @export
update.thetavec.v1 <- function(state, datlist, fpar, ppar, model,
                               nm=c(lambda0="lambda0",
                                    gamma="gamma",
                                    thetavec="thetavec",
                                    vvec="vvec",
                                    alpha="alpha",
                                    beta="beta",
                                    nu="nu",
                                    phi="phi",
                                    uvec="uvec",
                                    wvec="wvec",
                                    lambda0="lambda0")) {
  last.seed <- .Random.seed
  equal <- as.logical(rbinom(1,1,ppar$pequal))
  if(ppar$verbose && equal) cat("(equal):")
  if(equal) {
    k <- sample(fpar$kmax, 1)
  } else {
    k <- sample(fpar$kmax, 1, prob=state$uvec)
  }
  state.old <- state
  theta.old <- state[[nm["thetavec"]]][k]
  theta.new <- exp( rnorm(1,log(theta.old),ppar$sd.log.theta) )
  state[[nm["thetavec"]]][k] <- theta.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- state$llike - state.old$llike
  log.r <- (log.r - state[[nm["phi"]]]*(theta.new-theta.old)
                  + state[[nm["nu"]]]*log(theta.new/theta.old))
  if(ppar$verbose) {
    cat(sprintf("thetavec[%d]: (%g;%g) %g->%g: logr=%g\n",
                k, state.old$llike, state$llike,
                theta.old, theta.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("thetavec[%d]: %g->%g: logr=%g\n",
                k, theta.old, theta.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["thetavec"]
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed, k,
         file=paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata"))
  }
    
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["thetavec"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["thetavec"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

#' gamma V1 (IFR/DFR)
#'
#' @export
update.gamma.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(lambda0="lambda0",
                                 gamma="gamma",
                                 thetavec="thetavec",
                                 vvec="vvec",
                                 alpha="alpha",
                                 beta="beta",
                                 nu="nu",
                                 phi="phi",
                                 uvec="uvec",
                                 wvec="wvec",
                                 lambda0="lambda0")) {
  last.seed <- .Random.seed
  state.old <- state
  gamma.old <- state.old[[nm["gamma"]]]
  gamma.new <- exp( rnorm(1, log(gamma.old), ppar$sd.log.gamma) )
  state[[nm["gamma"]]] <- gamma.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + state[[nm["alpha"]]]*log(gamma.new/gamma.old) 
             - state[[nm["beta"]]]*(gamma.new-gamma.old)  )
  if(ppar$verbose) {
    cat(sprintf("gamma: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                gamma.old, gamma.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("gamma: %g->%g: logr=%g\n",
                gamma.old, gamma.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["gamma"]
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata"))
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["gamma"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["gamma"]] <- 0
    if(ppar$verbose) cat("-")
  }

  return(state)  
}

#' vvec V1 (IFR/DFR/LWB/SBT/MBT/LCV)
#'
#' @export
update.vvec.v1 <- function(state, datlist, fpar, ppar, model,
                           nm=c(lambda0="lambda0",
                                gamma="gamma",
                                thetavec="thetavec",
                                vvec="vvec",
                                alpha="alpha",
                                beta="beta",
                                nu="nu",
                                phi="phi",
                                uvec="uvec",
                                wvec="wvec",
                                lambda0="lambda0")) {
  last.seed <- .Random.seed  
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
    if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
      cat(sprintf("vvec[%d]: %g->%g: logr=%g\n",
                  k, v.old, v.new, log.r))
      if(ppar$interactive) browser()
      update.name <- nm["vvec"]
      fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
      fname <- gsub(":","-",fname,fixed=TRUE)
      save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
           file=fname)
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

<<<<<<< HEAD
=======
#' alpha V1 (IFR/DFR/LWB/SBT/MBT/LCV)
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
  last.seed <- .Random.seed
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
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("alpha: %g->%g: logr=%g\n",
                alpha.old, alpha.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["alpha"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
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

#' beta V1 (IFR/DFR/LWB/SBT/MBT/LCV)
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
  last.seed <- .Random.seed
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

#' phi V1 (IFR/DFR/LWB/SBT/MBT/LCV)
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
                               f1="f1", f2="f2")) {
  last.seed <- .Random.seed
  phi.old <- state[[nm["phi"]]]
  llike.old <- state$llike
  f1star <- fpar[[nm["f1"]]] + fpar$kmax
  f2star <- fpar[[nm["f2"]]] + sum(state[[nm["thetavec"]]])
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

>>>>>>> 64d170090dc20b09a0472af762c925c55df1618e
#' wvec V1 (IFR/DFR/LWB/SBT/MBT/LCV)
#'
#' @export
update.wvec.v1 <- function(state, datlist, fpar, ppar, model,
                           nm=c(lambda0="lambda0",
                                gamma="gamma",
                                thetavec="thetavec",
                                vvec="vvec",
                                alpha="alpha",
                                beta="beta",
                                nu="nu",
                                phi="phi",
                                uvec="uvec",
                                wvec="wvec",
                                lambda0="lambda0")) {
  last.seed <- .Random.seed
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
    ##uvec.new <- pmax(.Machine$double.neg.eps, uvec.new) # needed to stabilise
    uvec.new <- uvec.new/sum(uvec.new)
    
    # v values only change for indices up to and including k
    # but recalculate all due to potential numerical rounding errors, and avoid 1.0 also
    vvec.new <- pmin(1-fpar$epsilon, uvec.new/rev(cumsum(rev(uvec.new))))
    # recalculate uvec.new to ensure consistency
    uvec.new <- c(vvec.new[-fpar$kmax],1)*c(1,cumprod(1-vvec.new[-fpar$kmax]))
    state[[nm["gamma"]]] <- gamma.new
    state[[nm["vvec"]]] <- vvec.new
    state[[nm["uvec"]]] <- uvec.new
    state[[nm["wvec"]]] <- state[[nm["gamma"]]]*state[[nm["uvec"]]]
    if(model %in% c("IFR","DFR","LWB","SBT","MBT") && !is.na(nm["lambda0"])) {
       state[[nm["lambda0"]]] <- state[[nm["gamma"]]]*state[[nm["eta"]]]
    } # no need to do anything in the LCV case where gamma doesn't scale other parameters
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
      log.r <- (log.r + log((w.new/(gamma.new-wkmax))/(w.old/(gamma.old-wkmax.old))))
    }
    
    if(ppar$verbose) {
      cat(sprintf("wvec[%d]: (%g;%g) %g->%g: logr=%g\n",
                  k, state.old$llike, state$llike,
                  w.old, w.new, log.r))
    }
    if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==} || (log.r>0.9 && ppar$verbose)) { ##!!==
      cat(sprintf("wvec[%d]: %g->%g: logr=%g\n",
                  k, w.old, w.new, log.r))
      cat("model:"); cat(model); cat("\n")
      cat("names(state):"); cat(names(state)); cat("\n")
      cat("nm:\n"); print(nm)
      cat("k, llike, llike.old, gamma.new, gamma.new, sum(vv/vv), w.new, w.old, wkmax, wkmax.old)\n")
      cat(c(k, state$llike, state.old$llike, gamma.new, gamma.old, sum(log(vvec.new[1:k]/vvec.old[1:k])),
            w.new, w.old, wkmax, wkmax.old)); cat("\n")
      cat("logr: LR, gamma, sumvv, logw/wold, log.r:\n")
      cat(c(state$llike - state.old$llike,
            -state[[nm["beta"]]]*(gamma.new-gamma.old),
            sum(log(vvec.new[1:k]/vvec.old[1:k])),
            log((w.new/(gamma.new-wkmax))/(w.old/(gamma.old-wkmax.old))),
            log.r
            )); cat("\n")
      if(ppar$interactive) browser()
      update.name <- nm["wvec"]
      fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
      fname <- gsub(":","-",fname,fixed=TRUE)
      save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
           file=fname)
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
    
#' a V1 (LWB)
#'
#' @export
update.a.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(lambda0="lambda0",
                                 gamma="gamma",
                                 thetavec="thetavec",
                                 vvec="vvec",
                                 alpha="alpha",
                                 beta="beta",
                                 nu="nu",
                                 phi="phi",
                                 uvec="uvec",
                                 wvec="wvec",
                                 lambda0="lambda0")) {
  last.seed <- .Random.seed
  state.old <- state
  a.old <- state.old[[nm["a"]]]
  a.new <- exp( rnorm(1, log(a.old), ppar$sd.log.a) )
  state[[nm["a"]]] <- a.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + fpar$c1*log(a.new/a.old) 
             - fpar$c2*(a.new-a.old) )
  if(ppar$verbose) {
    cat(sprintf("a: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                a.old, a.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("a: %g->%g: logr=%g\n",
                a.old, a.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["a"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["a"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["a"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

#' lambda0 V2 (Gibbs) (LCV)
#'
#' @export
update.lambda0.v2 <- function(state, datlist, fpar, ppar, model,
                              nm=c(lambda0="lambda0",
                                   w0="w0",
                                   gamma="gamma",
                                   thetavec="thetavec",
                                   vvec="vvec",
                                   alpha="alpha",
                                   beta="beta",
                                   nu="nu",
                                   phi="phi",
                                   uvec="uvec",
                                   wvec="wvec")) {
  last.seed <- .Random.seed
  llike.old <- state$llike
  temp.model.list <- list(model=model, kmax=fpar$kmax, 
                          lambda0=1, ## Note: unscaled
                          w0=state[[nm["w0"]]],
                          thetavec=state[[nm["thetavec"]]],
                          wvec=state[[nm["wvec"]]]) 
  cc <- sum(int.lambda.func(datlist$tvec, model.list=temp.model.list))
  
  lambda0.old <- state[[nm["lambda0"]]]
  s1star <- fpar$s1 + datlist$nobs
  s2star <- fpar$s2 + cc
  lambda0.new <- rgamma(1, s1star, s2star)
  state[[nm["lambda0"]]] <- lambda0.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  if(ppar$verbose) {
    cat(sprintf("a: (%g;%g) %g->%g: logr=Gibbs\n",
                llike.old, state$llike,
                lambda0.old, lambda0.new))
  }
  state$accepted[nm["lambda0"]] <- 1

  return(state)
}

#' w0 V1 (LCV)
#'
#' @export
update.w0.v1 <- function(state, datlist, fpar, ppar, model,
                              nm=c(lambda0="lambda0",
                                   w0="w0",
                                   gamma="gamma",
                                   thetavec="thetavec",
                                   vvec="vvec",
                                   alpha="alpha",
                                   beta="beta",
                                   nu="nu",
                                   phi="phi",
                                   uvec="uvec",
                                   wvec="wvec")) {
  last.seed <- .Random.seed
  state.old <- state
  w0.old <- state.old[[nm["w0"]]]
  w0.new <- rnorm(1, w0.old, ppar$sd.w0) 
  state[[nm["w0"]]] <- w0.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r - 0.5*(w0.new^2-w0.old^2)/(fpar$sigmap.w0^2) ) 
  if(ppar$verbose) {
    cat(sprintf("w0: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                w0.old, w0.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("w0: %g->%g: logr=%g\n",
                w0.old, w0.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["w0"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["w0"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["w0"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

#' gamma V2 (MH) (SBT,MBT,LCV)
#'
#' @export
update.gamma.v2 <- function(state, datlist, fpar, ppar, model,
                              nm=c(lambda0="lambda0", # NA in first component of SBT
                                   w0="w0", #ignored in SBT, MBT
                                   gamma="gamma",
                                   thetavec="thetavec",
                                   vvec="vvec",
                                   alpha="alpha",
                                   beta="beta",
                                   nu="nu",
                                   phi="phi",
                                   uvec="uvec",
                                   wvec="wvec")) {
  last.seed <- .Random.seed
  state.old <- state
  gamma.old <- state.old[[nm["gamma"]]]
  gamma.new <- exp( rnorm(1, log(gamma.old), ppar$sd.log.gamma) )
  state[[nm["gamma"]]] <- gamma.new
  if(model!="LCV" && !is.na(nm[["lambda0"]])) {  
    # LCV does not scale lambda0 by gamma
    # lambda0 will be missing for component 1 of SBT
    state[[nm["lambda0"]]] <- state[[nm["gamma"]]]*state[[nm["eta"]]]
  }
  state[[nm["wvec"]]] <- state[[nm["gamma"]]]*state[[nm["uvec"]]]
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + state[[nm["alpha"]]]*log(gamma.new/gamma.old) 
             - state[[nm["beta"]]]*(gamma.new-gamma.old) )
  if(ppar$verbose) {
    cat(sprintf("gamma: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                gamma.old, gamma.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("gamma: %g->%g: logr=%g\n",
                gamma.old, gamma.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["gamma"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["gamma"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["gamma"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

#' pival V1 (MBT)
#'
#' @export
update.pival.v1 <- function(state, datlist, fpar, ppar, model,
                            nm=c(pival="pival")) {
  last.seed <- .Random.seed
  state.old <- state
  pival.old <- state.old[[nm["pival"]]]
  pival.new <- expit( rnorm(1, logit(pival.old), ppar$sd.logit.pival) ) 
  state[[nm["pival"]]] <- pival.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + log(pival.new*(1-pival.new)/(pival.old*(1-pival.old))) ) 
  if(ppar$verbose) {
    cat(sprintf("pival: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                pival.old, pival.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("pival: %g->%g: logr=%g\n",
                pival.old, pival.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["pival"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["pival"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["pival"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}

## MEW functions

#' alpha (MH) MEW
#'
#' @export
update.alpha.mew <- function(state, datlist, fpar, ppar, model,
                              nm=c(alpha="alpha",
                                   beta="beta",
                                   mu="mu",
                                   nu="nu")) {
  last.seed <- .Random.seed
  state.old <- state
  alpha.old <- state.old[[nm["alpha"]]]
  alpha.new <- exp( rnorm(1, log(alpha.old), ppar$sd.log.alpha) )
  state[[nm["alpha"]]] <- alpha.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + fpar$a1*log(alpha.new/alpha.old) 
                   - fpar$a2*(alpha.new-alpha.old) )
  if(ppar$verbose) {
    cat(sprintf("alpha: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                alpha.old, alpha.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("alpha: %g->%g: logr=%g\n",
                alpha.old, alpha.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["alpha"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
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

#' beta (MH) MEW
#'
#' @export
update.beta.mew <- function(state, datlist, fpar, ppar, model,
                            nm=c(alpha="alpha",
                                 beta="beta",
                                 mu="mu",
                                 nu="nu")) {
  last.seed <- .Random.seed
  state.old <- state
  beta.old <- state.old[[nm["beta"]]]
  beta.new <- exp( rnorm(1, log(beta.old), ppar$sd.log.beta) )
  state[[nm["beta"]]] <- beta.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + fpar$b1*log(beta.new/beta.old) 
                   - fpar$b2*(beta.new-beta.old) )
  if(ppar$verbose) {
    cat(sprintf("beta: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                beta.old, beta.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("beta: %g->%g: logr=%g\n",
                beta.old, beta.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["beta"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["beta"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["beta"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}


#' mu (MH) MEW
#'
#' @export
update.mu.mew <- function(state, datlist, fpar, ppar, model,
                          nm=c(alpha="alpha",
                               beta="beta",
                               mu="mu",
                               nu="nu")) {
  last.seed <- .Random.seed
  state.old <- state
  mu.old <- state.old[[nm["mu"]]]
  mu.new <- exp( rnorm(1, log(mu.old), ppar$sd.log.mu) )
  state[[nm["mu"]]] <- mu.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + fpar$s1*log(mu.new/mu.old) 
                   - fpar$s2*(mu.new-mu.old) )
  if(ppar$verbose) {
    cat(sprintf("mu: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                mu.old, mu.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("mu: %g->%g: logr=%g\n",
                mu.old, mu.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["mu"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["mu"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["mu"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}


#' nu (MH) MEW
#'
#' @export
update.nu.mew <- function(state, datlist, fpar, ppar, model,
                          nm=c(alpha="alpha",
                               beta="beta",
                               mu="mu",
                               nu="nu")) {
  last.seed <- .Random.seed
  state.old <- state
  nu.old <- state.old[[nm["nu"]]]
  nu.new <- exp( rnorm(1, log(nu.old), ppar$sd.log.nu) )
  state[[nm["nu"]]] <- nu.new
  state$llike <- llikef(state, datlist, fpar, model)
  state$lprior <- lpriorf(state, fpar, model)
  log.r <- (state$llike - state.old$llike)
  log.r <- ( log.r + fpar$t1*log(nu.new/nu.old) 
                   - fpar$t2*(nu.new-nu.old) )
  if(ppar$verbose) {
    cat(sprintf("nu: (%g;%g) %g->%g: logr=%g\n",
                state.old$llike, state$llike,
                nu.old, nu.new, log.r))
  }
  if(is.nan(log.r) || is.na(log.r) || length(log.r)==0 || is.nan(state$lprior)) { ##!!==
    cat(sprintf("nu: %g->%g: logr=%g\n",
                nu.old, nu.new, log.r))
    if(ppar$interactive) browser()
    update.name <- nm["nu"]
    fname <- paste0("dump-",model,"-",update.name,"-",gsub(" ","-",date()),".Rdata")
    fname <- gsub(":","-",fname,fixed=TRUE)
    save(update.name, state.old, state, datlist, fpar, ppar, model, last.seed,
         file=fname)
  }
  if(runif(1)<exp(log.r)) {
    # accept
    state$accepted[nm["nu"]] <- 1
    if(ppar$verbose) cat("+")
  } else {
    # reject
    state <- state.old
    state$accepted[nm["nu"]] <- 0
    if(ppar$verbose) cat("-")
  }
  
  return(state)
}




