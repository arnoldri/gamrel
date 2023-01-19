# Updating functions

################################################################################
#' Update the current state
#' 
#' @export
update_state <- function(state, datlist, fpar, ppar, model) {
  # update the state
  old.state <- state
  if(model%in%c("IFR","DFR")) {
    state <- update_state.ifrdfr(state, datlist, fpar, ppar, model)    
  } else if(model=="LWB") {
    state <- update_state.lwb(state, datlist, fpar, ppar, model)    
  } else if(model=="SBT") {
    state <- update_state.sbt(state, datlist, fpar, ppar, model)    
  } else if(model=="MBT") {
    state <- update_state.mbt(state, datlist, fpar, ppar, model)    
  } else if(model=="LCV") {
    state <- update_state.lcv(state, datlist, fpar, ppar, model)    
  } else {
    stop("Specified model has not been implemented")
  }
  return(state)
}

update_state.ifrdfr <- function(state, datlist, fpar, ppar, model) { 
  # Update an IFR or DFR state
  state$count <- state$count + 1
  
  # update eta ##!!==  
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    state.old <- state
    eta.old <- state.old$eta
    eta.new <- exp( rnorm(1, log(eta.old), ppar$sd.log.eta) )
    state$eta <- eta.new
    state$lambda0 <- state$gamma*state$eta
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
      state$accepted["eta"] <- 1
      if(ppar$verbose) cat("+")
    } else {
      # reject
      state <- state.old
      state$accepted["eta"] <- 0
      if(ppar$verbose) cat("-")
    }
    if(ppar$verbose) cat("\n")
  }
  
  # update gamma ##!!==  
  if(ppar$update["gamma"]) {
    if(ppar$verbose) cat("gamma:")
    llike.old <- state$llike
    gamma.old <- state$gamma
    temp.model.list <- list(model=model, kmax=kmax, 
                            lambda0=state$eta, ## Note: unscaled
                            thetavec=state$thetavec,
                            wvec=state$uvec)   ## Note: unscaled
    cc <- sum(int.lambda.func(datlist$tvec, model.list=temp.model.list))
    g1star <- state$alpha+datlist$n0
    g2star <- state$beta + cc
    gamma.new <- rgamma(1, g1star, g2star)
    state$gamma <- gamma.new
    state$lambda0 <- state$gamma*state$eta
    state$wvec <- state$gamma*state$uvec
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    state$accepted["gamma"] <- 1
    if(ppar$verbose) {
      cat(sprintf("gamma: (%g;%g) %g->%g: logr=Gibbs\n",
                  llike.old, state$llike,
                  gamma.old, gamma.new))
      cat("+")
    }
    if(ppar$verbose) cat("\n")
  }
  
  # update thetavec ##!!==  
  if(ppar$update["thetavec"]) {
    if(ppar$verbose) cat("thetavec:")
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
      theta.old <- state$thetavec[k]
      theta.new <- exp( rnorm(1,log(theta.old),ppar$sd.log.theta) )
      state$thetavec[k] <- theta.new
      state$llike <- llikef(state, datlist, fpar, model)
      state$lprior <- lpriorf(state, fpar, model)
      log.r <- state$llike - state.old$llike
      log.r <- (log.r - state$phi*(theta.new-theta.old)
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
        if(ppar$verbose) cat("+"); cat(naccepted); cat(";")
      } else {
        # reject
        state <- state.old
        if(ppar$verbose) cat("-")
      }
    }
    state$accepted["thetavec"] <- naccepted/length(ksamplevec)
    if(ppar$verbose) cat("\n")
  }
  
  # update vvec ##!!==  
  if(ppar$update["vvec"]) {
    if(ppar$verbose) cat("vvec:")
    sweep <- as.logical(rbinom(1,1,ppar$psweep))
    if(ppar$verbose && sweep) cat("(sweep):")
    if(sweep || ppar$ksweep) {
      ksamplevec <- 1:(fpar$kmax-1)
    } else {
      ksamplevec <- sample(fpar$kmax-1, ppar$ksim,
                           prob=state$uvec[-fpar$kmax], replace=TRUE)
    }
    naccepted <- 0
    for(k in ksamplevec) {
      state.old <- state
      ukmax.old <- state$uvec[fpar$kmax]
      v.old <- state$vvec[k]
      v.new <- expit( rnorm(1,logit(v.old),ppar$sd.logit.v) )
      state$vvec[k] <- v.new
      cp <- cumprod(1-state$vvec[-fpar$kmax])
      state$uvec <- state$vvec*c(1,cp)
      state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
      #state$uvec[fpar$kmax] <- 1-sum(state$uvec[-fpar$kmax])
      state$wvec <- state$gamma * state$uvec
      ukmax.new <- state$uvec[fpar$kmax]
      state$llike <- llikef(state, datlist, fpar, model)
      state$lprior <- lpriorf(state, fpar, model)
      log.r <- state$llike - state.old$llike
      if(sweep || ppar$ksweep) {
        # updating in sequence
        log.r <- (log.r + state$alpha*log((1-v.new)/(1-v.old)) 
                  + log(v.new/v.old)
        )
      } else {
        # updating at random
        log.r <- (log.r + state$alpha*log((1-v.new)/(1-v.old)) 
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
    state$accepted["vvec"] <- naccepted/length(ksamplevec)
    if(ppar$verbose) cat("\n")
  }
  
  # update alpha ##!!==  
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    state.old <- state
    alpha.old <- state.old$alpha
    alpha.new <- exp( rnorm(1, log(alpha.old), ppar$sd.log.alpha) )
    state$alpha <- alpha.new
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    log.ukmax <- log(state$uvec[fpar$kmax]) 
    #log.ukmax - sum(log(1-state$vvec[-fpar$kmax])) 
    a1star <- fpar$a1 + fpar$kmax - 1
    a2star <- fpar$a2 - log.ukmax - log(state$beta*state$gamma)
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
      state$accepted["alpha"] <- 1
      if(ppar$verbose) cat("+")
    } else {
      # reject
      state <- state.old
      state$accepted["alpha"] <- 0
      if(ppar$verbose) cat("-")
    }
    if(ppar$verbose) cat("\n")
  }
  
  # update beta ##!!==  
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    beta.old <- state$beta
    llike.old <- state$llike
    b1star <- fpar$b1 + state$alpha
    b2star <- fpar$b2 + state$gamma
    beta.new <- rgamma(1, b1star, b2star)
    state$beta <- beta.new
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    state$accepted["beta"] <- 1
    if(ppar$verbose) {
      cat(sprintf("beta: (%g;%g) %g->%g: logr=Gibbs\n",
                  llike.old, state$llike,
                  beta.old, beta.new))
      cat("+")
    }
    if(ppar$verbose) cat("\n")
  }
  
  # update phi ##!!==  
  if(ppar$update["phi"]) {
    if(ppar$verbose) cat("phi:")
    phi.old <- state$phi
    llike.old <- state$llike
    f1star <- fpar$f1 + fpar$kmax
    f2star <- fpar$f2 + sum(state$thetavec)
    phi.new <- rgamma(1, f1star, f2star)
    state$phi <- phi.new
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    state$accepted["phi"] <- 1
    if(ppar$verbose) {
      cat(sprintf("phi: (%g;%g) %g->%g: logr=Gibbs\n",
                  llike.old, state$llike,
                  phi.old, phi.new))
      cat("+")
    }
    if(ppar$verbose) cat("\n")
  }
  
  # update wvec ##!!==  
  if(ppar$update["wvec"]) {
    if(ppar$verbose) cat("wvec:")
    # Do not update wvec[kmax]
    #sweep <- as.logical(rbinom(1,1,ppar$psweep))
    sweep <- FALSE # do not do a global update of every wvec[k]
    if(ppar$verbose && sweep) cat("(sweep):")
    if(sweep || ppar$ksweep) {
      ksamplevec <- 1:(fpar$kmax-1)
    } else {
      ksamplevec <- sample(fpar$kmax-1, ppar$ksim,
                           prob=state$uvec[-fpar$kmax], replace=TRUE)
    }
    naccepted <- 0
    for(k in ksamplevec) {
      state.old <- state
      vvec.old <- state$vvec
      gamma.old <- state$gamma
      wkmax.old <- state$wvec[fpar$kmax]
      w.old <- state$wvec[k]
      w.new <- exp( rnorm(1,log(w.old),ppar$sd.log.w) )
      state$wvec[k] <- w.new
      gamma.new <- sum(state$wvec)
      uvec.new <- state$wvec/gamma.new
      uvec.new <- pmax(.Machine$double.neg.eps, uvec.new) # needed to stabilise
      uvec.new <- uvec.new/sum(uvec.new)
      
      # v values only change for indices up to and including k
      vvec.new <- vvec.old
      vvec.new[1:k] <- wvec[1:k]/(c(1, 1-cumsum(wvec[-kmax]))[1:k])
      state$gamma <- gamma.new
      state$vvec <- vvec.new
      state$uvec <- uvec.new
      state$wvec <- state$gamma*state$uvec
      state$lambda0 <- state$gamma*state$eta
      wkmax <- state$wvec[fpar$kmax]
      
      state$llike <- llikef(state, datlist, fpar, model)
      state$lprior <- lpriorf(state, fpar, model)

      log.r <- state$llike - state.old$llike
      log.r <- (log.r - state$beta*(gamma.new-gamma.old)
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
    state$accepted["wvec"] <- naccepted/length(ksamplevec)
    if(ppar$verbose) cat("\n")
  }
  
  # swap elements of thetavec ##!!==  
  if(ppar$update["thetaswap"]) {
    if(ppar$verbose) cat("thetawsap:")
    ksamplevec1 <- sample(fpar$kmax, ppar$kswap,
                          prob=state$uvec, replace=TRUE)
    ksamplevec2 <- sample(fpar$kmax, ppar$kswap, replace=TRUE)
    naccepted <- 0
    for(i in 1:ppar$kswap) {
      k1 <- ksamplevec1[i]
      k2 <- ksamplevec2[i]
      state.old <- state
      state$thetavec[c(k1,k2)] <- state.old$thetavec[c(k2,k1)]
      state$llike <- llikef(state, datlist, fpar, model)
      state$lprior <- lpriorf(state, fpar, model)
      log.r <- state$llike - state.old$llike
      log.r <- (log.r + log(state$wvec[k2]/state$wvec[k1]))
      if(ppar$verbose) {
        cat(sprintf("thetaswap[%d,%d]: (%g;%g) %g->%g: logr=%g\n",
                    k1, k2, state.old$llike, state$llike,
                    state.old$thetavec[k1], state.old$thetavec[k2], log.r))
      }
      if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
        cat(sprintf("thetaswap[%d,%d]: %g->%g: logr=%g\n",
                    k1, k2, state.old$thetavec[k1], state.old$thetavec[k2], log.r))
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
    state$accepted["thetaswap"] <- naccepted/length(ksamplevec)
    if(ppar$verbose) cat("\n")
  }
  
  return(state)
}


update_state.lwb <- function(state, datlist, fpar, ppar, model) { 
  # Update a LWB state
  state$count <- state$count + 1
  stop("Not yet implemented")
  
  # update eta !!==
  # update a !!==
  # update gamma !!==
  # update thetavec !!==
  # update vvec !!==
  # update alpha !!==
  # update beta !!==
  # update phi !!==
  # update wvec !!==
  # thetaswap !!==
  
  return(state)
}
update_state.sbt <- function(state, datlist, fpar, ppar, model) { 
  # Update a SBT state
  state$count <- state$count + 1
  stop("Not yet implemented")

  # update eta !!==
  
  # update gamma1 !!==
  # update thetavec1 !!==
  # update vvec1 !!==
  # update alpha1 !!==
  # update beta1 !!==
  # update phi1 !!==
  # update wvec1 !!==
  # thetaswap1 !!==
  
  # update gamma2 !!==
  # update thetavec2 !!==
  # update vvec2 !!==
  # update alpha2 !!==
  # update beta2 !!==
  # update phi2 !!==
  # update wvec2 !!==
  # thetaswap2 !!==
  
  return(state)
}
update_state.mbt <- function(state, datlist, fpar, ppar, model) { 
  # Update a MBT state
  state$count <- state$count + 1
  stop("Not yet implemented")
  
  # update pival !!==
  
  # update eta1 !!==
  # update gamma1 !!==
  # update thetavec1 !!==
  # update vvec1 !!==
  # update alpha1 !!==
  # update beta1 !!==
  # update phi1 !!==
  # update wvec1 !!==
  # thetaswap1 !!==
  
  # update eta2 !!==
  # update gamma2 !!==
  # update thetavec2 !!==
  # update vvec2 !!==
  # update alpha2 !!==
  # update beta2 !!==
  # update phi2 !!==
  # update wvec2 !!==
  # thetaswap2 !!==
  
  return(state)
}
update_state.lcv <- function(state, datlist, fpar, ppar, model) { 
  # Update a LCV state
  state$count <- state$count + 1
  stop("Not yet implemented")
  
  # update lambda0 !!==
  # update w0 !!==
  # update gamma !!==
  # update thetavec !!==
  # update vvec !!==
  # update alpha !!==
  # update beta !!==
  # update phi !!==
  # update wvec !!==
  # thetaswap !!==
  
  return(state)
}

####################################################################################
