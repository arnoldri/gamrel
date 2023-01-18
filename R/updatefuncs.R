# Updating functions

###############################################################################
#' Initialise objects for an MCMC chain
#' 
#' @export
init.objects <- function(tvec, obs,
                         kmax=100,
                         prior.par=NULL,
                         update.par=NULL,
                         model="IFR",
                         use.Cpp=FALSE,
                         seed=NULL,
                         generate="random") {  # generate can be "fixed" or "random"
  if(!is.null(seed)) set.seed(seed)
  if(model%in%c("IFR","DFR")) {
    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        a1=1, a2=1,
                        b1=1, b2=1,
                        f1=1, f2=1)
    }
    if(is.null(update.par)) {
      update.par <- list(psweep=0.1, # probability of the sweep move
                         sd.log.eta=0.1,
                         sd.log.theta=0.1,
                         sd.logit.v=0.1,
                         sd.log.w=0.1,
                         sd.log.alpha=0.1)
    }
    # fixed parameters
    parnames <- c("eta","gamma","thetavec","vvec","alpha","beta","phi")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 nu=prior.par$nu,                  # prior for eta (lambda0/gamma)
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 use.Cpp=use.Cpp) 
    # parameters to update
    update_parnames <- c(parnames,"wvec","thetaswap")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 psweep=update.par$psweep,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.w=update.par$sd.log.w,
                 sd.log.alpha=update.par$sd.log.alpha,
                 update=update,
                 verbose=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(eta=1/prior.par$nu,
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=prior.par$a1/prior.par$a2,
                   beta=prior.par$b1/prior.par$b2,
                   phi=prior.par$f1/prior.par$b2)
      epar$gamma <- epar$alpha/epar$beta
    } else {
      epar <- list(eta=rexp(1,prior.par$nu),
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   beta=rgamma(1,prior.par$b1,prior.par$b2),
                   phi=rgamma(1,prior.par$f1,prior.par$b2))
      epar$gamma <- rgamma(1,epar$alpha,epar$beta)
    }
    epar$thetavec <- rexp(kmax, epar$phi)
    epar$vvec <- rbeta.t(kmax, 1, epar$alpha)
  } else {
    stop("Specified model has not been implemented")
  }
  # data
  datlist <- list(n=length(tvec),      # size of data
                  n0=sum(obs),         # number of uncensored observations
                  tvec=tvec,           # failure/censoring times
                  obs=obs)             # vector of indicators: TRUE=observed, FALSE=censored
  # make the state
  state <- make.state(epar, datlist, fpar, ppar, model)
  # augment fpar
  fpar <- augment.fpar(state, fpar, model)
  # return all these objects
  return(list(epar=epar, datlist=datlist,
              fpar=fpar, ppar=ppar, model=model,
              state=state))
}

#' Construct the state object from the parameter vectors and data
#' 
#' @export
make.state <- function(epar, datlist, fpar, ppar, model) {
  if(model%in%c("IFR","DFR")) {
    state <- epar
    # complete state with useful quantities
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec[-fpar$kmax])
    state$uvec <- state$vvec*c(1,cp)
    state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec <- state$gamma * state$uvec
    # compute lambda0
    state$lambda0 <- state$gamma*state$eta
    # log L of the current state
    state$llike <- llikef(state, datlist, fpar, model)
    # log prior of the current state
    state$lprior <- lpriorf(state, fpar, model)
    # was the last update accepted?
    state$accepted <- rep(0,length(ppar$update))
    names(state$accepted) <- names(ppar$update)
    # counter
    state$count <- 0
  } else {
    stop("Specified model has not been implemented")
  }
  return(state)
}

#' Log likelihood of the current state
#' 
#' @export
llikef <- function(state, datlist, fpar, model) {
  # log likelihood of observations tvec
  lambda.vec <- lambda.func(tvec=datlist$tvec[datlist$obs],
                            model.list=c(list(model=model,kmax=fpar$kmax),
                                         state),
                            use.Cpp=fpar$use.Cpp)
  int.lambda.vec <- int.lambda.func(tvec=datlist$tvec,
                                    model.list=c(list(model=model,kmax=fpar$kmax),state),
                                    use.Cpp=fpar$use.Cpp)
  
  retval <- ( sum(log(lambda.vec)) - sum(int.lambda.vec) )
  return(retval)
}
#' Log prior of the current state
#' 
#' @export
lpriorf <- function(state, fpar, model) {
  # log prior of the state
  retval <- sum(lpriorf.vector(state, fpar, model))
  return(retval)
}
#' Components of the log prior of the current state - return as a vector
#' 
#' @export
lpriorf.vector <- function(state, fpar, model) {
  # log prior of the state
  lprior.vec <- rep(NA,length=length(fpar$parnames))
  names(lprior.vec) <- fpar$parnames
  if(model%in%c("IFR","DFR")) {
    # prior for eta
    lprior.vec["eta"] <- dexp(state$eta, fpar$nu)
    # prior for gamma
    lprior.vec["gamma"] <- dgamma(state$gamma, state$alpha, state$beta, log=TRUE)
    # prior for thetavec
    lprior.vec["thetavec"] <- sum( (log(state$phi)-state$phi*state$thetavec) )
    # prior for vvec
    lprior.vec["vvec"] <-  sum( (log(state$alpha) + (state$alpha-1)*log(1-state$vvec[-fpar$kmax])) )
    # prior for alpha
    lprior.vec["alpha"] <- dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
    # prior for beta
    lprior.vec["beta"] <- dgamma(state$beta, fpar$b1, fpar$b2, log=TRUE)
    # prior for phi
    lprior.vec["phi"] <- dgamma(state$phi, fpar$f1, fpar$f2, log=TRUE)
  } else {
    stop("Specified model has not been implemented")
  }
  return(lprior.vec)
}
#' Augment the fixed parameter object (fpar)
#' 
#' @export
augment.fpar <- function(state, fpar, model) {
  # add some extra stuff to fpar
  fpar$statenames <- names(state)
  fpar$stackind <- stack(state)$ind
  fpar$statevnames <- names(unlist(sapply(sapply(state,length),
                                          function(i) 1:i)))
  return(fpar)
}

################################################################################
#' Plot a representation of the current state
#' 
#' @description Plot the current state
#' 
#' @param state The current state
#' @param datlist The data list
#' @param fpar The fixed parameters
#' @param ppar Proposal parameters
#' @param model Type of model (IFR or DFR)
#' @param type Type of plot (ignored at the moment)
#' @param add Logical: add the state to an existing plot?
#' 
#' @export
plot_state <- function(state, datlist, fpar, ppar, model,
                       type=1, main=NULL, add=FALSE, ...) {
  # plot the state
  if(model=="IFR") {
    if(!add) {
      # start a new plot
      if(is.null(main)) {
        main <- sprintf("LP=%.3f; LL=%.3f; LPost=%.3f",
                        state$lprior,state$llike,state$llike+state$lprior)        
      }
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    points(state$thetavec, state$wvec, ...)
  } else if(model=="DFR") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    points(state$thetavec, state$wvec, ...)
  } else {
    stop(paste0("Model ",model," not recognised"))
  }
  invisible()
}
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
update_state.lwb <- function(state, datlist, fpar, ppar, model) { 
  # Update a LWB state
  state$count <- state$count + 1
  stop("Not yet implemented")
  return(state)
}
update_state.sbt <- function(state, datlist, fpar, ppar, model) { 
  # Update a SBT state
  state$count <- state$count + 1
  stop("Not yet implemented")
  return(state)
}
update_state.mbt <- function(state, datlist, fpar, ppar, model) { 
  # Update a MBT state
  state$count <- state$count + 1
  stop("Not yet implemented")
  return(state)
}
update_state.lcv <- function(state, datlist, fpar, ppar, model) { 
  # Update a LCV state
  state$count <- state$count + 1
  stop("Not yet implemented")
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
        if(ppar$verbose) cat("+")
      } else {
        # reject
        state <- state.old
        if(ppar$verbose) cat("-")
      }
      state$accepted["thetavec"] <- naccepted/length(ksamplevec)
    }
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
      state$accepted["vvec"] <- naccepted/length(ksamplevec)
    }
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
    sweep <- as.logical(rbinom(1,1,ppar$psweep))
    if(ppar$verbose && sweep) cat("(sweep):")
    if(sweep || ppar$ksweep) {
      ksamplevec <- 1:(fpar$kmax-1)
    } else {
      ksamplevec <- sample(fpar$kmax-1, ppar$ksim,
                           prob=state$uvec[-fpar$kmax], replace=TRUE)
    }
    naccepted <- 0
    #errcount <- 0 ##!!==
    for(k in ksamplevec) {
      #errcount<-errcount+1; cat(sprintf("%d:%d;\n",k,errcount)) ##!!==
      #options(warn=2) ##!!==
      
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
      vvec.new <- makev(uvec.new)
      vvec.new[vvec.new>=1 | vvec.new<=0] <- 0.5 ## numerical stabilisation
      state$gamma <- gamma.new
      state$vvec <- vvec.new
      state$uvec <- uvec.new
      state$wvec <- state$gamma*state$uvec
      state$lambda0 <- state$gamma*state$eta
      wkmax <- state$wvec[fpar$kmax]
      
      #errcount<-errcount+1; cat(sprintf("llike:%d:%d;\n",k,errcount)) ##!!==
      
      state$llike <- llikef(state, datlist, fpar, model)
      
      #errcount<-errcount+1; cat(sprintf("lprior:%d:%d;\n",k,errcount)) ##!!==
      
      state$lprior <- lpriorf(state, fpar, model)
      #state$lprior <- state$lprior ##!!==
      
      #errcount<-errcount+1; cat(sprintf("log.r:%d:%d;\n",k,errcount)) ##!!==
      
      #browser() ##!!==
      
      log.r <- state$llike - state.old$llike
      log.r <- (log.r - state$beta*(gamma.new-gamma.old)
                + sum(log(vvec.new[-fpar$kmax]/vvec.old[-fpar$kmax])))
      
      #errcount<-errcount+1; cat(sprintf("%d:%d;\n",k,errcount)) ##!!==
      
      if(sweep || ppar$ksweep) {
        # updating sequentially
        log.r <- log.r
      } else {
        # updating at random
        log.r <- (log.r + log(w.new/w.old)
                  + log((gamma.old-wkmax)/(gamma.new-wkmax))
        )
      }
      
      #errcount<-errcount+1; cat(sprintf("%d:%d;\n",k,errcount)) ##!!==
      
      if(ppar$verbose) {
        cat(sprintf("wvec[%d]: (%g;%g) %g->%g: logr=%g\n",
                    k, state.old$llike, state$llike,
                    w.old, w.new, log.r))
      }
      
      #errcount<-errcount+1; cat(sprintf("%d:%d;\n",k,errcount)) ##!!==
      
      if(is.nan(log.r) || is.na(log.r) || length(log.r)==0) { ##!!==
        cat(sprintf("wvec[%d]: %g->%g: logr=%g\n",
                    k, w.old, w.new, log.r))
        browser()
      }
      
      #errcount<-errcount+1; cat(sprintf("%d:%d;\n",k,errcount)) ##!!==
      
      if(runif(1)<exp(log.r)) {
        # accept
        naccepted <- naccepted + 1
        if(ppar$verbose) cat("+")
      } else {
        # reject
        state <- state.old
        if(ppar$verbose) cat("-")
      }
      state$accepted["thetavec"] <- naccepted/length(ksamplevec)
      
      #errcount<-errcount+1; cat(sprintf("%d:%d;\n",k,errcount)) ##!!==
      
    }
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
      state$accepted["thetaswap"] <- naccepted/length(ksamplevec)
    }
    if(ppar$verbose) cat("\n")
  }
  
  return(state)
}



####################################################################################
