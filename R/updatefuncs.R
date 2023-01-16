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
                        a1=1, a2=0.1,
                        b1=1, b2=0.1,
                        f1=1, f2=0.1)
    }
    if(is.null(update.par)) {
      update.par <- list(sd.log.eta=0.01,
                         sd.log.theta=0.01,
                         sd.logit.v=0.01,
                         sd.log.w=0.01,
                         sd.log.alpha=0.01)
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
    update_parnames <- c(parnames,"wvec")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=TRUE, # = all support points are updated each time
                 ksim=kmax,   # only used if *not* doing a sweep update
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
    epar$vvec <- rbeta(kmax, 1, epar$alpha)
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
    # update the weights (uvec, wvec) and compute lambda0
    state <- augment.state(state, fpar)
    # log L of the current state
    state$llike <- llikef(state, datlist, fpar, model)
    # log prior of the current state
    state$lprior <- lpriorf(state, fpar, model)
    # was the last update accepted?
    state$accepted <- rep(0,length(ppar$update))
    names(state$accepted) <- names(ppar$update)
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
augment.state <- function(state, fpar) {
  # given vvec, compute unscaled weights uvec [with sum(uvec)=1]
  # and scaled weights wvec=gamma*uvec
  # [note that vvec[kmax] exists, but is never used for anything]
  state$uvec <- state$vvec*c(1,cumprod(1-state$vvec[-fpar$kmax]))
  state$uvec[fpar$kmax] <- 1-sum(state$uvec[-fpar$kmax])
  state$wvec <- state$gamma * state$uvec
  # compute lambda0
  state$lambda0 <- state$gamma*state$eta
  return(state)
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
  if(model=="IFR") {
    state <- update_state.ifr(state, datlist, fpar, ppar)    
  } else if(model=="DFR") {
    invisible()
  } else {
    stop("Specified model has not been implemented")
  }
  return(state)
}

update_state.ifr <- function(state, datlist, fpar, ppar) { 
  # Update an IFR state
  model <- "IFR"
  
  # update eta ##!!==  
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    eta.old <- state$eta
    # log normal random walk proposal
    eta.new <- exp( rnorm(1,log(eta.old),ppar$sd.log.eta))
    llike.old <- state$llike
    state$eta <- eta.new
    llike.new <- llikef(state, datlist, fpar, model)
    log.r <- ( llike.new - llike.old 
               +log(eta.new/eta.old) - fpar$nu*(eta.new-eta.old) )
    if(ppar$verbose) {
      cat(sprintf("eta: %g->%g: logr=%g\n",
                  eta.old, eta.new, log.r))
    }
    if(runif(1)<exp(log.r)) {
      # accepted
      state$accepted["eta"] <- 1
      state$llike <- llike.new
      state$lprior <- lpriorf(state, fpar, model)
      state$lambda0 <- state$gamma*state$eta
    } else {
      # reject
      state$eta <- eta.old
      state$accepted["eta"] <- 0
    }
    if(ppar$verbose) cat("\n")
  }
  
  # update gamma !!==: Gibbs update
  if(ppar$update["gamma"]) {
    if(ppar$verbose) cat("gamma:")
    gamma.old <- state$gamma
    cc <- sum(int.lambda.func(datlist$tvec, 
                              model.list=list(model=model, 
                                              kmax=fpar$kmax,
                                              lambda0=state$eta, # NB - unscaled lambda0 used here
                                              thetavec=state$thetavec, 
                                              wvec=state$uvec),  # NB - unscaled weights used here
                              fpar$use.Cpp))
    g1star <- state$alpha+datlist$n0
    g2star <- state$beta+cc
    state$gamma <- rgamma(1, g1star, g2star)
    state$accepted["gamma"] <- 1
    state$wvec <- state$gamma*state$uvec
    state$lambda0 <- state$gamma*state$eta
    state$llike <- llikef(state, datlist, fpar, model)
    state$lprior <- lpriorf(state, fpar, model)
    if(ppar$verbose) {
      cat(sprintf("gamma: %g->%g: logr=Gibbs\n",
                  gamma.old, state$gamma))
    }
    #cat(sprintf("Gamma: %f->%f: r=%f\n",gamma.old, gamma.new, exp(log.r)))
    #cat(sprintf("q(.|a,b) = %f %f\n",
    #            state$alpha+datlist$n0, state$beta+cc))
    #cat(sprintf("  (llike,prior,lq): (%.3f %.3f %.3f)-> (%.3f %.3f %.3f)\n",
    #            llike.old, lprior.old, lq.new,
    #            llike.new, lprior.new, lq.old))
    if(ppar$verbose) cat("\n")
  }
  
  # update thetavec ##!!==
  if(ppar$update["thetavec"]) {
    if(ppar$verbose) cat("thetavec:")
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
    llike.old <- state$llike
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
        state$llike <- llike.new
        state$lprior <- lpriorf(state, fpar, model)
      } else {
        # reject
        state$thetavec[k] <- theta.old
      }
    }
    state$accepted["thetavec"] <- naccepted
    if(ppar$verbose) cat("\n")
  }
  
  # update vvec ##!!==
  if(ppar$update["vvec"]) {
    if(ppar$verbose) cat("vvec:")
    # logistic Normal proposal
    if(ppar$ksweep) {
      # sweep update
      ksamplevec <- 1:(fpar$kmax-1)
    } else {
      # random support point update: only select from 1:(kmax-1)
      ksamplevec <- sample(fpar$kmax-1, ppar$ksim,
                           prob=state$uvec[-fpar$kmax], replace=TRUE)
    }
    naccepted <- 0
    for(k in ksamplevec) {
      llike.old <- state$llike
      ukmax.old <- state$uvec[fpar$kmax]
      v.old <- state$vvec[k]
      v.new <- expit( rnorm(1,logit(v.old),ppar$sd.logit.v) )
      state$vvec[k] <- v.new
      state <- augment.state(state, fpar)
      ukmax.new <- state$uvec[fpar$kmax]
      llike.new <- llikef(state, datlist, fpar, model)
      
      log.r <- (state$alpha)*log((1-v.new)/(1-v.old))
      if(ppar$ksweep) {  # adjust for whether:
        # the component is being updated in sequence
        log.r <- log.r + log(v.new/v.old)
      } else {
        # the component was selected at random
        log.r <- log.r + 2*log(v.new/v.old) + log((1-ukmax.old)/(1-ukmax.new))
      }
      log.r <- log.r + llike.new - llike.old
      if(ppar$verbose) cat(k)
      if(ppar$verbose) cat(sprintf("[%g;%g;%g;%g;%g]\n",
                                   v.old,v.new,llike.old,llike.new,
                                   exp(log.r)))
      if(runif(1)<exp(log.r)) {
        # accepted
        naccepted <- naccepted + 1
        state$llike <- llike.new
        state$lprior <- lpriorf(state, fpar, model)
        if(ppar$verbose) cat("+")
      } else {
        # reject
        state$vvec[k] <- v.old
        state <- augment.state(state, fpar)
        if(ppar$verbose) cat("-")
      }
    }
    state$accepted["vvec"] <- naccepted
    if(ppar$verbose) cat("\n")
  }
  
  # update alpha ##!!==
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    alpha.old <- state$alpha
    a1star <- fpar$a1+fpar$kmax-1
    a2star <- ( fpar$a2-log(state$beta*state$gamma)
                -log(state$uvec[fpar$kmax]) )
    # log normal random walk proposal
    alpha.new <- exp( rnorm(1,log(alpha.old),ppar$sd.log.alpha) )
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
      state$lprior <- lpriorf(state, fpar, model)
      state$accepted["alpha"] <- 1
    } else {
      # reject
      state$alpha <- alpha.old
      state$accepted["alpha"] <- 0
    }
    if(ppar$verbose) cat("\n")
  }
  
  # update beta [OK]
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    # Gibbs update
    b1star <- fpar$b1+state$alpha
    b2star <- fpar$b2+state$gamma
    state$beta <- rgamma(1, b1star, b2star)
    state$lprior <- lpriorf(state, fpar, model)
    state$accepted["beta"] <- 1
    if(ppar$verbose) cat("\n")
  }
  
  # update phi [OK]
  if(ppar$update["phi"]) {
    if(ppar$verbose) cat("phi:")
    # Gibbs update
    f1star <- fpar$f1+fpar$kmax
    f2star <- fpar$f2+sum(state$thetavec)
    state$phi <- rgamma(1, f1star, f2star)
    state$lprior <- lpriorf(state, fpar, model)
    state$accepted["phi"] <- 1
    if(ppar$verbose) cat("\n")
  }
  
  # update wvec ##!!==
  if(ppar$update["wvec"]) {
    if(ppar$verbose) cat("wvec:")
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
    if(ppar$verbose) cat("wvec:")
    for(k in ksamplevec) {
      state.old <- state
      llike.old <- state$llike
      gamma.old <- state$gamma
      vvec.old <- state$vvec
      w.old <- state$wvec[k]
      w.new <- exp( rnorm(1,log(w.old),ppar$sd.log.w) )
      state$wvec[k] <- w.new
      gamma.new <- gamma.old-w.old+w.new
      state$gamma <- gamma.new
      state$uvec <- state$wvec/state$gamma
      state$vvec <- state$wvec/(rev(cumsum(rev(state$wvec))))
      state$lambda0 <- state$gamma*state$eta
      
      llike.new <- llikef(state, datlist, fpar, model)
      log.r <- ( (state$beta)*(gamma.new-gamma.old)
                 +sum( log(state$vvec[1:(fpar$kmax-1)]
                           /vvec.old[1:(fpar$kmax-1)])  
                       +(state$alpha-1)*log(
                         (1-state$vvec[1:(fpar$kmax-1)])
                         /(1-vvec.old[1:(fpar$kmax-1)]))  
                 )
      )
      if(ppar$ksweep) {  # adjust for whether:
        # the component is being updated in sequence
        log.r <- (log.r 
                  + (state$alpha-1)*log(gamma.new/gamma.old)
                  - log(w.new/w.old)
        )
      } else {
        # the component was selected at random
        log.r <- (log.r 
                  + (state$alpha-2)*log(gamma.new/gamma.old)
                  + log(w.new/w.old)
        )
      }
      log.r <- log.r + llike.new - llike.old
      if(ppar$verbose) cat(k)
      if(runif(1)<exp(log.r)) {
        # accepted
        naccepted <- naccepted + 1
        state$llike <- llike.new
        state$lprior <- lpriorf(state, fpar, model)
        if(ppar$verbose) cat("+")
      } else {
        # reject
        state <- state.old
        if(ppar$verbose) cat("-")
      }
    }
    if(ppar$verbose) cat("\n")
    state$accepted["wvec"] <- naccepted
  }
  
  return(state)
}

####################################################################################
