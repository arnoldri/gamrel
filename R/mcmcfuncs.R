#############################################################################
# MCMC functions
#############################################################################
makew <- function(vvec) {
  # make the weight vector
  kmax <- length(vvec)
  wvec <- c(vvec[-kmax],1)*c(1,cumprod(1-vvec[-kmax]))
  return(wvec)
}
makev <- function(wvec) {
  # make the v vector
  kmax <- length(wvec)
  vvec <- wvec/c(1, 1-cumsum(wvec[-kmax]))
  return(vvec)
}
###############################################################################
#' Burn an MCMC chain
#' 
#' @export
burn.chain <- function(nburn, state, datlist, fpar, ppar, model) {
  if(nburn>0) {
    for(i in 1:nburn) state <- update.state(state, datlist, fpar, ppar, model)
  }
  return(state)
}
#' Run an MCMC chain
#' 
#' @export
run.chain <- function(nstore, nthin=1, nburn=0,
                      state, datlist, fpar, ppar, model,
                      do.plot=FALSE, show.progress=NULL) {
  tstart <- Sys.time()
  if(nburn>0) state <- burn.chain(nburn, state, datlist, fpar, ppar, model)
  svec <- state.as.vector(state, fpar, model)
  smat <- array(NA, dim=c(nstore, length(svec)))
  smat[1,] <- svec
  dimnames(smat)[[2]] <- fpar$statevnames
  dimnames(smat)[[2]][grep("accepted",dimnames(smat)[[2]])] <- (
    paste0("accepted.",names(ppar$update)))

  if(do.plot) plot.state(state, datlist, fpar, ppar, model, add=FALSE)
  for(i in 2:nstore) {
    if(!is.null(show.progress)) {
      if(i%%show.progress==0) {
        cat(".")
        if(i%%(20*show.progress)==0) cat("\n")
      }
    }
    state <- burn.chain(nthin, state, datlist, fpar, ppar, model)
    if(do.plot) plot.state(state, datlist, fpar, ppar, model, add=FALSE)
    smat[i,] <- state.as.vector(state, fpar, model)
  }
  if(!is.null(show.progress)) cat("\n")
  telapsed <- Sys.time()-tstart
  print(telapsed)
  return(smat)
}
###############################################################################
#' Initialise objects for an MCMC chain
#' 
#' @export
init.objects <- function(tvec, n1, tau,
                         kmax=100,
                         priorpar=list(a1=1, a2=0.1,
                                       b1=1, b2=0.1,
                                       d1=1, d2=0.1,
                                       f1=1, f2=0.1,
                                       rexp.rate=1),
                         update.par=list(sd.log.alpha=0.01,
                                         sd.logit.v=0.01,
                                         sd.log.theta=0.01),
                         model="IFR",
                         use.Cpp=FALSE,
                         seed=NULL,
                         generate="random") {  # generate can be "fixed" or "random"
  if(!is.null(seed)) set.seed(seed)
  if(model%in%c("IFR","DFR")) {
    # fixed parameters
    fpar <- list(kmax=kmax,                      # sum truncation point
                 a1=priorpar$a1, a2=priorpar$a2, # prior for alpha
                 b1=priorpar$b1, b2=priorpar$b2, # prior for beta
                 d1=priorpar$d1, d2=priorpar$d2, # prior for lambda0
                 f1=priorpar$f1, f2=priorpar$f2, # prior for rexp.rate=phi
                 use.Cpp=use.Cpp) 
    # proposal parameters for updates
    parnames <- c("alpha","beta","phi","lambda0","gamma","vvec","thetavec")
    update <- rep(TRUE, length(parnames))
    names(update) <- parnames
    ppar <- list(ksweep=TRUE,# = all support points are updated each time
                 ksim=kmax,  # only used if *not* doing a sweep update
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.logit.v=update.par$sd.logit.v,
               sd.log.theta=update.par$sd.log.theta,
               update=update,
               verbose=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(alpha=priorpar$a1/priorpar$a2,
                   beta=priorpar$b1/priorpar$b2,
                   lambda0=priorpar$d1/priorpar$d2,
                   phi=priorpar$f1/priorpar$b2)
      epar$gamma <- epar$alpha/epar$beta
    } else {
      epar <- list(alpha=rgamma(priorpar$a1,priorpar$a2),
                   beta=rgamma(priorpar$b1,priorpar$b2),
                   lambda0=rgamma(priorpar$d1,priorpar$d2),
                   phi=rgamma(priorpar$f1,priorpar$b2))
      epar$gamma <- rgamma(epar$alpha,epar$beta)
    }
    epar$vvec <- rbeta(kmax, 1, epar$alpha)
    epar$thetavec <- rexp(kmax, epar$phi)
  } else {
    stop("Specified model has not been implemented")
  }
  # data
  datlist <- list(tvec=tvec,           # uncensored observations
                  n1=n1,               # number of censored observations
                  tau=tau,             # censoring time
                  n=length(tvec)+n1,   # total number of observations
                  n0=length(tvec))     # number of uncensored observations
  # return all these objects
  return(list(epar=epar, datlist=datlist,
              fpar=fpar, ppar=ppar, model=model))
}
#' Log likelihood of the current state
#' 
#' @export
llikef <- function(state, datlist, fpar, model) {
  # log likelihood of observations tvec (censored at tau)
  lambda.vec <- lambda.func(tvec=c(datlist$tvec,datlist$tau),
                            model.list=c(list(model=model),state),
                            use.Cpp=fpar$use.Cpp)
  int.lambda.vec <- int.lambda.func(tvec=c(datlist$tvec,datlist$tau),
                                    model.list=c(list(model=model),state),
                                    use.Cpp=fpar$use.Cpp)

  retval <- ( sum(log(lambda.vec[1:datlist$n0]))
              - sum(int.lambda.vec[1:datlist$n0])
              - datlist$n1*int.lambda.vec[datlist$n0+1]
              )
  return(retval)
}
#' Log prior of the current state
#' 
#' @export
lpriorf <- function(state, fpar, model) {
  # log prior of the state
  if(model%in%c("IFR","DFR")) {
    retval <- 0
    # prior for vvec
    retval <- retval + sum(log(state$alpha) + (state$alpha-1)*log(1-state$vvec))
    # prior for thetavec
    retval <- retval + sum(log(state$phi)-state$phi*state$thetavec)
    # prior for gamma
    retval <- retval + dgamma(state$gamma, state$alpha, state$beta, log=TRUE)
    # prior for alpha
    retval <- retval + dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
    # prior for beta
    retval <- retval + dgamma(state$beta, fpar$b1, fpar$b2, log=TRUE)
    # prior for lambda0
    retval <- retval + dgamma(state$lambda0, fpar$d1, fpar$d2, log=TRUE)
    # prior for phi
    retval <- retval + dgamma(state$phi, fpar$f1, fpar$f2, log=TRUE)
  } else {
    stop("Specified model has not been implemented")
  }
  return(retval)
}
#' Components of the log prior of the current state - return as a vector
#' 
#' @export
lpriorf.vector <- function(state, fpar, model) {
  # log prior of the state
  if(model%in%c("IFR","DFR")) {
    # prior for vvec
    lprior.vvec <-  sum( (log(state$alpha) + (state$alpha-1)*log(1-state$vvec)) )
    # prior for thetavec
    lprior.thetavec <- sum( (log(state$phi)-state$phi*state$thetavec) )
    # prior for gamma
    lprior.gamma <- dgamma(state$gamma, state$alpha, state$beta, log=TRUE)
    # prior for alpha
    lprior.alpha <- dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
    # prior for beta
    lprior.beta <- dgamma(state$beta, fpar$b1, fpar$b2, log=TRUE)
    # prior for lambda0
    lprior.lambda0 <- dgamma(state$lambda0, fpar$d1, fpar$d2, log=TRUE)
    # prior for phi
    lprior.phi <- dgamma(state$phi, fpar$f1, fpar$f2, log=TRUE)
    retval <- c(vvec=lprior.vvec, thetavec=lprior.thetavec,
                gamma=lprior.gamma, alpha=lprior.alpha, beta=lprior.beta,
              lambda0=lprior.lambda0, phi=lprior.phi)
    retval <- c(retval, lprior=sum(retval))
  } else {
    stop("Specified model has not been implemented")
  }
  return(retval)
}
update.vw <- function(state, fpar) {
  # given vvec, compute unscaled weights uvec [with sum(uvec)=1]
  # and scaled weights wvec=gamma*uvec
  # [note that vvec[kmax] exists, but is never used for anything]
  state$uvec <- state$vvec*c(1,cumprod(1-state$vvec[-fpar$kmax]))
  state$uvec[fpar$kmax] <- 1-sum(state$uvec[-fpar$kmax])
  state$wvec <- state$gamma * state$uvec
  return(state)
}
#' Construct the state object from the parameter vectors and data
#' 
#' @export
make.state <- function(epar, datlist, fpar, ppar, model) {
  if(model%in%c("IFR","DFR")) {
    state <- epar
    # complete state with useful quantities
    # update the weights (uvec, wvec)
    state <- update.vw(state, fpar)
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
################################################################################
#' Initialise the fixed parameters (fpar)
#' 
#' @export
init.fpar <- function(state, fpar, model) {
  # add some extra stuff to fpar
  fpar$statenames <- names(state)
  fpar$stackind <- stack(state)$ind
  fpar$statevnames <- names(unlist(sapply(sapply(state,length),
                                          function(i) 1:i)))
  return(fpar)
}
#' Convert the state from a list to a vector representation
#' 
#' @export
state.as.vector <- function(state, fpar, model) {
  statevec <- stack(state)$values
  return(statevec)
}
#' Convert the state from a vector to a list representation
#' 
#' @export
vector.as.state <- function(statevec, datlist, fpar, ppar, model) {
  if(model%in%c("IFR","DFR")) {
    state <- unstack(data.frame(values=statevec,ind=fpar$stackind))
    names(state$accepted) <- names(ppar$update)
  } else {
    stop(paste0("Model ",model," not recognised"))
  }
  return(state)
}
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
                       type=1, add=FALSE) {
  # plot the state
  if(model=="IFR") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w),
           main="")
    }
    points(state$thetavec, state$wvec)
  } else if(model=="DFR") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w),
           main="")
    }
    points(state$thetavec, state$wvec)
  } else {
    stop(paste0("Model ",model," not recognised"))
  }
  invisible()
}
#######################################################################
