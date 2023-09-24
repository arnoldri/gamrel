#############################################################################
# Generic MCMC functions
#
# state = list representation of the current state
# datlist = data
# fpar = fixed parameters associated with the full model (including all priors)
# ppar = updating parameters associated with the samplers
# model = name of the model (character string)
#
###############################################################################
#' Burn an MCMC chain
#' 
#' @export
burn.chain <- function(nburn, state, datlist, fpar, ppar, model) {
  if(nburn>0) {
    for(i in 1:nburn) state <- update_state(state, datlist, fpar, ppar, model)
  }
  return(state)
}
#' Run multiple MCMC chains - with starting states specified in state.list
#' 
#' @export
run.multiple.chains <- function(state.list, nstore, nthin=1, nburn=0,
                                datlist, fpar, ppar, model,
                                do.plot=FALSE, show.progress=NULL) {
  chain.list <- list()
  for(i in 1:length(state.list)) {
    chain.list[[i]] <- run.chain(nstore=nstore, nthin=nthin, nburn=nburn,
                                 state=state.list[[i]], datlist=datlist, fpar=fpar,
                                 ppar=ppar, model=model, do.plot=do.plot,
                                 show.progress=show.progress)
  }
  return(chain.list)
}

#' Run an MCMC chain
#' 
#' @export
run.chain <- function(nstore, nthin=1, nburn=0,
                      state, datlist, fpar, ppar, model,
                      do.plot=FALSE, show.progress=NULL) {
  tstart <- Sys.time()
  nthin <- max(1,nthin)
  if(nburn>0) state <- burn.chain(nburn, state, datlist, fpar, ppar, model)
  svec <- state.as.vector(state, fpar, model)
  smat <- array(NA, dim=c(nstore, length(svec)))
  smat[1,] <- svec
  dimnames(smat)[[2]] <- fpar$statevnames
  dimnames(smat)[[2]][grep("accepted",dimnames(smat)[[2]])] <- (
    paste0("accepted.",names(ppar$update)))

  if(do.plot) plot_state(state, datlist, fpar, ppar, model, add=FALSE)
  for(i in 2:nstore) {
    if(!is.null(show.progress)) {
      if(i%%show.progress==0) {
        cat(".")
        if(i%%(20*show.progress)==0) cat("\n")
      }
    }
    state <- burn.chain(nthin, state, datlist, fpar, ppar, model)
    if(do.plot) plot_state(state, datlist, fpar, ppar, model, add=FALSE)
    smat[i,] <- state.as.vector(state, fpar, model)
  }
  if(!is.null(show.progress)) cat("\n")
  telapsed <- Sys.time()-tstart
  print(telapsed)
  return(smat)
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
  state <- unstack(data.frame(values=statevec,ind=fpar$stackind))
  names(state$accepted) <- names(ppar$update)
  return(state)
}

#' Evaluate the WAIC value from a chain
#' 
#' @param smat Matrix output of an MCMC chain
#' @param datlist Data object
#' @param fpar Fixed parameters of the model
#' @param ppar Parameters used in running the sampler
#' @param model Model name
#' 
#' @export
waicfunc <- function(smat, datlist, fpar, ppar, model) {
  llmat <- t(sapply(1:nrow(smat), 
                    function(j) {
                      state <- vector.as.state(smat[j,], datlist, fpar, ppar, model)
                      sapply(1:datlist$n, function(i) {
                        subdatlist <- list(n=1, n0=as.numeric(datlist$obs[i]),
                                           tvec=datlist$tvec[i], obs=datlist$obs[i])
                        return(llikef(state, subdatlist, fpar, model))
                      })
                    }))
  # Check
  #range(apply(llmat,1,sum)-smat[,"llike"])
  
  waic <- -2*sum(log(apply(exp(llmat),2,mean))) +2*sum(apply(llmat,2,var))
  return(waic)
}
