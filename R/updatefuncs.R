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
  
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  
  # update eta ##!!==  
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    state <- update.eta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update gamma ##!!==  
  if(ppar$update["gamma"]) {
    if(ppar$verbose) cat("gamma:")
    state <- update.gamma.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update thetavec ##!!==  
  if(ppar$update["thetavec"]) {
    if(ppar$verbose) cat("thetavec:")
    state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update vvec ##!!==  
  if(ppar$update["vvec"]) {
    if(ppar$verbose) cat("vvec:")
    state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update alpha ##!!==  
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update beta ##!!==  
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update phi ##!!==  
  if(ppar$update["phi"]) {
    if(ppar$verbose) cat("phi:")
    state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update wvec ##!!==  
  if(ppar$update["wvec"]) {
    if(ppar$verbose) cat("wvec:")
    state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # swap elements of thetavec ##!!==  
  if(ppar$update["thetaswap"]) {
    if(ppar$verbose) cat("thetawsap:")
    state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  return(state)
}


update_state.lwb <- function(state, datlist, fpar, ppar, model) { 
  # Update a LWB state
  state$count <- state$count + 1
  parnm <- names(state)
  names(parnm) <- parnm
  
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
  parnm <- names(state)
  names(parnm) <- parnm
  
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
  parnm <- names(state)
  names(parnm) <- parnm
  
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
  parnm <- names(state)
  names(parnm) <- parnm
  
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
