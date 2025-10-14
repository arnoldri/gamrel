# Updating functions

################################################################################
#' Update the current state
#' 
#' @export
update_state <- function(state, datlist, fpar, ppar, model) {
  # update the state
  old.state <- state
  if(model%in%c("CON")) {
    state <- update_state.con(state, datlist, fpar, ppar, model)    
  } else if(model%in%c("IFR","DFR")) {
    state <- update_state.ifrdfr(state, datlist, fpar, ppar, model)    
  } else if(model=="LWB") {
    state <- update_state.lwb(state, datlist, fpar, ppar, model)    
  } else if(model=="SBT") {
    state <- update_state.sbt(state, datlist, fpar, ppar, model)    
  } else if(model=="MBT") {
    state <- update_state.mbt(state, datlist, fpar, ppar, model)    
  } else if(model=="LCV") {
    state <- update_state.lcv(state, datlist, fpar, ppar, model)    
  } else if(model%in%c("CIR","CDR")) {
    state <- update_state.circdr(state, datlist, fpar, ppar, model)    
  } else if(model%in%c("CVX")) {
    state <- update_state.cvx(state, datlist, fpar, ppar, model)    
  } else if(model=="MEW") {
    state <- update_state.mew(state, datlist, fpar, ppar, model)    
  } else {
    stop("Specified model has not been implemented")
  }
  return(state)
}

update_state.con <- function(state, datlist, fpar, ppar, model) { 
  # Update a CON state
  state$count <- state$count + 1
  
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm

  # update lambda0 ##!!== OK 
  if(ppar$update["lambda0"]) {
    if(ppar$verbose) cat("lambda0:")
    state <- update.lambda0.v0(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  return(state)
}

update_state.ifrdfr <- function(state, datlist, fpar, ppar, model) { 
  # Update an IFR or DFR state
  state$count <- state$count + 1
  
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update eta ##!!== OK 
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    state <- update.eta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update gamma ##!!== OK 
  if(ppar$update["gamma"]) {
    if(ppar$verbose) cat("gamma:")
    state <- update.gamma.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update thetavec ##!!== OK  
  if(ppar$update["thetavec"]) {
    if(ppar$verbose) cat("thetavec:")
    state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update vvec ##!!== OK
  if(ppar$update["vvec"]) {
    if(ppar$verbose) cat("vvec:")
    state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update alpha ##!!== OK  
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update beta ##!!== OK
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update phi ##!!==  OK
  if(ppar$update["phi"]) {
    if(ppar$verbose) cat("phi:")
    state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update wvec ##!!== OK
  if(ppar$update["wvec"]) {
    if(ppar$verbose) cat("wvec:")
    state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # swap elements of thetavec ##!!== OK  
  if(ppar$update["thetaswap"]) {
    if(ppar$verbose) cat("thetaswap:")
    state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  return(state)
}


update_state.lwb <- function(state, datlist, fpar, ppar, model) { 
  # Update a LWB state
  state$count <- state$count + 1

  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update eta ##!!== OK
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    state <- update.eta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update a ##!!== **
  if(ppar$update["a"]) {
    if(ppar$verbose) cat("a:")
    state <- update.a.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }

  # update gamma ##!!== OK
  if(ppar$update["gamma"]) {
    if(ppar$verbose) cat("gamma:")
    state <- update.gamma.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update thetavec ##!!== OK
  if(ppar$update["thetavec"]) {
    if(ppar$verbose) cat("thetavec:")
    state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update vvec ##!!== OK
  if(ppar$update["vvec"]) {
    if(ppar$verbose) cat("vvec:")
    state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update alpha ##!!== OK
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update beta ##!!== OK
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update phi ##!!== OK
  if(ppar$update["phi"]) {
    if(ppar$verbose) cat("phi:")
    state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update wvec ##!!== OK
  if(ppar$update["wvec"]) {
    if(ppar$verbose) cat("wvec:")
    state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # swap elements of thetavec ##!!== OK  
  if(ppar$update["thetaswap"]) {
    if(ppar$verbose) cat("thetaswap:")
    state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }

  return(state)
}
update_state.sbt <- function(state, datlist, fpar, ppar, model) { 
  # Update a SBT state
  state$count <- state$count + 1
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update eta ##!!== OK 
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    tparnm <- parnm
    tparnm["gamma"] <- "gamma2"
    state <- update.eta.v1(state, datlist, fpar, ppar, model, nm=tparnm)
    if(ppar$verbose) cat("\n")
  }

  for(j in c(1,2)) {  # j=1:DFR component; j=2:IFR component

    tparnm <- c(gamma="gamma",
                thetavec="thetavec",
                vvec="vvec",
                alpha="alpha",
                beta="beta",
                phi="phi",
                uvec="uvec",
                wvec="wvec",
                thetaswap="thetaswap",
                f1="f1", f2="f2")
    ntp <- names(tparnm)
    tparnm <- paste0(tparnm,j)
    names(tparnm) <- ntp
    tparnm <- c(eta="eta",lambda0=ifelse(j==1,NA,"lambda0"),tparnm)
    
    # update gamma ##!!== OK 
    gamma.name <- paste0("gamma",j)
    if(ppar$update[gamma.name]) {
      if(ppar$verbose) cat(gamma.name)
      state <- update.gamma.v2(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update thetavec ##!!== OK  
    thetavec.name <- paste0("thetavec",j)
    if(ppar$update[thetavec.name]) {
      if(ppar$verbose) cat(thetavec.name)
      state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update vvec ##!!== OK
    vvec.name <- paste0("vvec",j)
    if(ppar$update[vvec.name]) {
      if(ppar$verbose) cat(vvec.name)
      state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update alpha ##!!== OK  
    alpha.name <- paste0("alpha",j)
    if(ppar$update[alpha.name]) {
      if(ppar$verbose) cat(alpha.name)
      state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update beta ##!!== OK
    beta.name <- paste0("beta",j)
    if(ppar$update[beta.name]) {
      if(ppar$verbose) cat(beta.name)
      state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update phi ##!!==  OK
    phi.name <- paste0("phi",j)
    if(ppar$update[phi.name]) {
      if(ppar$verbose) cat(phi.name)
      state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update wvec ##!!== OK
    wvec.name <- paste0("wvec",j)
    if(ppar$update[wvec.name]) {
      if(ppar$verbose) cat(wvec.name)
      state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # swap elements of thetavec ##!!== OK
    thetaswap.name <- paste0("thetaswap",j)
    if(ppar$update[thetaswap.name]) {
      if(ppar$verbose) cat(thetaswap.name)
      state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
  }
  
  return(state)
}

update_state.mbt <- function(state, datlist, fpar, ppar, model) { 
  # Update a MBT state
  state$count <- state$count + 1
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update pival ##!!== OK 
  if(ppar$update["pival"]) {
    if(ppar$verbose) cat("pival:")
    state <- update.pival.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }

  for(j in c(1,2)) {  # j=1:DFR component; j=2:IFR component
    
    tparnm <- c(eta="eta",
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
                f1="f1", f2="f2")
    ntp <- names(tparnm)
    tparnm <- paste0(tparnm,j)
    names(tparnm) <- ntp

    # update eta ##!!== OK 
    eta.name <- paste0("eta",j)
    if(j==2 && ppar$update[eta.name]) { # only update eta2 (the IFR component)
      if(ppar$verbose) cat(eta.name)
      state <- update.eta.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update gamma ##!!== OK 
    gamma.name <- paste0("gamma",j)
    if(ppar$update[gamma.name]) {
      if(ppar$verbose) cat(gamma.name)
      state <- update.gamma.v2(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update thetavec ##!!== OK  
    thetavec.name <- paste0("thetavec",j)
    if(ppar$update[thetavec.name]) {
      if(ppar$verbose) cat(thetavec.name)
      state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update vvec ##!!== OK
    vvec.name <- paste0("vvec",j)
    if(ppar$update[vvec.name]) {
      if(ppar$verbose) cat(vvec.name)
      state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update alpha ##!!== OK  
    alpha.name <- paste0("alpha",j)
    if(ppar$update[alpha.name]) {
      if(ppar$verbose) cat(alpha.name)
      state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update beta ##!!== OK
    beta.name <- paste0("beta",j)
    if(ppar$update[beta.name]) {
      if(ppar$verbose) cat(beta.name)
      state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update phi ##!!==  OK
    phi.name <- paste0("phi",j)
    if(ppar$update[phi.name]) {
      if(ppar$verbose) cat(phi.name)
      state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update wvec ##!!== OK
    wvec.name <- paste0("wvec",j)
    if(ppar$update[wvec.name]) {
      if(ppar$verbose) cat(wvec.name)
      state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # swap elements of thetavec ##!!== OK
    thetaswap.name <- paste0("thetaswap",j)
    if(ppar$update[thetaswap.name]) {
      if(ppar$verbose) cat(thetaswap.name)
      state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
  }
  
  return(state)
}

update_state.lcv <- function(state, datlist, fpar, ppar, model) { 
  # Update a LCV state
  state$count <- state$count + 1
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update lambda0 ##!!== 
  if(ppar$update["lambda0"]) {
    if(ppar$verbose) cat("lambda0:")
    state <- update.lambda0.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update w0 ##!!== 
  if(ppar$update["w0"]) {
    if(ppar$verbose) cat("w0:")
    state <- update.w0.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update gamma ##!!== 
  if(ppar$update["gamma"]) {
    if(ppar$verbose) cat("gamma:")
    state <- update.gamma.v2(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update thetavec ##!!== OK
  if(ppar$update["thetavec"]) {
    if(ppar$verbose) cat("thetavec:")
    state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update vvec ##!!== OK
  if(ppar$update["vvec"]) {
    if(ppar$verbose) cat("vvec:")
    state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update alpha ##!!== OK
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update beta ##!!== OK
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update phi ##!!== OK
  if(ppar$update["phi"]) {
    if(ppar$verbose) cat("phi:")
    state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update wvec ##!!== OK
  if(ppar$update["wvec"]) {
    if(ppar$verbose) cat("wvec:")
    state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # swap elements of thetavec ##!!== OK  
  if(ppar$update["thetaswap"]) {
    if(ppar$verbose) cat("thetaswap:")
    state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }

  return(state)
}

update_state.circdr <- function(state, datlist, fpar, ppar, model) { 
  # Update an CIR or CDR state
  state$count <- state$count + 1
  
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update eta ##!!== OK 
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    state <- update.eta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update gamma ##!!== OK 
  if(ppar$update["gamma"]) {
    if(ppar$verbose) cat("gamma:")
    state <- update.gamma.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update thetavec ##!!== OK  
  if(ppar$update["thetavec"]) {
    if(ppar$verbose) cat("thetavec:")
    state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update vvec ##!!== OK
  if(ppar$update["vvec"]) {
    if(ppar$verbose) cat("vvec:")
    state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update alpha ##!!== OK  
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update beta ##!!== OK
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update phi ##!!==  OK
  if(ppar$update["phi"]) {
    if(ppar$verbose) cat("phi:")
    state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update wvec ##!!== OK
  if(ppar$update["wvec"]) {
    if(ppar$verbose) cat("wvec:")
    state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # swap elements of thetavec ##!!== OK  
  if(ppar$update["thetaswap"]) {
    if(ppar$verbose) cat("thetaswap:")
    state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  return(state)
}

update_state.cvx <- function(state, datlist, fpar, ppar, model) { 
  # Update a CVX state
  state$count <- state$count + 1
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update eta ##!!== OK 
  if(ppar$update["eta"]) {
    if(ppar$verbose) cat("eta:")
    tparnm <- parnm
    tparnm["gamma"] <- "gamma2"
    state <- update.eta.v1(state, datlist, fpar, ppar, model, nm=tparnm)
    if(ppar$verbose) cat("\n")
  }
  
  for(j in c(1,2)) {  # j=1:CDR component; j=2:CIR component
    
    tparnm <- c(gamma="gamma",
                thetavec="thetavec",
                vvec="vvec",
                alpha="alpha",
                beta="beta",
                phi="phi",
                uvec="uvec",
                wvec="wvec",
                thetaswap="thetaswap",
                f1="f1", f2="f2")
    ntp <- names(tparnm)
    tparnm <- paste0(tparnm,j)
    names(tparnm) <- ntp
    tparnm <- c(eta="eta",lambda0=ifelse(j==1,NA,"lambda0"),tparnm)
    
    # update gamma ##!!== OK 
    gamma.name <- paste0("gamma",j)
    if(ppar$update[gamma.name]) {
      if(ppar$verbose) cat(gamma.name)
      state <- update.gamma.v2(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update thetavec ##!!== OK  
    thetavec.name <- paste0("thetavec",j)
    if(ppar$update[thetavec.name]) {
      if(ppar$verbose) cat(thetavec.name)
      state <- update.thetavec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update vvec ##!!== OK
    vvec.name <- paste0("vvec",j)
    if(ppar$update[vvec.name]) {
      if(ppar$verbose) cat(vvec.name)
      state <- update.vvec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update alpha ##!!== OK  
    alpha.name <- paste0("alpha",j)
    if(ppar$update[alpha.name]) {
      if(ppar$verbose) cat(alpha.name)
      state <- update.alpha.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update beta ##!!== OK
    beta.name <- paste0("beta",j)
    if(ppar$update[beta.name]) {
      if(ppar$verbose) cat(beta.name)
      state <- update.beta.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update phi ##!!==  OK
    phi.name <- paste0("phi",j)
    if(ppar$update[phi.name]) {
      if(ppar$verbose) cat(phi.name)
      state <- update.phi.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # update wvec ##!!== OK
    wvec.name <- paste0("wvec",j)
    if(ppar$update[wvec.name]) {
      if(ppar$verbose) cat(wvec.name)
      state <- update.wvec.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
    
    # swap elements of thetavec ##!!== OK
    thetaswap.name <- paste0("thetaswap",j)
    if(ppar$update[thetaswap.name]) {
      if(ppar$verbose) cat(thetaswap.name)
      state <- update.thetaswap.v1(state, datlist, fpar, ppar, model, nm=tparnm)
      if(ppar$verbose) cat("\n")
    }
  }
  
  return(state)
}

update_state.mew <- function(state, datlist, fpar, ppar, model) { 
  # Update a MEW state
  state$count <- state$count + 1
  parnm <- unique(c(names(state),names(ppar$update)))
  names(parnm) <- parnm
  parnm <- c(parnm,c(f1="f1",f2="f2"))
  
  # update alpha ##!!== 
  if(ppar$update["alpha"]) {
    if(ppar$verbose) cat("alpha:")
    state <- update.alpha.mew(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update beta ##!!== 
  if(ppar$update["beta"]) {
    if(ppar$verbose) cat("beta:")
    state <- update.beta.mew(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update mu ##!!== 
  if(ppar$update["mu"]) {
    if(ppar$verbose) cat("mu:")
    state <- update.mu.mew(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  # update nu ##!!== 
  if(ppar$update["nu"]) {
    if(ppar$verbose) cat("nu:")
    state <- update.nu.mew(state, datlist, fpar, ppar, model, nm=parnm)
    if(ppar$verbose) cat("\n")
  }
  
  return(state)
}

####################################################################################
